function [U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit,e_k_fit)

% generate a divergence-free collocated initial field with a desired energy
% spectrum

dx = L/N;
dk = 2*pi/L;

U_hat = zeros(N,N,N);
V_hat = zeros(N,N,N);
W_hat = zeros(N,N,N);

% generate the modified wave-number of the finite-difference approximation
% of the first derivative with respect to which the field should be
% divergence-free
kmod = 2i*sin(pi*m/N)/dx; % = 1i*sin(k*dx/2)/(dx/2)
% NOTE: for a cell-centered method this should be different


% set the Fourier modes in one half of the plane
% the - N/2 is not set, because it is not possible to combine with the N/2
% mode (which is not captured) to get a real solution.
for im = 2:N/2+1
    for jm = 2:N
        for km = 2:N
            % determine the wave number of this m and set the energy (note
            % that 2*pi*r^2 is the area of a sphere, the (N/L)^3 takes
            % continuous transforms to the matlab FFTN
            k_mag = sqrt( k(im)^2 + k(jm)^2 + k(km)^2 );
            % the dx takes f_hat(k) to FFTN(k), the dk takes energy per k
            % to energy per m, and the area in m space instead of k space
            U_hat_mag = 1/dx^3*sqrt( 2*energy(k_mag,k_fit,e_k_fit)*dk/(4*pi*(k_mag/dk)^2) );
            
            % set the wave number for the divergence-free condition
            kdiv = [kmod(im),kmod(jm),kmod(km)];
            
            % generate the random angles
            theta_rand = 2*pi*rand(1,3);
            
            % generate two unit basis vectors normal to kdiv and some
            % "random" vector
            u1 = cross([1.123,-0.982,1.182],kdiv);
            u1 = u1/norm(u1);
            u2 = cross(u1,kdiv)/norm(kdiv);
                        
            % randomize the real and imaginary parts of u_hat_mag
            U_hat_re = U_hat_mag*cos(theta_rand(1));
            U_hat_im = U_hat_mag*sin(theta_rand(1));
            
            % generate two different unit vectors normal to kdiv
            r1 = cos(theta_rand(2))*u1 + sin(theta_rand(2))*u2;
            r2 = cos(theta_rand(3))*u1 + sin(theta_rand(3))*u2;
            
            % generate the field
            U_hat(im,jm,km) = r1(1)*U_hat_re + 1i*r2(1)*U_hat_im;
            V_hat(im,jm,km) = r1(2)*U_hat_re + 1i*r2(2)*U_hat_im;
            W_hat(im,jm,km) = r1(3)*U_hat_re + 1i*r2(3)*U_hat_im;
        end
    end
end

% make the field real-valued
for im = N/2+2:N
    for jm = 2:N
        for km = 2:N
            U_hat(im,jm,km) = U_hat(N+2-im,N+2-jm,N+2-km)';
            V_hat(im,jm,km) = V_hat(N+2-im,N+2-jm,N+2-km)';
            W_hat(im,jm,km) = W_hat(N+2-im,N+2-jm,N+2-km)';
        end
    end
end

im = N/2+1;
for jm = N/2+2:N
    for km = 2:N
        U_hat(im,jm,km) = U_hat(N+2-im,N+2-jm,N+2-km)';
        V_hat(im,jm,km) = V_hat(N+2-im,N+2-jm,N+2-km)';
        W_hat(im,jm,km) = W_hat(N+2-im,N+2-jm,N+2-km)';
    end
end

im = N/2+1;
jm = N/2+1;
for km = N/2+2:N
    U_hat(im,jm,km) = U_hat(N+2-im,N+2-jm,N+2-km)';
    V_hat(im,jm,km) = V_hat(N+2-im,N+2-jm,N+2-km)';
    W_hat(im,jm,km) = W_hat(N+2-im,N+2-jm,N+2-km)';
end

% set the constant mode to 0
im = N/2+1;
jm = N/2+1;
km = N/2+1;
U_hat(im,jm,km) = 0;
V_hat(im,jm,km) = 0;
W_hat(im,jm,km) = 0;

end