function [U_hat,V_hat,W_hat] = kangrescale(U_hat,V_hat,W_hat,N,L,m,k,k_fit,e_k_fit)

% rescale the modes by the ratio of their energy and the desired energy, do
% not make divergence free, the field is assumed to be already

dx = L/N;
dk = 2*pi/L;

% set the Fourier modes in one half of the plane
% the - N/2 is not set, because it is not possible to combine with the N/2
% mode (which is not captured) to get a real solution.
for im = 2:N/2+1
    for jm = 2:N
        for km = 2:N
            % determine the wave number of this m the desired energy (note
            % that 2*pi*r^2 is the area of a sphere, the (N/L)^3 takes
            % continuous transforms to the matlab FFTN
            k_mag = sqrt( k(im)^2 + k(jm)^2 + k(km)^2 );
            % the dx takes f_hat(k) to FFTN(k), the dk takes energy per k
            % to energy per m, and the area in m space instead of k space
            U_hat_mag_des = 1/dx^3*sqrt( 2*energy(k_mag,k_fit,e_k_fit)*dk/(4*pi*(k_mag/dk)^2) );
            
            % compute the actual magnitude of the mode
            U_hat_mag = sqrt(abs(U_hat(im,jm,km))^2 + abs(V_hat(im,jm,km))^2 + abs(W_hat(im,jm,km))^2);
            
            % rescale the field
            rratio = U_hat_mag_des/U_hat_mag;
            U_hat(im,jm,km) = rratio*U_hat(im,jm,km);
            V_hat(im,jm,km) = rratio*V_hat(im,jm,km);
            W_hat(im,jm,km) = rratio*W_hat(im,jm,km);
            
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