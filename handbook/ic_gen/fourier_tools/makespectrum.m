function [k_mag,e_k_mag] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res)

% generate a spectrum of the Fourier transforms with a given resolution

% the resolution of this procedure (n.o. points per |k|) is small for lower
% wave numbers and larger for high wave numbers. Therefore res should be
% quite large to not make considerable errors in the k = 1 energy.

% make the physical Fourier transforms
U_hat_phys = (L/N)^3*U_hat;
V_hat_phys = (L/N)^3*V_hat;
W_hat_phys = (L/N)^3*W_hat;

% make the refined wavenumber grid on [-1/2,1/2] (to compute |k| more 
% accurately)
dr = 1/res;
gridim = zeros(res,res,res);
gridjm = zeros(res,res,res);
gridkm = zeros(res,res,res);
for im = 1:res
    for jm = 1:res
        for km = 1:res
            gridim(im,jm,km) = (-1/2 + im*dr - dr/2);
            gridjm(im,jm,km) = (-1/2 + jm*dr - dr/2);
            gridkm(im,jm,km) = (-1/2 + km*dr - dr/2);
        end
    end
end

% compute the energy of the field
e_phys = 1/2*(U_hat_phys.*conj(U_hat_phys) + V_hat_phys.*conj(V_hat_phys) + W_hat_phys.*conj(W_hat_phys));

% make the energy spectrum
m_max = ceil(sqrt(3)*N/2);
e_k_mag = zeros(1,m_max+1);


for im = 1:N
    for jm = 1:N
        for km = 1:N
            % make a grid of |m| values around each m = (im,jm,km)
            mgrid = reshape(round(sqrt( (gridim + m(im)).^2 + (gridjm + m(jm)).^2 + (gridkm + m(km)).^2 )),res^3,1);
            % add the energy for each |m|
            e_k_per_cell = e_phys(im,jm,km)/res^3;
            for ii = 1:res^3
                e_k_mag(mgrid(ii)+1) = e_k_mag(mgrid(ii)+1) + e_k_per_cell;
            end
        end
    end
end

% compute the energy per physical wave number k (instead of m)
k_mag = 2*pi/L*[0:m_max];
e_k_mag = e_k_mag*L/(2*pi);

end