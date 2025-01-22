function [U_hat,V_hat,W_hat] = fftcenttovert(U_hat,V_hat,W_hat,N,m)

% perform a cell center to vertex transformation in Fourier space
phfac = exp(1i*pi*m/N);

% Shift the phase of the Fourier transform so that evaluation gives the
% value at halve a grid cell later
for im = 1:N
    for jm = 1:N
        for km = 1:N
            trf = phfac(im)*phfac(jm)*phfac(km);
            U_hat(im,jm,km) = trf*U_hat(im,jm,km);
            V_hat(im,jm,km) = trf*V_hat(im,jm,km);
            W_hat(im,jm,km) = trf*W_hat(im,jm,km);
        end
    end
end

end