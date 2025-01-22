function [U_hat,V_hat,W_hat] = fftcoltostag(U_hat,V_hat,W_hat,N,m)

% perform a collocated to staggered transformation in Fourier space
phfac = exp(1i*pi*m/N);

% Shift the phase of the Fourier transform so that evaluation gives the
% value at halve a grid cell later
for im = 1:N
    for jm = 1:N
        for km = 1:N
            U_hat(im,jm,km) = phfac(im)*U_hat(im,jm,km);
            V_hat(im,jm,km) = phfac(jm)*V_hat(im,jm,km);
            W_hat(im,jm,km) = phfac(km)*W_hat(im,jm,km);
        end
    end
end

end