function [U_hat,V_hat,W_hat] = fftboxfilter(U_hat,V_hat,W_hat,N,L,k,alph)

dx = L/N;

% perform a box filter operation is spectral space
ffac = sin(alph*dx*k/2)./(alph*dx*k/2);
ffac(N/2+1) = 1; % this gives a NaN otherwise

% Apply the box filter to each component
for im = 1:N
    for jm = 1:N
        for km = 1:N
            boxfac = ffac(im)*ffac(jm)*ffac(km);
            U_hat(im,jm,km) = boxfac*U_hat(im,jm,km);
            V_hat(im,jm,km) = boxfac*V_hat(im,jm,km);
            W_hat(im,jm,km) = boxfac*W_hat(im,jm,km);
        end
    end
end

end