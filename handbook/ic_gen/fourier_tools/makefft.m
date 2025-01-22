function [U_hat,V_hat,W_hat] = makefft(U,V,W)

% makes the Fourier transform and shift the wave numbers

% perform the Fourier transform
U_hat = fftn(U);
V_hat = fftn(V);
W_hat = fftn(W);

% shift the Fourier transform
U_hat = fftshift(U_hat);
V_hat = fftshift(V_hat);
W_hat = fftshift(W_hat);

end