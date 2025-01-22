function [U,V,W] = makeifft(U_hat,V_hat,W_hat)

% Perform the inverse Fourier transform

% shift the Fourier transform
U_hat = ifftshift(U_hat);
V_hat = ifftshift(V_hat);
W_hat = ifftshift(W_hat);

% perform the Fourier transform
U = ifftn(U_hat);
V = ifftn(V_hat);
W = ifftn(W_hat);

end