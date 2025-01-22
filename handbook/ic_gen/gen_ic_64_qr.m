% makes a 64^3 LES initial condition for the qr model

format long

addpath 'fourier_tools'

% do not use Jane's numbering
jnumb = 0;

% set grid
N = 64;
L = 1;

% load the filtered energy spectrum
[k_fit,e_k_fit] = loadspec('CBC_512');

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% generate initial field (do not use the k_fit including 0)
lkfit = length(k_fit);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));

% put this solution on a staggered grid (Wybe indexing)
[U_hat,V_hat,W_hat] = fftcoltostag(U_hat,V_hat,W_hat,N,m);

% transform to physical space
[U,V,W] = makeifft(U_hat,V_hat,W_hat);

% give the physical and spectral energy
enphys = enifft(U,V,W,N,L)
enspec = enfft(U_hat,V_hat,W_hat,N,L)

% check if real-valued
imagU = max(max(max(abs(imag(U)))))
imagV = max(max(max(abs(imag(V)))))
imagW = max(max(max(abs(imag(W)))))

% check the divergence-free condition
divmax = ifftdivmax(U,V,W,N,L)

% save the field to a file
savefield(U,V,W,N,'qr/CBC_64_qr',jnumb)

%transform back to cell-centered field
[U_hat,V_hat,W_hat] = fftstagtocol(U_hat,V_hat,W_hat,N,m);

% make a spectrum
res = 10;
[k_mag,e_k_mag] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% give the energy of the spectrum
e_sum = sum(e_k_mag)*2*pi/L

% save the spectrum
savespec(k_mag,e_k_mag,'qr/CBC_64_qr');

% draw spectra
loglog(k_mag,e_k_mag,'r',k_fit,e_k_fit,'b')

