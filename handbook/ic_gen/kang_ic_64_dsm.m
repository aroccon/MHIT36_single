% Kang rescale the dynamic Smagorinsky model at t' = 42

format long

addpath 'fourier_tools'

% do not use Jane's numbering
jnumb = 0;

% set grid
N = 64;
L = 1;

% load the filtered 512 energy spectrum for Kang rescaling
[k_fit,e_k_fit] = loadspec('CBC_512_filt');

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% load the instantaneous conditions
[U,V,W] = loadfield(N,'dsm/CBC_64_dsm_42',jnumb);

enphys = enifft(U,V,W,N,L)

% and Fourier transform them
[U_hat,V_hat,W_hat] = makefft(U,V,W);

% take the loaded field back to cell-center data
[U_hat,V_hat,W_hat] = fftstagtocol(U_hat,V_hat,W_hat,N,m);

% Apply Kang rescaling
res = 10;
[k_mag_kang,e_k_mag_kang] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% Kang-rescale the Fourier modes
lkfit = length(k_fit);
[U_hat,V_hat,W_hat] = kangrescale(U_hat,V_hat,W_hat,N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));

% put the rescaled field on a staggered grid (Wybe indexing)
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
savefield(U,V,W,N,'dsm/CBC_64_dsm_42_kang',jnumb)

%transform back to cell-centered field
[U_hat,V_hat,W_hat] = fftstagtocol(U_hat,V_hat,W_hat,N,m);

% make a spectrum
res = 10;
[k_mag,e_k_mag] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% give the energy of th spectrum
e_sum = sum(e_k_mag)*2*pi/L

% save the spectrum 
savespec(k_mag,e_k_mag,'dsm/CBC_64_dsm_42_kang');

% draw spectra
loglog(k_mag,e_k_mag,'r',k_fit,e_k_fit,'b',k_mag_kang,e_k_mag_kang,'g')

