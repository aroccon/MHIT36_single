% makes a 512^3 initial condition for an (optional) dns and compute the
% spectra of a box-filtered field for to generate the initial condition for
% the dynamic Smagorinsky model

format long

addpath 'fourier_tools'

% do not use Jane numbering
jnumb = 0;

% set grid
N = 512;
L = 1;
alph = 8;

% load experimental data
load('CBC_exp.mat')
% Nondimensinalization
M = 5.08; % in cm
U0 = 1000; % in cm/s
L_ref = 11*M; % in cm
u_ref = sqrt(3/2)*22.2; % in cm/s
k_42 = k_42 * L_ref;
E_42 = E_42 / (u_ref^2*L_ref);
k_98 = k_98 * L_ref;
E_98 = E_98 / (u_ref^2*L_ref);
k_171 = k_171 * L_ref;
E_171 = E_171 / (u_ref^2*L_ref);
% set data to fit
k_fit = k_42;
e_k_fit = E_42;

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% generate initial field
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit,e_k_fit);

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
savefield(U,V,W,N,'CBC_512',jnumb)

%transform back to cell-centered field
[U_hat,V_hat,W_hat] = fftstagtocol(U_hat,V_hat,W_hat,N,m);

% make a spectrum
res = 10;
[k_mag,e_k_mag] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% give the energy of th spectrum
e_sum = sum(e_k_mag)*2*pi/L

% save the spectrum
savespec(k_mag,e_k_mag,'CBC_512');

% filter the Fourier transform of the initial condition
[U_hat,V_hat,W_hat] = fftboxfilter(U_hat,V_hat,W_hat,N,L,k,alph);

% make a spectrum
res = 10;
[k_mag_f,e_k_mag_f] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_f,e_k_mag_f,'CBC_512_filt');

% draw spectra
loglog(k_mag,e_k_mag,'r',k_fit,e_k_fit,'b',k_mag_f,e_k_mag_f,'g')

