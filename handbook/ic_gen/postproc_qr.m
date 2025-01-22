% postprocess the qr data

format long

addpath 'fourier_tools'

% do not use Jane's numbering
jnumb = 0;

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

% load the filtered dns data
[k_fdns_42,e_k_fdns_42] = loadspec('CBC_exp_42');
[k_fdns_98,e_k_fdns_98] = loadspec('CBC_exp_98');
[k_fdns_171,e_k_fdns_171] = loadspec('CBC_exp_171');

% set grid
N = 64;
L = 1;

% load the energy spectrum
[k_fit,e_k_fit] = loadspec('CBC_512');

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% load the instantaneous conditions
[U1,V1,W1] = loadfield(N,'qr/CBC_64_qr_42_kang',jnumb);
[U2,V2,W2] = loadfield(N,'qr/CBC_64_qr_98',jnumb);
[U3,V3,W3] = loadfield(N,'qr/CBC_64_qr_171',jnumb);

enphys1 = enifft(U1,V1,W1,N,L)
enphys2 = enifft(U2,V2,W2,N,L)
enphys3 = enifft(U3,V3,W3,N,L)

% and Fourier transform them
[U1_hat,V1_hat,W1_hat] = makefft(U1,V1,W1);
[U2_hat,V2_hat,W2_hat] = makefft(U2,V2,W2);
[U3_hat,V3_hat,W3_hat] = makefft(U3,V3,W3);

% take the loaded field back to cell-center data
[U1_hat,V1_hat,W1_hat] = fftstagtocol(U1_hat,V1_hat,W1_hat,N,m);
[U2_hat,V2_hat,W2_hat] = fftstagtocol(U2_hat,V2_hat,W2_hat,N,m);
[U3_hat,V3_hat,W3_hat] = fftstagtocol(U3_hat,V3_hat,W3_hat,N,m);

% make a spectrum
res = 10;
[k_mag_qr1,e_k_mag_qr1] = makespectrum(U1_hat,V1_hat,W1_hat,N,L,m,res);
[k_mag_qr2,e_k_mag_qr2] = makespectrum(U2_hat,V2_hat,W2_hat,N,L,m,res);
[k_mag_qr3,e_k_mag_qr3] = makespectrum(U3_hat,V3_hat,W3_hat,N,L,m,res);

% save the spectrum (for use in Kang rescaling)
savespec(k_mag_qr1,e_k_mag_qr1,'qr/CBC_64_qr_42');
savespec(k_mag_qr2,e_k_mag_qr2,'qr/CBC_64_qr_98');
savespec(k_mag_qr3,e_k_mag_qr3,'qr/CBC_64_qr_171');

% draw spectra
loglog(k_fit,e_k_fit,'b',k_mag_qr1,e_k_mag_qr1,'r',k_mag_qr2,e_k_mag_qr2,'r',k_mag_qr3,e_k_mag_qr3,'r',k_42,E_42,'ks',k_98,E_98,'ks',k_171,E_171,'ks',k_fdns_42,e_k_fdns_42,'g',k_fdns_98,e_k_fdns_98,'g',k_fdns_171,e_k_fdns_171,'g')

