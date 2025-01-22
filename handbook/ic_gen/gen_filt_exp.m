% makes spectra of the (filtered) experimental data

format long

addpath 'fourier_tools'

% do not use Jane numbering
jnumb = 0;

% set grid
N = 256;
L = 1;
alph = 4;

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

% make a spectrum
res = 10;
[k_mag_42,e_k_mag_42] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_42,e_k_mag_42,'CBC_exp_42');

% filter the field
[U_hat,V_hat,W_hat] = fftboxfilter(U_hat,V_hat,W_hat,N,L,k,alph);

% make a spectrum
res = 10;
[k_mag_42,e_k_mag_42] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_42,e_k_mag_42,'CBC_exp_fil_42');


% set data to fit
k_fit = k_98;
e_k_fit = E_98;

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% generate initial field
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit,e_k_fit);

% make a spectrum
res = 10;
[k_mag_98,e_k_mag_98] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_98,e_k_mag_98,'CBC_exp_98');

% filter the field
[U_hat,V_hat,W_hat] = fftboxfilter(U_hat,V_hat,W_hat,N,L,k,alph);

% make a spectrum
res = 10;
[k_mag_98,e_k_mag_98] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_98,e_k_mag_98,'CBC_exp_fil_98');

% set data to fit
k_fit = k_171;
e_k_fit = E_171;

% make the grid
[n,m,x,k] = makefftgrid(N,L);

% generate initial field
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit,e_k_fit);

% make a spectrum
res = 10;
[k_mag_171,e_k_mag_171] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_171,e_k_mag_171,'CBC_exp_171');

% filter the field
[U_hat,V_hat,W_hat] = fftboxfilter(U_hat,V_hat,W_hat,N,L,k,alph);

% make a spectrum
res = 10;
[k_mag_171,e_k_mag_171] = makespectrum(U_hat,V_hat,W_hat,N,L,m,res);

% save the spectrum
savespec(k_mag_171,e_k_mag_171,'CBC_exp_fil_171');


loglog(k_mag_42,e_k_mag_42,'r',k_mag_98,e_k_mag_98,'r',k_mag_171,e_k_mag_171,'r',k_42,E_42,'rs',k_98,E_98,'rs',k_171,E_171,'rs')