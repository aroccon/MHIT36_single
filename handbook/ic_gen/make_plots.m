% draws the spectra of the qr and the dsm method

close all

addpath 'fourier_tools'

% load the experimental data points
load('CBC_exp.mat')
M = 5.08; % in cm
U0 = 1000; % in cm/s
L_ref = 11*M; % in cm
u_ref = sqrt(3/2)*22.2; % in cm/s
k_exp_42 = k_42 * L_ref;
e_k_exp_42 = E_42 / (u_ref^2*L_ref);
k_exp_98 = k_98 * L_ref;
e_k_exp_98 = E_98 / (u_ref^2*L_ref);
k_exp_171 = k_171 * L_ref;
e_k_exp_171 = E_171 / (u_ref^2*L_ref);

% load the fitted spectra
[k_fexp_42,e_k_fexp_42] = loadspec('CBC_exp_42');
[k_fexp_98,e_k_fexp_98] = loadspec('CBC_exp_98');
[k_fexp_171,e_k_fexp_171] = loadspec('CBC_exp_171');

% load the fitted filtered spectra
[k_ffexp_42,e_k_ffexp_42] = loadspec('CBC_exp_fil_42');
[k_ffexp_98,e_k_ffexp_98] = loadspec('CBC_exp_fil_98');
[k_ffexp_171,e_k_ffexp_171] = loadspec('CBC_exp_fil_171');

% load the spectra from the qr model
[k_qr_42,e_k_qr_42] = loadspec('qr/CBC_64_qr_42');
[k_qr_98,e_k_qr_98] = loadspec('qr/CBC_64_qr_98');
[k_qr_171,e_k_qr_171] = loadspec('qr/CBC_64_qr_171');

% load the spectra from the dsm model
[k_dsm_42,e_k_dsm_42] = loadspec('dsm/CBC_64_dsm_42');
[k_dsm_98,e_k_dsm_98] = loadspec('dsm/CBC_64_dsm_98');
[k_dsm_171,e_k_dsm_171] = loadspec('dsm/CBC_64_dsm_171');

% load the energy decay
[en_dsm] = load('dsm/bulk01.dat');
[en_qr] = load('qr/bulk01.dat');

% ----------------------------------------------------
% fit a 64^3 field to the (filtered) experimental data to compare energy

en = zeros(3,3);
en(:,1) = [42,98,171]';
N = 64;
L = 1;
[n,m,x,k] = makefftgrid(N,L);

% unfiltered
k_fit = k_fexp_42;
e_k_fit = e_k_fexp_42;
lkfit = length(k_fexp_42);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(1,2) = enfft(U_hat,V_hat,W_hat,N,L);

k_fit = k_fexp_98;
e_k_fit = e_k_fexp_98;
lkfit = length(k_fexp_98);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(2,2) = enfft(U_hat,V_hat,W_hat,N,L);

k_fit = k_fexp_171;
e_k_fit = e_k_fexp_171;
lkfit = length(k_fexp_171);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(3,2) = enfft(U_hat,V_hat,W_hat,N,L);

% filtered
k_fit = k_ffexp_42;
e_k_fit = e_k_ffexp_42;
lkfit = length(k_ffexp_42);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(1,3) = enfft(U_hat,V_hat,W_hat,N,L);

k_fit = k_ffexp_98;
e_k_fit = e_k_ffexp_98;
lkfit = length(k_ffexp_98);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(2,3) = enfft(U_hat,V_hat,W_hat,N,L);

k_fit = k_ffexp_171;
e_k_fit = e_k_ffexp_171;
lkfit = length(k_ffexp_171);
[U_hat,V_hat,W_hat] = makefield(N,L,m,k,k_fit(2:lkfit),e_k_fit(2:lkfit));
en(3,3) = enfft(U_hat,V_hat,W_hat,N,L);


% ----------------------------------------------------

% make the plots

% - spectra
figure(1)

loglog(k_dsm_42,e_k_dsm_42,'b',k_qr_42,e_k_qr_42,'r',k_fexp_42,e_k_fexp_42,'g',k_ffexp_42,e_k_ffexp_42,'c',k_exp_42,e_k_exp_42,'sk',k_dsm_98,e_k_dsm_98,'b',k_qr_98,e_k_qr_98,'r',k_fexp_98,e_k_fexp_98,'g',k_ffexp_98,e_k_ffexp_98,'c',k_exp_98,e_k_exp_98,'sk',k_dsm_171,e_k_dsm_171,'b',k_qr_171,e_k_qr_171,'r',k_fexp_171,e_k_fexp_171,'g',k_ffexp_171,e_k_ffexp_171,'c',k_exp_171,e_k_exp_171,'sk')
legend('dyn. Sm.','QR','fit. exp.','fit. filt. exp.','exp.','Location','SouthWest')
xlabel('k')
ylabel('E(k)')

xlim([2*pi,sqrt(3*32^2)*2*pi])
ylim([1e-5,2e-2])

title('Spectra at t = 42, 98, 171')



% - energy decay

figure(2)

% to express time in units from paper (see notes_nondim.pdf)
tfac = 1/(2.471818e-3);
toff = 0.103816;

plot(tfac*(en_dsm(:,2)+toff),en_dsm(:,3),'r',tfac*(en_qr(:,2)+toff),en_qr(:,3),'b',en(:,1),en(:,3),'sr',en(:,1),en(:,2),'sb')
legend('dsm','qr','fit. box-filt. exp.','fit. exp.')
xlabel('t U_0 / M')
ylabel('E_k')





title('energy decay')