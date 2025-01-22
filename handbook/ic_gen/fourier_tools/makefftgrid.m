function [n,m,x,k] = makefftgrid(N,L)

% make the grids in physical and Fourier space to connect the discrete
% Fourier transform to its physical meaning

dx = L/N;
dk = 2*pi/L;

% make the index vectors
n = [0:N-1];
m = [-N/2:N/2-1];
x = n*dx;
k = m*dk;

end