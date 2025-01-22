function [U_c,V_c,W_c] = ifftgetstagfilt(U,V,W,N,alph)

% extracts the staggered values for a alph times coarser grid
N_c = N/alph;

U_c = zeros(N_c,N_c,N_c);
V_c = zeros(N_c,N_c,N_c);
W_c = zeros(N_c,N_c,N_c);

% Apply the box filter to each component
for in = 1:N_c
    for jn = 1:N_c
        for kn = 1:N_c
            U_c(in,jn,kn) = U(alph*in,alph/2+(jn-1)*alph,alph/2+(kn-1)*alph);
            V_c(in,jn,kn) = V(alph/2+(in-1)*alph,alph*jn,alph/2+(kn-1)*alph);
            W_c(in,jn,kn) = W(alph/2+(in-1)*alph,alph/2+(jn-1)*alph,alph*kn);
        end
    end
end

end