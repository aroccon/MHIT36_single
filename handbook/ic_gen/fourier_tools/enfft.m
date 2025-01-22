function en = enfft(U_hat,V_hat,W_hat,N,L)

% Compute the total energy of a transformed field

dx = L/N;

en = 0;
for im = 1:N
    for jm = 1:N
        for km = 1:N
            en = en + dx^6/L^3*(abs(U_hat(im,jm,km))^2 + abs(V_hat(im,jm,km))^2 + abs(W_hat(im,jm,km))^2)/2;
        end     
    end
end

end