function en = enifft(U,V,W,N,L)

% Compute the total energy of a field

dx = L/N;

en = 0;
for in = 1:N
    for jn = 1:N
        for kn = 1:N
            en = en + dx^3*(U(in,jn,kn)^2 + V(in,jn,kn)^2 + W(in,jn,kn)^2)/2;
        end     
    end
end

end