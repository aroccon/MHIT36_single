function divmax = ifftdivmax(U,V,W,N,L)

% compute the maximum divergence

divmax = 0;
dx = L/N;
for in = 2:N
    for jn = 2:N
        for kn = 2:N
            divijk = (U(in,jn,kn) - U(in-1,jn,kn))/dx + (V(in,jn,kn) - V(in,jn-1,kn))/dx + (W(in,jn,kn) - W(in,jn,kn-1))/dx;
            if (abs(divijk) > divmax)
                divmax = abs(divijk);
            end
        end
    end
end

end