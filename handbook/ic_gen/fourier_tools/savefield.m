function savefield(U,V,W,N,fname,jnumb)

% save velocity fields to a binary as ${fname}_{U,V,W}.bin. If jnumb is set
% to 1, the Jane numbering transformation is performed

% optionally perform the Wybe to Jane numbering transformation (TODO Check)
if (jnumb == 1)
    U = U([N,1:N-1],:,:);
    V = V(:,[N,1:N-1],:);
    W = W(:,:,[N,1:N-1]);
end

% save the u field
floc = [fname,'_U.bin'];
fid = fopen(floc,'w');
fwrite(fid,U,'double');
fclose(fid);

% save the v field
floc = [fname,'_V.bin'];
fid = fopen(floc,'w');
fwrite(fid,V,'double');
fclose(fid);

% save the w field
floc = [fname,'_W.bin'];
fid = fopen(floc,'w');
fwrite(fid,W,'double');
fclose(fid);

end