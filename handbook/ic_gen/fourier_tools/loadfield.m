function [U,V,W] = loadfield(N,fname,jnumb)

% load staggered velocity fields from a binary from the files 
% ${fname}_{U,V,W}.bin. If jnumb is set to 1, then the Jane numbering 
% transformation is performed.

% load the u field
floc = [fname,'_U.bin'];
fid = fopen(floc,'r');
U = reshape(fread(fid,N*N*N,'double'),N,N,N);
fclose(fid);

% load the v field
floc = [fname,'_V.bin'];
fid = fopen(floc,'r');
V = reshape(fread(fid,N*N*N,'double'),N,N,N);
fclose(fid);

% load the w field
floc = [fname,'_W.bin'];
fid = fopen(floc,'r');
W = reshape(fread(fid,N*N*N,'double'),N,N,N);
fclose(fid);

% optionally perform the Jane to Wybe numbering transformation (TODO Check)
if (jnumb == 1)
    U = U([2:N,1],:,:);
    V = V(:,[2:N,1],:);
    W = W(:,:,[2:N,1]);
end

end