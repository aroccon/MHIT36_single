function savespec(k_mag,e_k_mag,fname)

% save energy spectrum to a binary as ${fname}_spec.mat. 

floc = [fname,'_spec.mat'];

save(floc,'k_mag','e_k_mag');

end