function [k_mag,e_k_mag] = loadspec(fname)

% load energy spectrum from ${fname}_spec.mat. 

floc = [fname,'_spec.mat'];

load(floc);

end