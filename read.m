clear all
nx=128;
fid = fopen('u_00001000.dat');
u = fread(fid, nx*nx*nx, 'double');
fclose(fid);
u=reshape(u,[nx nx nx]);

fid = fopen('v_00001000.dat');
v = fread(fid, nx*nx*nx, 'double');
fclose(fid);
v=reshape(v,[nx nx nx]);

fid = fopen('w_00001000.dat');
w = fread(fid, nx*nx*nx, 'double');
fclose(fid);
u=reshape(w,[nx nx nx]);

vtkwrite('u.vtk', 'structured_points', 'u',u)
vtkwrite('v.vtk', 'structured_points', 'v',v)
vtkwrite('w.vtk', 'structured_points', 'w',w)
