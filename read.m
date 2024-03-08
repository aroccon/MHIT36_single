clear
nx=256;
np=0;
fid = fopen('u_00020000.dat');
u = fread(fid, nx*nx*nx, 'double');
fclose(fid);
u=reshape(u,[nx nx nx]);

fid = fopen('v_00020000.dat');
v = fread(fid, nx*nx*nx, 'double');
fclose(fid);
v=reshape(v,[nx nx nx]);

fid = fopen('w_00020000.dat');
w = fread(fid, nx*nx*nx, 'double');
fclose(fid);
w=reshape(w,[nx nx nx]);

fid = fopen('phi_00000000.dat');
phi = fread(fid, np*nx*nx, 'double');
fclose(fid);
phi=reshape(phi,[nx nx nx]);
vtkwrite('phi.vtk', 'structured_points', 'phi',phi)


vtkwrite('u.vtk', 'structured_points', 'u',u)
vtkwrite('v.vtk', 'structured_points', 'v',v)
vtkwrite('w.vtk', 'structured_points', 'w',w)
vtkwrite('w.vtk', 'structured_points', 'w',w)

