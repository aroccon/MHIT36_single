clear
nx=128;
np=8;
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

fid = fopen('xp_00000000.dat');
xp = fread(fid, np*3, 'double');
fclose(fid);
xp=reshape(xp,[3 np]);

vtkwrite('u.vtk', 'structured_points', 'u',u)
vtkwrite('v.vtk', 'structured_points', 'v',v)
vtkwrite('w.vtk', 'structured_points', 'w',w)
vtkwrite('w.vtk', 'structured_points', 'w',w)

