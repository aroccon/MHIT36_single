program read_to_paraview
use commondata
implicit none

integer :: i,j,k
logical :: check
character(len=40) :: namefile

!! overide (if there are problems)
nstart=0
dump=100
nend=1700
nx=64

ny=nx
nz=nx

! read particles data
do i=nstart,nend,dump
 call read_particles(i)
enddo


end program read_to_paraview
