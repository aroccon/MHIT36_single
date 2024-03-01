subroutine generate_output(nstep)

use commondata

integer :: i,k,j
character(len=40) :: namefile
character(len=10) :: lf
character(len=16) :: str4
character(len=80) :: buffer
character(len=3) :: setnum
integer :: nstep

! end of line character
lf=achar(10)

write(str4(1:16),'(i16)') nptot
!write(setnum,'(i3.3)') 1
write(namefile,'(a,i8.8,a)') './output/part_',nstep,'.vtk'
write(*,*) namefile

open(66,file=trim(namefile),status='new',form='unformatted',access='stream',convert='big_endian')

! start writing vtk file

! write header
buffer='# vtk DataFile Version 3.0'//lf
write(66) trim(buffer)
buffer='Particle data'//lf
write(66) trim(buffer)
buffer='BINARY'//lf
write(66) trim(buffer)
buffer='DATASET UNSTRUCTURED_GRID'//lf
write(66) trim(buffer)
buffer='POINTS '//str4//' float'//lf
write(66) trim(buffer)
do i=1,nptot
 write(66) real(pp(i,1)),real(pp(i,2)),real(pp(i,3))
enddo

buffer='CELL_TYPES'//str4//lf
write(66) trim(buffer)
do i=1,nptot
 write(66) 1 ! cell type 1, VTK_VERTEX
enddo

buffer='POINT_DATA '//str4//lf
write(66) trim(buffer)
buffer='SCALARS position integer 1'//lf
write(66) trim(buffer)
buffer='LOOKUP_TABLE default'//lf
write(66) trim(buffer)
do i=1,nptot
 write(66) 1 ! empty buffer
enddo


buffer=lf
write(66) trim(buffer)
close(66,status='keep')

return
end
