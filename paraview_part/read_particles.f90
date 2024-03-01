subroutine read_particles(nstep)

use commondata

integer :: nstep, idp
integer :: dummy
character(len=40) :: namedir,namefile
character(len=8) :: numfile
character(len=3) :: setnum
logical :: check



namedir='../src/output/'
write(numfile,'(i8.8)') nstep


namefile=trim(namedir)//'xp_'//numfile//'.dat'
write(*,*) namefile
inquire(file=trim(namefile),exist=check)

write(*,*) "Reading"

allocate(pp(nptot,3))

if(check.eqv..true.)then
write(*,*) 'Reading step ',nstep,' out of ',nend,' , particles'
  !reading particles position
  namefile=trim(namedir)//'xp_'//numfile//'.dat'
  open(666,file=trim(namefile),form='unformatted',access='stream',status='old',convert='little_endian')
  !write(*,*) "dummy", dummy !debug only
  do i=1,nptot
     !read id particle + position
     read(666) pp(i,1:3)
  enddo
  close(666,status='keep')


 ! write(*,*) pp(2,3)
  ! generate paraview output file
  call generate_output(nstep)
endif

deallocate(pp)

return
end
