
!##########################################################################
!###########################################################################
subroutine readinput
use velocity
use phase
use fastp
use param
implicit none

open(unit=55,file='input.dat',form='formatted',status='old')
!Time step parameters
read(55,*) restart
read(55,*) tstart
read(55,*) tfin
read(55,*) dump
!Flow parameters
read(55,*) inflow
read(55,*) dt
read(55,*) mu
read(55,*) rho
! forcing parameters
read(55,*) f1
read(55,*) f2
read(55,*) f3
read(55,*) k0
! phase-field parameters
read(55,*) radius
read(55,*) sigma
read(55,*) eps   


! compute pre-defined constant 
pi=4.d0*datan(1.d0)
lx=2.d0*pi
dx=lx/(dble(nx)-1)
dxi=1.d0/dx
ddxi=1.d0/dx/dx
rhoi=1.d0/rho

!for debug only, check parameters
write(*,*) "******MHIT36******"
write(*,*) "Restart ", restart
write(*,*) "Tstart  ", tstart
write(*,*) "Tfin    ", tfin
write(*,*) "Dump    ", dump
write(*,*) "inflow", inflow
write(*,*) "dt", dt
write(*,*) "mu", mu
write(*,*) "rho", rho
write(*,*) "f1,f2,f3,k0", f1,f2,f3,k0
write(*,*) "radius", radius
write(*,*) "sigma", sigma
write(*,*) "eps", eps
write(*,*) "lx", lx
write(*,*) "dx", dx
write(*,*) "dxi", dxi
write(*,*) "ddxi", ddxi
write(*,*) "rhoi", rhoi


end subroutine


