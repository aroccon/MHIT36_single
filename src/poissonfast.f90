subroutine poissonfast
use openacc
use cuFFT
use iso_c_binding
use fastp
use param
implicit none
integer :: gerr,i,j,k
double precision :: pm

! A. Roccon 28/02/2024
! 3D FFT-based solution of p_xx + p_yy  + p_zz = rhs;
! with periodic boundary conditions along all directions


!Perform FFT3D forward of the rhsp
!$acc host_data use_device(rhsp,pc)
gerr = gerr + cufftExecD2Z(cudaplan_fwd,rhsp,pc)
!$acc end host_data

!$acc parallel loop collapse(3) present(pc,delsq)
do i=1,nx/2+1
  do j=1,nx
    do k=1,nx
        pc(i,j,k) = pc(i,j,k) / delsq(i,j,k)
    enddo
  enddo
enddo

!$acc parallel present(pc)
pc(1,1,1) = 0.d0
!$acc end parallel

!$acc host_data use_device(pc,p)
gerr = gerr + cufftExecZ2D(cudaplan_bwd,pc,p)
!$acc end host_data

! scale p
!$acc parallel loop present(p)
do k=1,nx
  do j=1,nx
    do i=1,nx
      p(i,j,k) = p(i,j,k) / (nx*nx*nx)
    enddo
  enddo
enddo


return
end subroutine





subroutine init_cufft
use openacc
use cuFFT
use iso_c_binding
use fastp
use param
implicit none
integer :: gerr,i,j,k
integer(kind=int_ptr_kind()) :: workSize(1)
double precision :: kx, ky, kz

! create plans
! Creata plans (forth and back)
!Plan forward
gerr=0
gerr=gerr+cufftCreate(cudaplan_fwd)
gerr=gerr+cufftMakePlan3d(cudaplan_fwd,nx,nx,nx,CUFFT_D2Z,workSize)
if (gerr.ne.0) write(*,*) "Error in cuFFT plan FWD:", gerr

!Plan backward
gerr=0
gerr=gerr+cufftCreate(cudaplan_bwd)
gerr=gerr+cufftMakePlan3d(cudaplan_bwd,nx,nx,nx,CUFFT_Z2D,workSize)
if (gerr.ne.0) write(*,*) "Error in cuFFT plan BWD:", gerr


!create wavenumber
! wavenumbers 
do i=1,nx/2
    kk(i)=i-1
enddo
do i=nx/2+1,nx
    kk(i)=-nx + i -1
enddo

do i=1,nx
    kx = 2.d0*pi*kk(i)/lx
    do j=1,nx
        ky = 2.d0*pi*kk(j)/lx
        do k=1,nx
            kz = 2.d0*pi*kk(k)/lx
            delsq(i,j,k) = ( 2.d0*(cos(kx*dx) - 1.d0)  + 2.d0*(cos(ky*dx) - 1.d0)    + 2.d0*(cos(kz*dx) - 1.d0) ) / dx**2
        enddo
    enddo
enddo




!create delsq
!do i=1,nx
!    do j=1,nx
!        do k=1,nx
!            delsq(i,j,k) = -(kk(i)**2d0 + kk(j)**2d0 + kk(k)**2d0)
!        enddo
!    enddo
!enddo

delsq(1,1,1) = 1.d0

end subroutine
