program mhit
! A. Roccon 08/02/2024
! Homogenous isotropic turbulence solver
! Constant density and viscosity
! 2nd order finite difference + fastPoisson3D solver for pressure 
! ABC forcing scheme (see Comparison of forcing schemes to sustain
! homogeneous isotropic turbulence)
! Must be validated, qualitative results seems fine
! Runs on Nvidia GPU (cuFFT), FFTW to be implemented 
! Eulero explicit in time (to be ubgraded to AB for NS)


use openacc
use fastp
use param

#define phiflag 0
#define partflag 0
#define openaccflag 1

implicit none
double precision :: lx,dx,dxi,ddxi
double precision :: dt
double precision :: rho,mu
double precision :: f1,f2,f3,k0
double precision :: tstart,tend
double precision :: uc,vc,wc
double precision :: h11,h12,h13,h21,h22,h23,h31,h32,h33,cou
integer :: tfin,i,j,k,t,im,jm,km,ip,jp,kp
double precision :: x(nx),pi
double precision, allocatable :: div(:,:,:)
double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
double precision, allocatable :: ustar(:,:,:), vstar(:,:,:), wstar(:,:,:)
double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
double precision, allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
!PFM variables
double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:)
double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
!particle variables
double precision, allocatable :: xp(:,:), vp(:,:), ufp(:,:), fp(:,:) 

call acc_set_device_num(1,acc_device_nvidia)

! initialize parameters
tfin=10
dt=0.001d0
pi=4.d0*datan(1.d0)
lx=2.d0*pi
dx=lx/(dble(nx)-1)
dxi=1.d0/dx
ddxi=1.d0/dx/dx
rho=1.d0
mu=0.03d0
! forcing parameters (ABC)
f1=1.d0
f2=1.d0
f3=1.d0
k0=2.d0


!allocate variables
!NS variables
allocate(u(nx,nx,nx),v(nx,nx,nx),w(nx,nx,nx))
allocate(ustar(nx,nx,nx),vstar(nx,nx,nx),wstar(nx,nx,nx))
allocate(rhsu(nx,nx,nx),rhsv(nx,nx,nx),rhsw(nx,nx,nx))
allocate(fx(nx,nx,nx),fy(nx,nx,nx),fz(nx,nx,nx))
allocate(div(nx,nx,nx))
allocate(delsq(nx,nx,nx))
allocate(kk(nx),kx(nx,nx,nx),ky(nx,nx,nx),kz(nx,nx,nx))
!PFM variables
#if phiflag == 1
allocate(phi(nx,nx,nx),rhsphi(nx,nx,nx))
allocate(normx(nx,nx,nx),normy(nx,nx,nx),normz(nx,nx,nx))
#endif
!particles arrays
#if partflag == 1
allocate(xp(np,3),vp(np,3),ufp(np,3),fp(np,3)
#endif

x(1)=0.0d0
do i=2,nx
    x(i)=x(i-1)+dx
enddo    

write(*,*) "Initialize velocity field"
!$acc kernels
do k = 1,nx
    do j= 1,nx
        do i = 1,nx
            u(i,j,k) =  sin(x(i))*cos(x(j))*cos(x(k))
            v(i,j,k) =  cos(x(i))*sin(x(j))*cos(x(k))
            w(i,j,k) =  0.d0
       enddo
    enddo
enddo
!$acc end kernels
!write(*,*) 'Initialize phase field (WIP)'

!initialize the plan for cuFFT
call init_cufft

!Save initial fields
!call write_output(t)

! Start temporal loop
do t=1,tfin

    call cpu_time(tstart)
    write(*,*) "Time step",t,"of",tfin

    ! Advance marker function
    ! to be implemented, second-order PFM, CAC

    ! Projection step, convective terms
    !Convective terms NS
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                ip=i+1
                jp=j+1
                kp=k+1
                im=i-1
                jm=j-1
                km=k-1
                if (ip .gt. nx) ip=1
                if (jp .gt. nx) jp=1
                if (kp .gt. nx) kp=1   
                if (im .lt. 1) im=nx
                if (jm .lt. 1) jm=nx
                if (km .lt. 1) km=nx 
                ! compute the products (conservative form)
                h11 = (u(ip,j,k)+u(i,j,k))*(u(ip,j,k)+u(i,j,k))     - (u(i,j,k)+u(im,j,k))*(u(i,j,k)+u(im,j,k))
                h12 = (u(i,jp,k)+u(i,j,k))*(v(i,jp,k)+v(im,jp,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
                h13 = (u(i,j,kp)+u(i,j,k))*(w(i,j,kp)+w(im,j,kp))   - (u(i,j,k)+u(i,j,km))*(w(i,j,k)+w(im,j,k))
                h21 = (u(ip,j,k)+u(ip,jm,k))*(v(ip,j,k)+v(i,j,k))   - (u(i,j,k)+u(i,jm,k))*(v(i,j,k)+v(im,j,k))
                h22 = (v(i,jp,k)+v(i,j,k))*(v(i,jp,k)+v(i,j,k))     - (v(i,j,k)+v(i,jm,k))*(v(i,j,k)+v(i,jm,k))
                h23 = (w(i,j,kp)+w(i,jm,kp))*(v(i,j,kp)+v(i,j,k))   - (w(i,j,k)+w(i,jm,k))*(v(i,j,k)+v(i,j,km))
                h31 = (w(ip,j,k)+w(i,j,k))*(u(ip,j,k)+u(ip,j,km))   - (w(i,j,k)+w(im,j,k))*(u(i,j,k)+u(i,j,km))
                h32 = (v(i,jp,k)+v(i,jp,km))*(w(i,jp,k)+w(i,j,k))   - (v(i,j,k)+v(i,j,km))*(w(i,j,k)+w(i,jm,k))
                h33 = (w(i,j,kp)+w(i,j,k))*(w(i,j,kp)+w(i,j,k))     - (w(i,j,k)+w(i,j,km))*(w(i,j,k)+w(i,j,km))
                ! compute the derivative
                h11=h11*0.25*dxi
                h12=h12*0.25*dxi
                h13=h13*0.25*dxi
                h21=h21*0.25*dxi
                h22=h22*0.25*dxi
                h23=h23*0.25*dxi
                h31=h31*0.25*dxi
                h32=h32*0.25*dxi
                h33=h33*0.25*dxi
                ! add to the rhs
                rhsu(i,j,k)=-(h11+h12+h13)
                rhsv(i,j,k)=-(h21+h22+h23)
                rhsw(i,j,k)=-(h31+h32+h33)
            enddo
        enddo
    enddo
    !$acc end kernels
  
    ! Compute viscous terms
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                ip=i+1
                jp=j+1
                kp=k+1
                im=i-1
                jm=j-1
                km=k-1
                if (ip .gt. nx) ip=1
                if (jp .gt. nx) jp=1
                if (kp .gt. nx) kp=1
                if (im .lt. 1) im=nx
                if (jm .lt. 1) jm=nx
                if (km .lt. 1) km=nx 
                h11 = mu*(u(ip,j,k)-2*u(i,j,k)+u(im,j,k))*ddxi
                h12 = mu*(u(i,jp,k)-2*u(i,j,k)+u(i,jm,k))*ddxi
                h13 = mu*(u(i,j,kp)-2*u(i,j,k)+u(i,j,km))*ddxi
                h21 = mu*(v(ip,j,k)-2*v(i,j,k)+v(im,j,k))*ddxi
                h22 = mu*(v(i,jp,k)-2*v(i,j,k)+v(i,jm,k))*ddxi
                h23 = mu*(v(i,j,kp)-2*v(i,j,k)+v(i,j,km))*ddxi
                h31 = mu*(w(ip,j,k)-2*w(i,j,k)+w(im,j,k))*ddxi
                h32 = mu*(w(i,jp,k)-2*w(i,j,k)+w(i,jm,k))*ddxi
                h33 = mu*(w(i,j,kp)-2*w(i,j,k)+w(i,j,km))*ddxi
                rhsu(i,j,k)= rhsu(i,j,k)+(h11+h12+h13)/rho
                rhsv(i,j,k)= rhsv(i,j,k)+(h21+h22+h23)/rho
                rhsw(i,j,k)= rhsw(i,j,k)+(h31+h32+h33)/rho
            enddo
        enddo
    enddo
    !$acc end kernels

    ! forcing term (always x because is the same axis)
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                fx(i,j,k)=f1*sin(k0*x(k))+f3*sin(k0*x(j))
                fy(i,j,k)=f2*sin(k0*x(i))+f1*sin(k0*x(k))
                fz(i,j,k)=f3*sin(k0*x(j))+f2*sin(k0*x(i))
                rhsu(i,j,k)= rhsu(i,j,k) + fx(i,j,k)
                rhsv(i,j,k)= rhsv(i,j,k) + fy(i,j,k)
                rhsw(i,j,k)= rhsw(i,j,k) + fz(i,j,k)
            enddo
        enddo
    enddo
    !$acc end kernels

    ! find u, v and w star (explicit Eulero)
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                ustar(i,j,k) = u(i,j,k) + dt*rhsu(i,j,k)
                vstar(i,j,k) = v(i,j,k) + dt*rhsv(i,j,k)
                wstar(i,j,k) = w(i,j,k) + dt*rhsw(i,j,k)
            enddo
        enddo
    enddo
    !$acc end kernels
    !write(*,*) "maxval u*", maxval(ustar)
    !write(*,*) "maxval v*", maxval(vstar)
    !write(*,*) "maxval w*", maxval(wstar)
    ! Compute rhs of Poisson equation div*ustar
    ! Compute divergence at the cell center (nx+1 avilable on u,v,w)
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                ip=i+1
                jp=j+1
                kp=k+1
                if (ip > nx) ip=1
                if (jp > nx) jp=1
                if (kp > nx) kp=1
                rhsp(i,j,k) = (rho*dxi/dt)*(ustar(ip,j,k)-ustar(i,j,k))
                rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(vstar(i,jp,k)-vstar(i,j,k))
                rhsp(i,j,k) = rhsp(i,j,k) + (rho*dxi/dt)*(wstar(i,j,kp)-wstar(i,j,k))
            enddo
        enddo
    enddo
    !$acc end kernels

    !write(*,*) "max rhsp", maxval(rhsp)

    ! call Poisson solver (3DFastPoissonsolver, periodic BCs)
    call poissonfast

    ! Correct velocity 
   !$acc kernels 
    do i=1,nx
        do j=1,nx
            do k=1,nx
                im=i-1
                jm=j-1
                km=k-1
                if (im < 1) im=nx
                if (jm < 1) jm=nx
                if (km < 1) km=nx   
                u(i,j,k)=ustar(i,j,k) - dt/rho*(p(i,j,k)-p(im,j,k))*dxi
                v(i,j,k)=vstar(i,j,k) - dt/rho*(p(i,j,k)-p(i,jm,k))*dxi
                w(i,j,k)=wstar(i,j,k) - dt/rho*(p(i,j,k)-p(i,j,km))*dxi
            enddo
        enddo
    enddo
   !$acc end kernels 
 
   ! Check divergence (can be skipped)
   !$acc kernels 
   do i=1,nx
        do j=1,nx
            do k=1,nx
                ip=i+1
                jp=j+1
                kp=k+1
                if (ip .gt. nx) ip=1
                if (jp .gt. nx) jp=1
                if (kp .gt. nx) kp=1   
                div(i,j,k) = dxi*(u(ip,j,k)-u(i,j,k) + v(i,jp,k)-v(i,j,k) + w(i,j,kp)-w(i,j,k))
            enddo
        enddo
    enddo
    !$acc end kernels

    !write(*,*) "Maximum value of divergence", maxval(div)

    ! Advance particles
    #if partflag==1
    ! get velocity at particle position, trilinear?
    do i=1,np
      
    enddo
    #endif
    ! avdance velocity and position


    !Check before next time step
    !check courant number
    !$acc kernels
    cou=0.d0
    uc=maxval(u)
    vc=maxval(v)
    wc=maxval(w)
    cou=max(uc*dt*dxi,vc*dt*dxi)
    cou=max(cou,wc*dt*dxi)

    write(*,*) "Courant number             ", cou
    !$acc end kernels

    call cpu_time(tend)
    print '(" Time elapsed = ",f6.1," ms")',1000*(tend-tstart)

enddo


!deallocate
!NS variables
deallocate(u,v,w)
deallocate(ustar,vstar,wstar)
deallocate(rhsu,rhsv,rhsw)
!PFM variables
#if phiflag==1
deallocate(phi,rhsphi,normx,normy,normz)
#endif
!Partciles variables
#if partflag==1
deallocate(xp,vp,ufp,fp)
#endif

end program 
