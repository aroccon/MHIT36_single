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
use velocity
use phase
use particles

#define phiflag 1
#define partflag 0
#define openaccflag 1

implicit none
double precision :: dxi,ddxi
double precision :: rho,mu,rhoi
double precision :: pos,radius,sigma,eps,gamma,umax,normod
double precision :: f1,f2,f3,k0
double precision :: tstart,tend
double precision :: uc,vc,wc
double precision :: h11,h12,h13,h21,h22,h23,h31,h32,h33,cou
integer :: tfin,i,j,k,t,im,jm,km,ip,jp,kp,dump
double precision :: x(nx)


call acc_set_device_num(1,acc_device_nvidia)

! initialize parameters
tfin=2000
dump=100
dt=0.001d0
pi=4.d0*datan(1.d0)
lx=2.d0*pi
dx=lx/(dble(nx)-1)
dxi=1.d0/dx
ddxi=1.d0/dx/dx
rho=1.d0
rhoi=1.d0/rho
mu=0.01d0
radius=1.0d0
sigma=0.1d0
eps=dx
! forcing parameters (ABC)
f1=1.d0
f2=1.d0
f3=1.d0
k0=2.d0


!allocate variables
!NS variables
allocate(u(nx,nx,nx),v(nx,nx,nx),w(nx,nx,nx))
allocate(p(nx,nx,nx),rhsp(nx,nx,nx))
allocate(pc(nx/2+1,nx,nx),rhspc(nx/2+1,nx,nx))
allocate(ustar(nx,nx,nx),vstar(nx,nx,nx),wstar(nx,nx,nx))
allocate(rhsu(nx,nx,nx),rhsv(nx,nx,nx),rhsw(nx,nx,nx))
allocate(fx(nx,nx,nx),fy(nx,nx,nx),fz(nx,nx,nx),div(nx,nx,nx))
!allocate(div(nx,nx,nx))
allocate(delsq(nx,nx,nx))
allocate(kk(nx))!,kx(nx,nx,nx),ky(nx,nx,nx),kz(nx,nx,nx))
!PFM variables
#if phiflag == 1
allocate(phi(nx,nx,nx),rhsphi(nx,nx,nx))
allocate(normx(nx,nx,nx),normy(nx,nx,nx),normz(nx,nx,nx))
allocate(curv(nx,nx,nx),gradphix(nx,nx,nx),gradphiy(nx,nx,nx),gradphiz(nx,nx,nx))
allocate(fxst(nx,nx,nx),fyst(nx,nx,nx),fzst(nx,nx,nx))
#endif
!particles arrays
#if partflag == 1
allocate(xp(np,3),vp(np,3),ufp(np,3),fp(np,3))
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
uc=maxval(u)
vc=maxval(v)
wc=maxval(w)
umax=max(wc,max(uc,vc))
write(*,*) "umax", umax
!$acc end kernels
#if phiflag == 1
write(*,*) 'Initialize phase field'
do k = 1,nx
    do j= 1,nx
        do i = 1,nx
            pos=(x(i)-lx/2)**2d0 +  (x(j)-lx/2)**2d0 + (x(k)-lx/2)**2d0
            phi(i,j,k) = 0.5d0*(1.d0-tanh((sqrt(pos)-radius)/2/eps))
       enddo
    enddo
enddo
#endif

#if partflag == 1
write(*,*) 'Initialize particles'
call random_number(xp)
xp=xp*lx
#endif

!initialize the plan for cuFFT
call init_cufft

!Save initial fields
t=0
call writefield(t,1)
call writefield(t,2)
call writefield(t,3)
call writefield(t,4)
#if phiflag == 1
call writefield(t,5)
#endif
#if partflag == 1 
call writepart(t)
#endif

!use later for FFT (no need to load them afterwards)
!!!$acc enter data create(p,pc)
!!!$acc enter data create(rhsp,rhspc)


! Start temporal loop
do t=1,tfin

    call cpu_time(tstart)
    write(*,*) "Time step",t,"of",tfin


    ! Advance marker function
    ! Compute convective term (A)
    #if phiflag == 1
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
                rhsphi(i,j,k) = - (u(ip,j,k)*0.5d0*(phi(ip,j,k)+phi(i,j,k)) - u(i,j,k)*0.5d0*(phi(i,j,k)+phi(im,j,k)))*dxi  &
                                - (v(i,jp,k)*0.5d0*(phi(i,jp,k)+phi(i,j,k)) - v(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,jm,k)))*dxi  &
                                - (w(i,j,kp)*0.5d0*(phi(i,j,kp)+phi(i,j,k)) - w(i,j,k)*0.5d0*(phi(i,j,k)+phi(i,j,km)))*dxi
            enddo
        enddo
    enddo
    !$acc end kernels

    !Compute diffusive term 
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
                rhsphi(i,j,k)=rhsphi(i,j,k)+gamma*(eps*(phi(ip,j,k)-2.d0*phi(i,j,k)+phi(im,j,k))*ddxi + &
                                                    eps*(phi(i,jp,k)-2.d0*phi(i,j,k)+phi(i,jm,k))*ddxi + &         
                                                    eps*(phi(i,j,kp)-2.d0*phi(i,j,k)+phi(i,j,km))*ddxi)
            enddo
        enddo
    enddo
   !$acc end kernels

    !Compute Sharpening term
    ! Step 1: Compute gradients
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
                normx(i,j,k) = (phi(ip,j,k) - phi(im,j,k))
                normy(i,j,k) = (phi(i,jp,k) - phi(i,jm,k))
                normz(i,j,k) = (phi(i,j,kp) - phi(i,j,km)) 
            enddo
        enddo
    enddo 
    !$acc end kernels

    ! Step 2: Compute normals (1.e-16 is a numerical tolerance)
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                normod = 1.d0/(sqrt(normx(i,j,k)**2d0 + normy(i,j,k)**2d0 + normz(i,j,k)**2d0) + 1.0E-16)
                normx(i,j,k) = normx(i,j,k)*normod
                normy(i,j,k) = normy(i,j,k)*normod
                normz(i,j,k) = normz(i,j,k)*normod
            enddo
        enddo
    enddo
    !$acc end kernels

    !write(*,*) "umax", umax

    ! Compute sharpening term
    !$acc kernels
    gamma=1.d0*umax
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
                rhsphi(i,j,k)=rhsphi(i,j,k)+gamma*(((phi(ip,j,k)**2d0-phi(ip,j,k))*normx(ip,j,k)-(phi(im,j,k)**2d0-phi(im,j,k))*normx(im,j,k))*0.5d0*dxi + &
                                                   ((phi(i,jp,k)**2d0-phi(i,jp,k))*normy(i,jp,k)-(phi(i,jm,k)**2d0-phi(i,jm,k))*normy(i,jm,k))*0.5d0*dxi + &
                                                   ((phi(i,j,kp)**2d0-phi(i,j,kp))*normz(i,j,kp)-(phi(i,j,km)**2d0-phi(i,j,km))*normz(i,j,km))*0.5d0*dxi)
            enddo
        enddo
    enddo
    !$acc end kernels


    ! Compute new phase field n+1
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                phi(i,j,k) = phi(i,j,k) + dt*rhsphi(i,j,k)
            enddo
        enddo
    enddo
    !$acc end kernels
    !write(*,*) "maxvalphi", maxval(phi), "minvalphi", minval(phi)
    if (maxval(phi) .lt. 0.5d0) write(*,*) "Phi is gone"
    #endif


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
                h11=h11*0.25d0*dxi
                h12=h12*0.25d0*dxi
                h13=h13*0.25d0*dxi
                h21=h21*0.25d0*dxi
                h22=h22*0.25d0*dxi
                h23=h23*0.25d0*dxi
                h31=h31*0.25d0*dxi
                h32=h32*0.25d0*dxi
                h33=h33*0.25d0*dxi
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
                h11 = mu*(u(ip,j,k)-2.d0*u(i,j,k)+u(im,j,k))*ddxi
                h12 = mu*(u(i,jp,k)-2.d0*u(i,j,k)+u(i,jm,k))*ddxi
                h13 = mu*(u(i,j,kp)-2.d0*u(i,j,k)+u(i,j,km))*ddxi
                h21 = mu*(v(ip,j,k)-2.d0*v(i,j,k)+v(im,j,k))*ddxi
                h22 = mu*(v(i,jp,k)-2.d0*v(i,j,k)+v(i,jm,k))*ddxi
                h23 = mu*(v(i,j,kp)-2.d0*v(i,j,k)+v(i,j,km))*ddxi
                h31 = mu*(w(ip,j,k)-2.d0*w(i,j,k)+w(im,j,k))*ddxi
                h32 = mu*(w(i,jp,k)-2.d0*w(i,j,k)+w(i,jm,k))*ddxi
                h33 = mu*(w(i,j,kp)-2.d0*w(i,j,k)+w(i,j,km))*ddxi
                rhsu(i,j,k)=rhsu(i,j,k)+(h11+h12+h13)*rhoi
                rhsv(i,j,k)=rhsv(i,j,k)+(h21+h22+h23)*rhoi
                rhsw(i,j,k)=rhsw(i,j,k)+(h31+h32+h33)*rhoi
            enddo
        enddo
    enddo
    !$acc end kernels

    ! forcing term (always x because is the same axis)
    !$acc kernels
    do i=1,nx
        do j=1,nx
            do k=1,nx
                rhsu(i,j,k)= rhsu(i,j,k) + f1*sin(k0*x(k))+f3*sin(k0*x(j))
                rhsv(i,j,k)= rhsv(i,j,k) + f2*sin(k0*x(i))+f1*sin(k0*x(k))
                rhsw(i,j,k)= rhsw(i,j,k) + f3*sin(k0*x(j))+f2*sin(k0*x(i))
            enddo
        enddo
    enddo
    !$acc end kernels

    ! Surface tension forces
    !#if phiflag == 1
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
                curv(i,j,k)=0.5d0*(normx(ip,j,k)-normx(im,j,k))*dxi+0.5d0*(normy(i,jp,k)-normy(i,jm,k))*dxi+0.5d0*(normz(i,j,kp)-normz(i,j,km))*dxi
                gradphix(i,j,k)=0.5d0*(phi(ip,j,k)-phi(im,j,k))*dxi
                gradphiy(i,j,k)=0.5d0*(phi(i,jp,k)-phi(i,jm,k))*dxi
                gradphiz(i,j,k)=0.5d0*(phi(i,j,kp)-phi(i,j,km))*dxi
                fxst(i,j,k)=-6.d0*sigma*curv(i,j,k)*gradphix(i,j,k)*phi(i,j,k)*(1.d0-phi(i,j,k))
                fyst(i,j,k)=-6.d0*sigma*curv(i,j,k)*gradphiy(i,j,k)*phi(i,j,k)*(1.d0-phi(i,j,k))
                fzst(i,j,k)=-6.d0*sigma*curv(i,j,k)*gradphiz(i,j,k)*phi(i,j,k)*(1.d0-phi(i,j,k))
                rhsu(i,j,k)=rhsu(i,j,k) + 0.5d0*(fxst(im,j,k)+fxst(i,j,k))*rhoi
                rhsv(i,j,k)=rhsv(i,j,k) + 0.5d0*(fyst(i,jm,k)+fyst(i,j,k))*rhoi
                rhsw(i,j,k)=rhsw(i,j,k) + 0.5d0*(fzst(i,j,km)+fzst(i,j,k))*rhoi
            enddo
        enddo
    enddo
    !$acc end kernels
    !#endif

    !write(*,*) "rhsu", maxval(rhsu)
    !write(*,*) "rhsv", maxval(rhsv)
    !write(*,*) "rhsw", maxval(rhsw)


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

    ! Compute rhs of Poisson equation div*ustar: divergence at the cell center 
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



    ! call Poisson solver (3DFastPoissons + periodic BCs)
    call poissonfast

    !write(*,*) "maxp", maxval(p)

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
 
   ! Check divergence (can be skipped in production)
    !!$acc kernels 
    !do i=1,nx
    !    do j=1,nx
    !        do k=1,nx
    !            ip=i+1
    !            jp=j+1
    !            kp=k+1
    !            if (ip .gt. nx) ip=1
    !            if (jp .gt. nx) jp=1
    !            if (kp .gt. nx) kp=1   
    !            div(i,j,k) = dxi*(u(ip,j,k)-u(i,j,k) + v(i,jp,k)-v(i,j,k) + w(i,j,kp)-w(i,j,k))
    !        enddo
    !    enddo
    !enddo
    !!$acc end kernels

    !write(*,*) "maxdiv", maxval(div)

    ! Advance particles (get velocity and advance according to particle type)
    #if partflag==1
    call get_velocity
    call move_part
    #endif

    !Check before next time step
    !check courant number
    !$acc kernels
    cou=0.d0
    uc=maxval(u)
    vc=maxval(v)
    wc=maxval(w)
    umax=max(wc,max(uc,vc))
    cou=umax*dt*dxi

    write(*,*) "Courant number             ", cou
    !$acc end kernels

    call cpu_time(tend)
    print '(" Time elapsed = ",f6.1," ms")',1000*(tend-tstart)


     !output fields
     if (mod(t,dump) .eq. 0) then
        write(*,*) "Saving output files"
	call writefield(t,1)
	call writefield(t,2)
	call writefield(t,3)
        call writefield(t,4)
	#if phiflag == 1
	call writefield(t,5)
        #endif
        #if partflag == 1 
        call writepart(t)
        #endif
     endif

enddo


!deallocate
!NS variables
deallocate(u,v,w)
deallocate(p,pc,rhsp,rhspc)
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








