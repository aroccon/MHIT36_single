module param
    integer, parameter :: nx=128
    integer, parameter :: np=0
end module param


module fastp
    use param
    use cufft
    use iso_c_binding
    integer :: cudaplan_fwd,cudaplan_bwd
    double precision, allocatable :: delsq(:,:,:)
    double precision, allocatable :: kk(:)
    double precision, allocatable :: kx(:,:,:), ky(:,:,:), kz(:,:,:)
    real(c_double) :: p(nx,nx,nx), rhsp(nx,nx,nx)
end module fastp


module velocity
   double precision, allocatable :: div(:,:,:)
   double precision, allocatable :: u(:,:,:), v(:,:,:), w(:,:,:)
   double precision, allocatable :: ustar(:,:,:), vstar(:,:,:), wstar(:,:,:)
   double precision, allocatable :: rhsu(:,:,:), rhsv(:,:,:), rhsw(:,:,:)
   double precision, allocatable :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
end module velocity


module phase
   double precision, allocatable :: phi(:,:,:), rhsphi(:,:,:)
   double precision, allocatable :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
end module phase


module particles
   double precision, allocatable :: xp(:,:), vp(:,:), ufp(:,:), fp(:,:)
end module particles
