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
