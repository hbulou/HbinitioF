module global
  implicit none
  !------------------------------------------
  type t_point
     double precision::q(3) ! coordinate
     double precision::d     ! distance from the center of the cell
     double precision::val   ! value of the Hartree potential
  end type t_point
  !------------------------------------------
  type t_mesh
     integer :: Nx,Ny,Nz,N
     integer,allocatable :: list_neighbors(:,:),n_neighbors(:)
     integer,allocatable :: list_bound(:,:),n_bound(:)  ! list_bound is linked with bound(:) 
     double precision :: dx,dy,dz,dv
     double precision :: center(3)
     integer :: dim
     integer::nbound
     type(t_point),allocatable::bound(:)
  end type t_mesh
     !------------------------------------------
  type t_param
     logical::restart
     character(len=32)::prefix
     logical::init_wf
     logical::extrapol
     integer::extrap_add
     integer::ieof
     integer::loopmax
     integer::nvecmin
     integer::nvecmax
     integer::Nx
     integer::nvec_to_cvg
     integer::noccstate
     double precision :: ETA
     double precision::box_width
     integer:: dim !dimension of the mesh 1(1D), 2(2D) or 3(3D)
     double precision::Iperturb
     double precision::sigma
  end type t_param
  !------------------------------------------
  type t_wavefunction
     double precision,allocatable::S(:)
     double precision,allocatable :: Sprev(:),dS(:) ! eigenvalues
     integer :: nwfc,N
     double precision,allocatable::wfc(:,:)
     double precision,allocatable::rho(:)
  end type t_wavefunction

contains
end module global
