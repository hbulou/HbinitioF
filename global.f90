module global
  implicit none
  !------------------------------------------
  type t_nrj
     double precision::last
     double precision::previous
     double precision::dnrj
  end type t_nrj
  type tt_cvg
     double precision::nrj
     double precision::nrjold
     double precision::dnrj
     double precision:: resi
     logical::cvg
  end type tt_cvg
  type t_cvg
     integer::nwfc
     type(t_nrj)::total_nrj
     type(tt_cvg),allocatable::wfc(:)
     integer,allocatable::list_idx(:)
     double precision :: ETA
     integer :: nvec_to_cvg
     integer :: ncvg
     logical::cvg
  end type t_cvg
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
  type t_wavefunction
     integer :: nwfc
     double precision,allocatable::wfc(:,:)
     double precision,allocatable::eps(:)
     double precision,allocatable ::epsprev(:),deps(:) ! eigenvalues
  end type t_wavefunction
  !------------------------------------------
  type t_molecule
     type(t_wavefunction)::wf       ! wavefunctions of the molecule
!     integer :: N                             ! number of point in the mesh
     type(t_mesh)::mesh               !
     double precision,allocatable::rho(:)
  end type t_molecule

contains
end module global
