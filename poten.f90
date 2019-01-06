module poten
  use global
  use IO
  implicit none
  !------------------------------------------
  type t_potential
     double precision,allocatable :: ext(:) ! external potential
     double precision,allocatable :: hartree(:) ! hartreel potential
     double precision,allocatable :: Vx(:) ! exchange potential
     double precision,allocatable :: perturb(:) ! perturbation potential
     double precision,allocatable :: tot(:) ! perturbation potential
     double precision::EX,Ehartree
  end type t_potential
  
contains
    ! --------------------------------------------------------------------------------------
  !
  !             init_pot()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_pot(mesh,param,pot)
    implicit none
    type(t_mesh)::mesh
    type(t_potential)::pot
    type(t_param)::param


    allocate(pot%ext(mesh%N))
    allocate(pot%hartree(mesh%N))
    allocate(pot%Vx(mesh%N))
    pot%hartree=0.0
    pot%Vx=0.0
    allocate(pot%perturb(mesh%N))
    allocate(pot%tot(mesh%N))
    call Vext(mesh,pot%ext)
!    call Vperturb(mesh,pot,param)
    pot%tot=pot%ext+pot%hartree !+pot%perturb
  end subroutine init_pot
  ! --------------------------------------------------------------------------------------
  !
  !              Vext()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Vext(m,pot_ext)
    implicit none
    type(t_mesh) :: m
    double precision :: pot_ext(:)
    double precision :: pts(3),rsqr
    
!    character (len=1024) :: filename
    integer :: i,j,k,nn
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           3D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(m%dim.eq.3) then
       do k=1,m%Nz
          pts(3)=k*m%dz
          do i=1,m%Nx
             pts(1)=i*m%dx
             do j=1,m%Ny
                pts(2)=j*m%dy
                rsqr=(pts(1)-m%center(1))**2+(pts(2)-m%center(2))**2+(pts(3)-m%center(3))**2
                nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
!                pot_ext(nn)=0.5*1.0*rsqr
                pot_ext(nn)=-1.0/sqrt(rsqr)
             end do
          end do
       end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           2D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if(m%dim.eq.2) then
       open(unit=1,file="pot_ext.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          do j=1,m%Ny
             pts(2)=j*m%dy
             rsqr=(pts(1)-m%center(1))**2+(pts(2)-m%center(2))**2
             nn=j+(i-1)*m%Ny
             pot_ext(nn)=0.5*1.0*rsqr
             write(1,*) pts(1),pts(2),pot_ext(nn)
          end do
       end do
       close(1)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           1D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else    if(m%dim.eq.1) then
       open(unit=1,file="pot_ext.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          rsqr=(pts(1)-m%center(1))**2
          pot_ext(i)=.5*1.0*rsqr
          write(1,*) pts(1),pot_ext(i)
       end do
       close(1)
    else
       print *,' STOP in Vext(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
       !stop
  end subroutine Vext

end module poten
