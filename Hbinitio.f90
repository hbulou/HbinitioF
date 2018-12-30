  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Hbinitio
  !$ use OMP_LIB
  use time_tracking
  use global
  use poten
  use IO
  use param_mod
  use mesh_mod
  implicit none
  !  include 'mpif.h'
  !------------------------------------------
  type t_GramSchmidt
     integer :: nindep
     integer :: ndep ! number of linear dependencies discovered
  end type t_GramSchmidt
  type(t_mesh) :: mesh,mesh2
  !------------------------------------------
  type t_cvg
     integer,allocatable:: list_cvg(:)
     integer :: ncvg
     double precision :: ETA
     integer :: nvec_to_cvg
  end type t_cvg
  type(t_cvg) :: cvg
  type(t_time) :: time_spent
  type(t_param)::param,param2
  !------------------------------------------
  type t_perturb
     double precision,allocatable::coeff(:,:)
  end type t_perturb
  type(t_perturb)::perturb
  type (t_potential)::pot,pot2
  type(t_wavefunction):: wf,wf2
  !------------------------------------------
  integer :: i,j,k,n,l
  type t_pseudo
     integer::n
     integer :: npot,npotu
     double precision :: a,b,zval
     double precision,allocatable::r(:)
     double precision,allocatable::pot(:,:)
  end type t_pseudo
  type(t_pseudo) :: pp

  double precision::x,y,z
  !  integer::ierr,my_id,num_procs
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call time_tracking_init(time_spent)
!  call mpi_init(ierr )
!  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
!  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  open(unit=1,file="energy.dat",form='formatted',status='unknown')
  write(1,*) ;  close(1)

  call read_param(param)

!  call read_pp(pp)
  
  call new_mesh(mesh,param)  

  call init_pot(mesh,param,pot)
  


  open(unit=1,file="eigenvalues.dat",form='formatted',status='unknown'); write(1,*);  close(1)
  cvg%nvec_to_cvg=param%nvec_to_cvg
  allocate(perturb%coeff(cvg%nvec_to_cvg,cvg%nvec_to_cvg))
  cvg%ETA=param%ETA

  call new_wf(wf,param,mesh)
!  call numerov(wf,pot,mesh)

  call davidson(param,mesh,cvg,wf,pot,time_spent)


    call read_param(param2)
    param2%dim=param%dim
    param2%box_width=param%box_width
    param2%Nx=param%Nx+10
    param2%init_wf=.FALSE.
    call new_mesh(mesh2,param2)  
    call init_pot(mesh2,param2,pot2)
    call new_wf(wf2,param2,mesh2)

    do k=1,mesh2%Nz
       do i=1,mesh2%Nx
          do j=1,mesh2%Ny
             n=j+(i-1)*mesh2%Ny+(k-1)*mesh2%Ny*mesh2%Nx;             
             x=i*mesh2%dx
             y=j*mesh2%dY
             z=k*mesh2%dZ
             do l=1,param2%nvecmin
                wf2%wfc(n,l)=interpolate(x,y,z,mesh,wf,l)
             end do
          end do
       end do
    end do
!    call save_cube_3D(wf2%wfc(:,1),'essai.cube',mesh2)
 !   stop
    call davidson(param2,mesh2,cvg,wf2,pot2,time_spent)
    
  stop
  !--------------------------------------------------------------------------
  !
  ! perturbation theory
  !
  !--------------------------------------------------------------------------
  call calc_coeff(param,pot,mesh,wf,perturb)
  do i=1,cvg%nvec_to_cvg
     print *,perturb%coeff(i,i),perturb%coeff(i,i)+wf%S(i)
  end do

  do i=1,cvg%nvec_to_cvg
     print *,(perturb%coeff(i,j),j=1,cvg%nvec_to_cvg)
  end do
  !--------------------------------------------------------------------------
  !
  ! end of Hbinitio
  !
  !--------------------------------------------------------------------------
  deallocate(wf%Sprev)
  deallocate(wf%dS)
  deallocate(wf%S)
  deallocate(pot%ext)
  deallocate(pot%perturb)
  deallocate(pot%tot)
  deallocate(perturb%coeff)
  call free_mesh(mesh)
  call cpu_time(time_spent%end)
  if (cvg%ncvg.ge.cvg%nvec_to_cvg) print *,'Main > Job DONE !'
  print '("Main > Total Time = ",e16.6," seconds.")',time_spent%end-time_spent%start
!  call mpi_finalize(ierr)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  ! --------------------------------------------------------------------------------------
  !
  !             new_wf()
  !
  ! --------------------------------------------------------------------------------------
  subroutine new_wf(wf,param,mesh)
    implicit none
    type(t_wavefunction)::wf
    type(t_param)::param
    type(t_mesh)::mesh
    wf%nwfc=param%nvecmin
    wf%N=mesh%N
    allocate(wf%S(wf%nwfc))
    allocate(wf%Sprev(wf%nwfc))
    allocate(wf%dS(wf%nwfc))
    allocate(wf%wfc(wf%N,wf%nwfc))
    allocate(wf%rho(wf%N))
  end subroutine new_wf
  ! --------------------------------------------------------------------------------------
  !
  !             interpolate()
  !
  ! --------------------------------------------------------------------------------------
  function interpolate(x,y,z,mesh,wf,l)
    implicit none
    double precision::x,y,z
    integer :: i,j,k,n,i0,j0,k0,i1,j1,k1
    double precision::x0,y0,z0,x1,y1,z1
    type(t_mesh)::mesh
    type(t_wavefunction)::wf
    integer,parameter::ndim=8
    double precision::A(ndim,ndim),C(ndim)
    integer::info,ipiv(ndim),l
    double precision::interpolate
    
    i0=floor(x/mesh%dx)
    j0=floor(y/mesh%dy)
    k0=floor(z/mesh%dz)

    i1=i0+1
    j1=j0+1
    k1=k0+1
    x0=i0*mesh%dx
    y0=j0*mesh%dy
    z0=k0*mesh%dz
    x1=i1*mesh%dx
    y1=j1*mesh%dy
    z1=k1*mesh%dz

    i=1;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y0 ; A(i,4)=z0 ; A(i,5) = x0*y0 ; A(i,6)=x0*z0 ; A(i,7) = y0*z0 ; A(i,8) = x0*y0*z0
    i=2;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y0 ; A(i,4)=z0 ; A(i,5) = x1*y0 ; A(i,6)=x1*z0 ; A(i,7) = y0*z0 ; A(i,8) = x1*y0*z0
    i=3;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y1 ; A(i,4)=z0 ; A(i,5) = x0*y1 ; A(i,6)=x0*z0 ; A(i,7) = y1*z0 ; A(i,8) = x0*y1*z0
    i=4;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y1 ; A(i,4)=z0 ; A(i,5) = x1*y1 ; A(i,6)=x1*z0 ; A(i,7) = y1*z0 ; A(i,8) = x1*y1*z0
    i=5;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y0 ; A(i,4)=z1 ; A(i,5) = x0*y0 ; A(i,6)=x0*z1 ; A(i,7) = y0*z1 ; A(i,8) = x0*y0*z1
    i=6;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y0 ; A(i,4)=z1 ; A(i,5) = x1*y0 ; A(i,6)=x1*z1 ; A(i,7) = y0*z1 ; A(i,8) = x1*y0*z1
    i=7;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y1 ; A(i,4)=z1 ; A(i,5) = x0*y1 ; A(i,6)=x0*z1 ; A(i,7) = y1*z1 ; A(i,8) = x0*y1*z1
    i=8;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y1 ; A(i,4)=z1 ; A(i,5) = x1*y1 ; A(i,6)=x1*z1 ; A(i,7) = y1*z1 ; A(i,8) = x1*y1*z1

    C=0.0
    i=i0 ; j=j0 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(1)=wf%wfc(n,l)
    i=i1 ; j=j0 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(2)=wf%wfc(n,l)
    i=i0 ; j=j1 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(3)=wf%wfc(n,l)
    i=i1 ; j=j1 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(4)=wf%wfc(n,l)
    i=i0 ; j=j0 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(5)=wf%wfc(n,l)
    i=i1 ; j=j0 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(6)=wf%wfc(n,l)
    i=i0 ; j=j1 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(7)=wf%wfc(n,l)
    i=i1 ; j=j1 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%N) C(8)=wf%wfc(n,l)

    ! do i=1,8
    !    print *,(A(i,j),j=1,8)
    ! end do
    ! do i=1,8
    !    print *,C(i)
    ! end do


!    CALL DGETRS('N' ,  ndim ,  ndim, A , ndim , IPIV, C , ndim , INFO)
    CALL DGESV( ndim , 1, A , ndim , IPIV, C , ndim , INFO)
    ! print *,info
    ! do i=1,8
    !    print *,'C(',i,')=',C(i)
    ! end do
    interpolate=C(1)+C(2)*x+C(3)*y+C(4)*z+C(5)*x*y+C(6)*x*z+C(7)*y*z+C(8)*x*y*z

  end function interpolate
  ! --------------------------------------------------------------------------------------
  !
  !             Numerov()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov(wf,pot,mesh)
    type(t_wavefunction)::wf
    type(t_potential)::pot
    type(t_mesh)::mesh

    integer :: i,j,n_nodes_wanted,idxwfc,nwfc,l
    double precision,allocatable::r(:),rho(:)
    integer::iloop,iloopmax
    character (len=1024) :: filename
    double precision::EHartree
    
    allocate(r(wf%N))
    r(1)=1.0e-12
    do i=2,wf%N
       r(i)=(i-1)*mesh%dx
    end do


    ! V=rR=r^(l+1)*summation
    

    iloopmax=20
    do iloop=1,iloopmax
       l=0
       do i=1,wf%N
          pot%tot(i)=-1.0/r(i)+0.5*l*(l+1)/r(i)**2+pot%hartree(i)
       end do
       print *,minval(pot%tot),maxval(pot%tot)
       
       nwfc=cvg%nvec_to_cvg

       do i=1,nwfc
          n_nodes_wanted=i       
          idxwfc=i
          call numerov_step(n_nodes_wanted,wf,r,cvg,param,idxwfc,pot,l)
       end do
       

       

       open(unit=1,file="energy.dat",form='formatted',status='unknown',access='append')
       write(1,'(I0,F16.8,F16.8,F16.8)') iloop,wf%S(1),EHartree,2*wf%S(1)-EHartree
       close(1)


       
       call serie(wf,mesh,r)
       
       call Hartree_cg(wf,mesh,r,pot)

       ! Hartree energy
       EHartree=0.0
       do i=1,mesh%N-1
          EHartree=EHartree+&
               0.5*mesh%dx*(pot%hartree(i)*wf%wfc(i,1)**2+pot%hartree(i+1)*wf%wfc(i+1,1)**2)
       end do
       
       do i=1,nwfc
          print *,wf%S(i),' Ha=',27.211*wf%S(i),' eV'
       end do


       write(filename,'(a,i0,a)') 'wfc',iloop,'.dat'
       open(unit=1,file=filename,form='formatted',status='unknown')
       do i=1,mesh%N
          write(1,*) r(i),(wf%wfc(i,j),j=1,nwfc)
       end do
       close(1)



       write(filename,'(a,i0,a)') 'potential',iloop,'.dat'
       open(unit=1,file=filename,form='formatted',status='unknown')
       do i=1,mesh%N
          write(1,*) r(i),pot%tot(i),pot%hartree(i),wf%rho(i)
       end do
       close(1)
       
    end do
    
stop
    
    allocate(rho(wf%N))
    rho(1)=0
    do i=2,mesh%N
       rho(i)=-(wf%wfc(i,1)/r(i))**2
    end do

    open(unit=1,file='rho.dat',form='formatted',status='unknown')
    do i=1,mesh%N
       write(1,*) r(i),wf%wfc(i,1),rho(i)
    end do
    close(1)

    
    stop

    

  end subroutine numerov
  ! --------------------------------------------------------------------------------------
  !
  !             Hartree_cg() (conjugate gradient)
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Hartree_cg(wf,mesh,r,pot)
    implicit none
    type(t_potential)::pot
    double precision::r(:)
    type(t_wavefunction)::wf
    type(t_mesh)::mesh

    integer :: i
    double precision :: charge_inf
    double precision,allocatable::source(:)
    double precision,allocatable::U(:)
    double precision,allocatable::b(:)

    
    
    charge_inf=0.0
    do i=1,mesh%N-1
       charge_inf=charge_inf+0.5*mesh%dx*(wf%wfc(i,1)**2+wf%wfc(i+1,1)**2)
    end do
    print *,"charge @ infinite=",charge_inf 

    allocate(source(mesh%N))
    source(1)=0.0
    wf%rho(1)=0.0
    do i=2,mesh%N
       source(i)=-wf%wfc(i+1,1)**2/r(i)
       wf%rho(i)=source(i)/r(i)
    end do

    allocate(U(mesh%N))
    U(1)=0.0
    U(mesh%N)=charge_inf
    do i=2,mesh%N-1
       U(i)=U(1)+(i-1)*(U(mesh%N)-U(1))/(mesh%N-1)
!       U(i)=rand()
    end do

    

    
    allocate(b(mesh%N))
    do i=1,mesh%N
       b(i)=source(i)
    end do
    b(2)=b(2)-U(1)/mesh%dx**2
    b(mesh%N-1)=b(mesh%N-1)-U(mesh%N)/mesh%dx**2


    pot%hartree(1)=0.0
    pot%hartree(mesh%N)=charge_inf
    do i=2,mesh%N
       pot%hartree(i)=pot%hartree(i)*r(i)
    end do

    
    call Conjugate_gradient(-b(2:mesh%N-1),pot%hartree(2:mesh%N-1),mesh%N-2,mesh%dx)    
    
    do i=2,mesh%N
       pot%hartree(i)=pot%hartree(i)/r(i)
    end do
    
    deallocate(b)
    deallocate(source)
    deallocate(U)
  end subroutine Hartree_cg
  ! --------------------------------------------------------------------------------------
  !
  !             Conjugate_gradient()
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Conjugate_gradient(b,x,n,h)
    implicit none
    double precision::b(:),x(:)
    integer::n
    double precision::h

    double precision,allocatable::g(:)
    double precision,allocatable::d(:)
    double precision,allocatable::Ad(:)
    double precision,allocatable::Ax(:)
    double precision::alpha,beta,dAd,f,fold,cvg
    integer::i,iloop,iloopmax
    double precision, external :: ddot
    allocate(g(n))
    allocate(d(n))
    allocate(Ad(n))
    allocate(Ax(n))


    
    
    g(1)=(x(2)-2*x(1))/h**2+b(1)
    do i=2,n-1
       g(i)=(x(i+1)-2*x(i)+x(i-1))/h**2+b(i)
    end do
    g(n)=(x(n-1)-2*x(n))/h**2+b(n)


    call dcopy(n,g,1,d,1) ! g->d

!    alpha=-1.0 ;    call dscal(n,alpha,d,1) ! d-> -d
!    print *,ddot(n,d,1,d,1)
!    stop

    cvg=1.0e-8
    fold=0.0
    f=2*cvg
    iloopmax=10000
    iloop=1
    do while(.not.((iloop.le.iloopmax).and.(abs(f-fold).lt.cvg))) 
!       do iloop=1,iloopmax

       Ad(1)=(d(2)-2*d(1))/h**2
       do i=2,n-1
          Ad(i)=(d(i+1)-2*d(i)+d(i-1))/h**2
       end do
       Ad(n)=(d(n-1)-2*d(n))/h**2

       dAd=ddot(n,d,1,Ad,1)
       alpha=-ddot(n,g,1,d,1)/dAd
       
       call daxpy(n,alpha,d,1,x,1)  ! x+alpha*d -> x



       
       call daxpy(n,alpha,Ad,1,g,1) ! g+alpha*Ad -> g 


       fold=f
       f=ddot(n,b,1,x,1)
       Ax(1)=(x(2)-2*x(1))/h**2
       do i=2,n-1
          Ax(i)=(x(i+1)-2*x(i)+x(i-1))/h**2
       end do
       Ax(n)=(x(n-1)-2*x(n))/h**2
       f=f+0.5*ddot(n,x,1,Ax,1)
!       print *,iloop,f,abs(f-fold),ddot(n,g,1,g,1)


       beta=ddot(n,g,1,Ad,1)/dAd
       call dscal(n,beta,d,1) ! beta*d-> d
       alpha=-1.0
       call daxpy(n,alpha,g,1,d,1) ! -g+d -> d 
       iloop=iloop+1
       if(iloop.eq.iloopmax) then
          print *,"ERROR in COnjugtae_gradient"
          stop
       end if
    end do
    
    deallocate(g)
    deallocate(d)
    deallocate(Ad)
    deallocate(Ax)
  end subroutine Conjugate_gradient
  ! --------------------------------------------------------------------------------------
  !
  !             Conjugate_gradient_3D()
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Conjugate_gradient_3D(b,x,n,h,m)
    implicit none
    double precision::b(:),x(:)
    integer::n
    double precision::h,hsqr
    type(t_mesh)::m

    double precision,allocatable::g(:)
    double precision,allocatable::d(:)
    double precision,allocatable::Ad(:)
    double precision,allocatable::Ax(:)
    double precision::alpha,beta,dAd,f,fold,cvg
    integer::i,iloop,iloopmax,l
    character (len=1024) :: filecube

    double precision, external :: ddot
    allocate(g(n))
    allocate(d(n))
    allocate(Ad(n))
    allocate(Ax(n))
    hsqr=h**2

    
    
    do i=1,n
       g(i)=-2*m%dim*x(i)/hsqr
       do l=1,m%n_neighbors(i)
          g(i)=g(i)+x(m%list_neighbors(i,l))/hsqr
       end do
       g(i)=g(i)+b(i)
    end do



    call dcopy(n,g,1,d,1) ! g->d

!    alpha=-1.0 ;    call dscal(n,alpha,d,1) ! d-> -d
!    print *,ddot(n,d,1,d,1)
!    stop

    cvg=1.0e-8
    fold=0.0
    f=2*cvg
    iloopmax=10000
    iloop=1
    do while(.not.((iloop.le.iloopmax).and.(abs(f-fold).lt.cvg))) 
!       do iloop=1,iloopmax


       do i=1,n
          Ad(i)=-2*m%dim*d(i)/hsqr
          do l=1,m%n_neighbors(i)
             Ad(i)=Ad(i)+d(m%list_neighbors(i,l))/hsqr
          end do
        end do


       dAd=ddot(n,d,1,Ad,1)
       alpha=-ddot(n,g,1,d,1)/dAd
       
       call daxpy(n,alpha,d,1,x,1)  ! x+alpha*d -> x



       
       call daxpy(n,alpha,Ad,1,g,1) ! g+alpha*Ad -> g 


       fold=f
       f=ddot(n,b,1,x,1)
       print *,f-fold
       Ax(1)=(x(2)-2*x(1))/h**2
       do i=2,n-1
          Ax(i)=(x(i+1)-2*x(i)+x(i-1))/h**2
       end do
       Ax(n)=(x(n-1)-2*x(n))/h**2
       f=f+0.5*ddot(n,x,1,Ax,1)
!       print *,iloop,f,abs(f-fold),ddot(n,g,1,g,1)


       beta=ddot(n,g,1,Ad,1)/dAd
       call dscal(n,beta,d,1) ! beta*d-> d
       alpha=-1.0
       call daxpy(n,alpha,g,1,d,1) ! -g+d -> d 
       iloop=iloop+1
       if(iloop.eq.iloopmax) then
          print *,"ERROR in COnjugtae_gradient"
          stop
       end if
    end do

    write(filecube,'(a)') 'hartee.cube'
    call save_cube_3D(x,filecube,m)
    
    deallocate(g)
    deallocate(d)
    deallocate(Ad)
    deallocate(Ax)
  end subroutine Conjugate_gradient_3D
  ! --------------------------------------------------------------------------------------
  !
  !             serie()
  ! It allows us to compute the value of R(@r=0) from u(0)/0
  ! --------------------------------------------------------------------------------------
  subroutine serie(wf,mesh,r)
    implicit none
    type(t_mesh)::mesh
    type(t_wavefunction)::wf
    integer, parameter :: nmax=5
    double precision, dimension(nmax) :: a
    double precision:: r(:)
    double precision, dimension(nmax,nmax) :: b
    double precision :: som
    integer :: i, j,info, lda, ldb, nrhs
    integer, dimension(nmax) :: ipiv
      
      

    do i=1,nmax
       do j=1,nmax
          b(i,j)=r(1+i)**j  
       end do
       a(i)=wf%wfc(1+i,1)
      end do
      nrhs = 1 ! number of right hand sides in b
      lda = nmax  ! leading dimension of a
      ldb = nmax  ! leading dimension of b

      call dgesv(nmax, nrhs, b, lda, ipiv, a, ldb, info)
!      print *,'a=',a

!      print *,'-----------------------------------------'
!      print *,'r ','u ','som'
!      print *,'-----------------------------------------'
      do j=1,nmax
         som=0
         do i=1,nmax
            som=som+a(i)*(r(1+j)**i)
         end do
!         print *,r(1+j),wf%wfc(1+j,1),som
      end do


!    open(unit=1,file='R.dat',form='formatted',status='unknown')
!    write(1,*) 0.0,a(1)
!    do i=2,mesh%N
!       write(1,*) r(i),wf%wfc(i,1)/r(i)
!    end do
!    close(1)



      ! Note: the solution is returned in b
      ! and a has been changed.

    end subroutine serie
  ! --------------------------------------------------------------------------------------
  !
  !             numerov_step()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov_step(n_nodes_wanted,wf,r,cvg,param,idxwfc,pot,l)
    integer :: n_nodes_wanted
    type(t_potential)::pot
    type(t_param)::param
    type(t_cvg)::cvg
    double precision::r(:)
    type(t_wavefunction)::wf
    integer :: iloop
    integer::n_nodes,idxwfc
    double precision::eps
    double precision,allocatable::Q(:),Vin(:),Vout(:),sqrd
    integer :: impt,l
    double precision::emin,emax,dVin,dVout,Iout,Iin,deps,epsold,facsign,eta
    logical,parameter :: outward=.TRUE.,inward=.FALSE.
    ! first we search an eigenenergy close to the soution
    ! by considering the number of nodes of the wavefunction
    allocate(Q(wf%N))
    allocate(Vin(wf%N))
    allocate(Vout(wf%N))
    facsign=(0.5*mod(n_nodes_wanted,2)-1)
    sqrd=mesh%dx**2
    
    emax=min(0.0,maxval(pot%tot))
    emin=1e10
    print *,'------------------------------------------------------------------------------------------'
    print *,"Searching a wfc with ",n_nodes_wanted," nodes" 

    eps=pot%tot(wf%N/2)
    n_nodes=-1
    
    do while(.not.(n_nodes.eq.n_nodes_wanted))
       call compute_Q(Q,wf%N,eps,r,pot)
       Vout(1)=0.0; Vout(2)=0.001
       call numerov_integrate(outward,Q,Vout,wf%N,sqrd)
       n_nodes=count_nodes(Vout,wf%N)
!       write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(0) eps(',iloop,')=',eps,&
!            ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
       if(n_nodes.gt.n_nodes_wanted) then
          emax=eps
          facsign=eps-1.0
       else if(n_nodes.le.n_nodes_wanted) then
          emin=eps
          facsign=eps+1.0
       end if
       if(emin.lt.emax) then
          eps=0.5*(emin+emax)
       else
          eps=facsign
       end if
    end do
!    write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(1) eps(',iloop,')=',eps,&
!         ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
    do while(.not.(n_nodes.eq.(n_nodes_wanted+1)))
       eps=0.5*(emin+emax)
       call compute_Q(Q,wf%N,eps,r,pot)
       Vout(1)=0.0; Vout(2)=0.001
       call numerov_integrate(outward,Q,Vout,wf%N,sqrd)
       n_nodes=count_nodes(Vout,wf%N)
!       write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(2) eps(',iloop,')=',eps,&
!            ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
       if(n_nodes.gt.(n_nodes_wanted+1)) then
          emax=eps
          facsign=eps-1.0
       else if(n_nodes.lt.(n_nodes_wanted+1)) then
          emin=eps
          facsign=eps+1.0
       end if
       if(emin.lt.emax) then
          eps=0.5*(emin+emax)
       else
          eps=facsign
       end if
    end do
!    write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(3) eps(',iloop,')=',eps,&
!         ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
!    stop

    ! Energy levels of Hydrogen atom
    ! n=1 -13.59
    ! n=2   -3.40
    ! n=3  -1.51
    ! n=4  -0.85
    ! n=5  -0.54

    eta=2*cvg%ETA
    iloop=1

    do while((iloop.le.param%loopmax).and.(abs(eta).gt.cvg%ETA))
       eps=0.5*(emin+emax)
       call compute_Q(Q,wf%N,eps,r,pot)
       impt=matching_point(Q,wf%N)
       Vout(1)=0.0; Vout(2)=0.001
       call numerov_integrate(outward,Q,Vout,wf%N,sqrd)
       Vin(wf%N)=0.0; Vin(wf%N-1)=0.001
       call numerov_integrate(inward,Q,Vin,wf%N,sqrd)
       n_nodes=count_nodes(Vout,wf%N)

        Iout=0.0
        do i=1,impt-1
           Iout=Iout+0.5*mesh%dx*(Vout(i)**2+Vout(i+1)**2)
        end do
        Iin=0.0
        do i=impt,wf%N-1
           Iin=Iin+0.5*mesh%dx*(Vin(i)**2+Vin(i+1)**2)
        end do
        dVin=0.5*(Vin(impt+1)-Vin(impt-1))/mesh%dx
        dVout=0.5*(Vout(impt+1)-Vout(impt-1))/mesh%dx
        deps=dVin/Vin(impt)-dVout/Vout(impt)
        deps=deps/(Iout/Vout(impt)**2+Iin/Vin(impt)**2)


       if(iloop.gt.1) eta=eps-epsold
       epsold=eps
       !write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(2) eps(',iloop,')=',eps,&
       !     ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'

 !      write(*,'(A,I4,A,F8.4,A,F8.4,A,F8.4,A,E12.6,A,E12.6,A,I4,A,E12.6,A,E12.6)') 'epsmax(',iloop,')=',eps,&
  !          ' [emin,emax]=[',emin,',',emax,'] eta=',eta,' deps=',deps,' impt=',impt,&
   !         ' dVout=',dVout/Vout(impt),' dVin=',dVin/Vin(impt)


       if(n_nodes.le.(n_nodes_wanted)) then
          emin=eps
       else
          emax=eps
       end if



       iloop=iloop+1
       if(iloop.eq.param%loopmax) then
          print *,"ERROR in numerov_step"
          stop
       end if
    end do
    
!     stop
!     eps=0.5*(emin+emax)
!     call compute_Q(Q,wf%N,eps,r,pot)
!     impt=matching_point(Q,wf%N)
!     iloop=1
!     eta=2*cvg%ETA
!     do while((iloop.le.param%loopmax).and.(abs(eta).gt.cvg%ETA))
!        call compute_Q(Q,wf%N,eps,r,pot)
!        impt=matching_point(Q,wf%N)
!        Vout(1)=0.0; Vout(2)=0.001
!        call numerov_integrate(outward,Q,Vout,wf%N,sqrd)
!        Vin(mesh%N)=0.0; Vin(mesh%N-1)=0.001
!        call numerov_integrate(inward,Q,Vin,wf%N,sqrd)


!            open(unit=1,file="numerov.dat",form='formatted',status='unknown')
!            do i=1,wf%N
!               write(1,*) r(i),Vout(i),Vin(i)
!            end do
!            close(1)

    
!        Iout=0.0
!        do i=1,impt-1
!           Iout=Iout+0.5*mesh%dx*(Vout(i)**2+Vout(i+1)**2)
!        end do
!        Iin=0.0
!        do i=impt,wf%N-1
!           Iin=Iin+0.5*mesh%dx*(Vin(i)**2+Vin(i+1)**2)
!        end do
       
! !       deps=(dVin-dVout)/(Iout/Vout(impt)**2+Iin/Vin(impt)**2)
!        deps=dVin/Vin(impt)-dVout/Vout(impt)
!        deps=deps/(Iout/Vout(impt)**2+Iin/Vin(impt)**2)
!        eta=eps-epsold  
!        write(*,'(A,I4,A,F8.4,A,F8.4,A,F8.4,A,E12.6,A,E12.6,A,I4,A,E12.6,A,E12.6)') 'epsmax(',iloop,')=',eps,&
!             ' [emin,emax]=[',emin,',',emax,'] eta=',eta,' deps=',deps,' impt=',impt,&
!             ' dVout=',dVout/Vout(impt),' dVin=',dVin/Vin(impt)
! !       print *,'eps(',iloop,')=',eps,deps,eta
!        epsold=eps
!        if(deps.lt.0.0) then
!           emax=eps
!        else
!           emin=eps
!        end if
!        eps=0.5*(emin+emax)
!        eps=eps-deps
!        iloop=iloop+1
!     end do
!     print *,'eps(',iloop,')=',eps,deps,eta


    wf%S(idxwfc)=eps
    do i=1,impt
       wf%wfc(i,idxwfc)=Vout(i)/Vout(impt)
    end do
    do i=impt+1,mesh%N
       wf%wfc(i,idxwfc)=Vin(i)/Vin(impt)
    end do
    call norm(mesh,wf%wfc(:,idxwfc))
    
    
    deallocate(Q)
    deallocate(Vin)
    deallocate(Vout)
      

    end subroutine numerov_step

    ! --------------------------------------------------------------------------------------
    !
    !             matching_point()
    !
    ! --------------------------------------------------------------------------------------
    function matching_point(Q,N)
!      function matching_point(Q,N,eps,r,pot)
      implicit none
      integer::N
      double precision::Q(:)!,r(:),eps
      !type(t_potential)::pot
      integer::matching_point,i
      !      call compute_Q(Q,N,eps,r,pot)
      matching_point=-1
      do i=2,N-1
         if((Q(i)*Q(i+1)).le.0) then
            !     print *,i,pot%tot(i),Q(i),Q(i+1)
            if(matching_point.lt.0)          matching_point=i
         end if
      end do
    end function matching_point
    ! --------------------------------------------------------------------------------------
    !
    !             compute_Q()
    !
    ! --------------------------------------------------------------------------------------
    subroutine compute_Q(Q,N,eps,r,pot)
      double precision::Q(:),r(:),eps
      type(t_potential)::pot
    integer::N,i
    Q(1)=10000.0
    do i=2,N
       Q(i)=2*(eps-pot%tot(i))
    end do
    open(unit=1,file="Q.dat",form='formatted',status='unknown')
    do i=1,N
       write(1,*) r(i),Q(i)
    end do
    close(1)

  end subroutine compute_Q
  ! --------------------------------------------------------------------------------------
  !
  !             n_nodes()
  !
  ! --------------------------------------------------------------------------------------
  function count_nodes(V,N)
    double precision::V(:)
    integer::N
    integer::count_nodes
    count_nodes=0
    do i=1,N-1
       if((V(i)*V(i+1)).le.0) count_nodes=count_nodes+1
    end do
  end function count_nodes
  ! --------------------------------------------------------------------------------------
  !
  !             Numerov_integrate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov_integrate(outward,Q,y,N,sqrd)
    logical :: outward
    double precision::Q(:),y(:)
    double precision::t(3),sqrd
    integer :: N,i

    if(outward) then
!       print *,'outward'
       do i=2,N-1
          t(1)=1.0+sqrd*Q(i+1)/12.0
          t(2)=2*(1.0-5.0*sqrd*Q(i)/12.0)
          t(3)=1.0+sqrd*Q(i-1)/12.0
          y(i+1)=(t(2)*y(i)-t(3)*y(i-1))/t(1)
       end do
    else
!       print *,'inward'
       do i=N-1,2,-1
          t(1)=1.0+sqrd*Q(i-1)/12.0
          t(2)=2*(1.0-5.0*sqrd*Q(i)/12.0)
          t(3)=1.0+sqrd*Q(i+1)/12.0
          y(i-1)=(t(2)*y(i)-t(3)*y(i+1))/t(1)
       end do
    end if
  end subroutine numerov_integrate
  ! --------------------------------------------------------------------------------------
  !
  !             V_from_wfc()
  !
  ! --------------------------------------------------------------------------------------
       subroutine V_from_wfc(param,wf,mesh,V)
         double precision, external :: ddot
         integer :: lorb
         type(t_param)::param
         type(t_wavefunction)::wf
         type(t_mesh)::mesh
         double precision::V(:,:),normloc
          do lorb=1,param%nvecmin
             call dcopy(mesh%N,wf%wfc(:,lorb),1,V(:,lorb),1) ! g->d
             normloc=ddot(mesh%N,V(:,lorb),1,V(:,lorb),1)
             normloc=1.0/sqrt(normloc)
             call dscal(mesh%N,normloc,V(:,lorb),1)
          end do
!    GS%nindep=1
!    call GramSchmidt(V,Vdump,nvec,m,GS)

        end subroutine V_from_wfc
  ! --------------------------------------------------------------------------------------
  !
  !             davidson()
  !
  ! --------------------------------------------------------------------------------------
  subroutine davidson(param,mesh,cvg,wf,pot,time_spent)
    type(t_param)::param
    type(t_mesh)::mesh
    type(t_cvg)::cvg
    type(t_potential)::pot
    type(t_wavefunction)::wf
    type(t_time)::time_spent
    
    integer :: nvec
    double precision,allocatable :: V(:,:) ! wavefunctions
!    double precision::normloc
 !   double precision, external :: ddot
    integer :: iloop,i  !,lorb
    !  nvecmin=2
    nvec=param%nvecmin
    allocate(V(mesh%N,nvec))
    if (param%init_wf)   then
       if (.not.(param%restart))   then
          print *,"new calculation"
          call init_basis_set(V,nvec,mesh,wf,'rand')
       else
          print *,'restart an old calculation'
          call read_config(V,mesh,nvec)
       end if
    else
       call init_basis_set(V,nvec,mesh,wf,'wfc')
          !      call check_ortho(V,nvec,mesh)
          !       stop
    end if
    
    iloop=1
    cvg%ncvg=0
    wf%Sprev(:)=0.0
    wf%dS(:)=0.0
    do while((iloop.le.param%loopmax).and.(cvg%ncvg.lt.cvg%nvec_to_cvg))
       write(*,'(A)') 'Main > #######################################'     
       write(*,'(A,I4,A)') 'Main > ############ scf loop=',iloop,' ############'
       write(*,'(A)') 'Main > #######################################'     
       call davidson_step(nvec,V,mesh,param%nvecmin,iloop,cvg,pot,time_spent,wf)

       call compute_density(wf,param,V,mesh)
       call read_param(param)
       cvg%ETA=param%ETA
       
    end do
    call save_config(V,mesh,param%nvecmin)
    
    
    !--------------------------------------------------------------------------
    !
    ! save the wavefunctions
    !
    !--------------------------------------------------------------------------
    call save_wavefunction(param,mesh,V,wf)
    deallocate(V)
  end subroutine davidson
  ! --------------------------------------------------------------------------------------
  !
  !              DAVIDSON_STEP()
  !
  ! --------------------------------------------------------------------------------------
  subroutine davidson_step(nvec,V,m,nvecmin,iloop,cvg,pot,time_spent,wf)
    implicit none
    type(t_wavefunction)::wf
    integer :: nvec,nvecmin,iloop
    double precision,allocatable :: V(:,:)   !,pot_ext(:)
    type(t_mesh) :: m
    type(t_cvg)::cvg
    type(t_time)::time_spent
    type(t_potential)::pot
    
    double precision,allocatable :: S(:)          ! eigenvalues
    double precision,allocatable :: T(:,:)        ! reduced matrix T
    double precision,allocatable :: VRitz(:,:)    ! Ritz's vectors
    double precision,allocatable :: residual(:,:) ! residual
    double precision,allocatable :: delta(:,:)    ! delta vectors
    double precision,allocatable :: Vnew(:,:)     ! Vnew
    integer :: i
    integer :: ndelta
    integer :: newnvec
    type(t_GramSchmidt)::GS
    ! T (reduced matrix) computing
    allocate(T(nvec,nvec))
    call cpu_time(time_spent%start_loc)
    call compute_T(T,V,nvec,m,pot%tot)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> compute_T')
    
    ! Diagonatilzation of T
    allocate(S(nvec))

    call cpu_time(time_spent%start_loc)
    call diagonalization(S,T,nvec)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Diagonalization')

    wf%S(1:nvecmin)=S(1:nvecmin)
    wf%dS(:)=wf%S(1:nvecmin)-wf%Sprev(:)
    wf%Sprev(:)=wf%S(1:nvecmin)
    do i=1,nvecmin
       write(*,'(A,I6,A,F12.6,A,E12.2,A)') 'Main > Eigenvalue(',i,'): ',wf%S(i),'(',wf%dS(i),')'
    end do
    !call cpu_time(inter)
    open(unit=1,file="eigenvalues.dat",form='formatted',status='unknown',access='append')
    write(1,*) iloop,wf%S(1:nvecmin)
    close(1)
    ! computation of the Ritz's vectors
    allocate(VRitz(m%N,nvec))

    call cpu_time(time_spent%start_loc)
    call Ritz(VRitz,V,T,nvec)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Diagonalization')

    
    ! computation of residual
    allocate(residual(m%N,nvec))
    allocate(cvg%list_cvg(nvec))
    cvg%list_cvg(:)=0
    cvg%ncvg=0
    

    call cpu_time(time_spent%start_loc)
    call compute_residual(residual,VRitz,wf%S,nvec,cvg,m,pot%tot)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Residual')

    ! computation of delta
    allocate(delta(m%N,nvec))
    delta(:,:)=0.0
    !call cpu_time(inter)

    call cpu_time(time_spent%start_loc)
    call compute_delta(delta,residual,wf%S,nvec,cvg,m,ndelta,pot%tot)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Delta')

    
    deallocate(V)
    allocate(V(m%N,nvec+ndelta))
    V(:,1:nvec)=VRitz(:,:)
    print *,'ndelta=',ndelta
    print *,'nvec=',nvec
    V(:,nvec+1:nvec+ndelta)=delta(:,:ndelta)
    allocate(Vnew(m%N,nvec+ndelta))
    newnvec=nvec+ndelta
    GS%nindep=nvec
    !call cpu_time(inter)
    call GramSchmidt(Vnew,V,newnvec,m,GS)    
    ! call cpu_time(inter2);     call dbg(iloop,inter,inter2,'GS')
    
    print *,'Main > ',GS%ndep,newnvec
    
    deallocate(V)
    if(newnvec.le.param%nvecmax) then
       nvec=newnvec
    else
       print *,'Main > restart from nvecmin'
       nvec=nvecmin
    end if
    allocate(V(m%N,nvec))
    V(:,:)=Vnew(:,1:nvec)
    
    !call check_ortho(V,nvec,m)
    print *,'Main > New size of the basis ',nvec
    iloop=iloop+1
    deallocate(S)
    deallocate(T)
    deallocate(VRitz)
    deallocate(residual)
    deallocate(delta)
    deallocate(Vnew)
    deallocate(cvg%list_cvg)
  end subroutine davidson_step
  ! --------------------------------------------------------------------------------------
  !
  !              compute_density()
  !
  ! --------------------------------------------------------------------------------------
       subroutine compute_density(wf,param,V,mesh)
         implicit none
         type(t_wavefunction)::wf
         type(t_param)::param
         type(t_mesh)::mesh
         double precision::V(:,:)
         integer::i,j
         double precision::charge
         double precision,allocatable::b(:)

         do j=1,param%noccstate
            call dcopy(wf%N,V(:,j),1,wf%wfc(:,j),1)
            call norm(mesh,wf%wfc(:,j))
         end do



         wf%rho=0.0
         do j=1,param%noccstate
            do i=1,mesh%N
               wf%rho(i)=wf%rho(i)-wf%wfc(i,j)**2
            end do
         end do
         charge=mesh%dv*sum(wf%rho)
         print *,"Compute Density > Charge ",charge

         do i=1,mesh%nbound
            mesh%bound(i)%val=charge/mesh%bound(i)%d
            !print *,mesh%bound(i)%q,mesh%bound(i)%val
         end do
         

    allocate(b(mesh%N))
    do i=1,mesh%N
       b(i)=wf%rho(i)
       do j=1,mesh%n_bound(i)
          b(i)=b(i)+mesh%bound(mesh%list_bound(i,j))%val/mesh%dx**2
       end do
!       print *,b(i)
    end do
!    b(mesh%N-1)=b(mesh%N-1)-U(mesh%N)/mesh%dx**2
!    call Conjugate_gradient_3D(-b,pot%hartree,mesh%N,mesh%dx,mesh)    
  end subroutine compute_density
  ! --------------------------------------------------------------------------------------
  !
  !              calc_coeff()
  !
  ! --------------------------------------------------------------------------------------
  subroutine calc_coeff(param,pot,mesh,wf,perturb)
    implicit none
    type(t_perturb)::perturb
    type(t_wavefunction)::wf
    type(t_mesh)::mesh
    type(t_potential)::pot
    type(t_param)::param
    integer::i,j
    do i=1,param%nvec_to_cvg
       do j=1,param%nvec_to_cvg
          perturb%coeff(i,j)=mesh%dv*sum(wf%wfc(:,i)*pot%perturb*wf%wfc(:,j))
       end do
    end do
  end subroutine calc_coeff
  ! --------------------------------------------------------------------------------------
  !
  !                             save_config()
  !
  ! subroutine to save the configuration of the calculation in order to restart it
  ! later if necessary
  ! --------------------------------------------------------------------------------------
  subroutine save_config(V,m,nvecmin)
    implicit none
    type(t_mesh)::m
    double precision :: V(:,:)
    integer::nvecmin,i,j
    open(unit=1,file="evectors.dat",form='formatted',status='unknown')
    do i=1,m%N
       write(1,*) (V(i,j),j=1,nvecmin)
    end do
    close(1)
  end subroutine save_config
  ! --------------------------------------------------------------------------------------
  !
  !              read_config()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_config(V,m,nvecmin)
    implicit none
    type(t_mesh)::m
    double precision :: V(:,:)
    integer::nvecmin,i,j
    logical :: file_exists
    INQUIRE(FILE="evectors.dat", EXIST=file_exists)
    if(file_exists) then
       open(unit=1,file="evectors.dat",form='formatted',status='unknown')
       do i=1,m%N
          read(1,*) (V(i,j),j=1,nvecmin)
       end do
       close(1)
    else
       print *,"### ERROR ### evectors.dat doesn't exist"
       stop
    end if
  end subroutine read_config
  ! --------------------------------------------------------------------------------------
  !
  !              Vperturb()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Vperturb(m,pot,param)
    implicit none
    type(t_mesh) :: m
    type(t_param)::param
    type(t_potential)::pot
    double precision :: pts(3),rsqr

    double precision, parameter :: pi=3.1415927
    double precision :: invsig
    double precision :: facperturb
    integer :: i,j,nn

    facperturb=param%Iperturb/sqrt(2*pi*param%sigma**2)
    invsig=0.5/param%sigma**2
    if(m%dim.eq.2) then
       open(unit=1,file="pot_perturb.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          do j=1,m%Ny
             pts(2)=j*m%dy
             rsqr=(pts(1)-m%center(1))**2+(pts(2)-m%center(2))**2
             nn=j+(i-1)*m%Ny
             pot%perturb(nn)=facperturb*exp(-invsig*rsqr)
             write(1,*) pts(1),pts(2),pot%perturb(nn)
          end do
       end do
       close(1)
    else        if(m%dim.eq.1) then
       open(unit=1,file="pot_perturb.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          rsqr=(pts(1)-m%center(1))**2
          pot%perturb(i)=facperturb*exp(-invsig*rsqr)
          write(1,*) pts(1),pot%perturb(i)
       end do
       close(1)
    else
       print *,' STOP in Vperturb(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
  end subroutine Vperturb
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_DELTA()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_delta(delta,r,lambda,nvec,cvg,m,ndelta,pot_tot)
    ! INPUT: the residual |r>, the Ritz's vectors |VRitz>, the eigenvalues lambda
    ! OUTPUT : the correction |delta> to improve the  Ritz's vectors so that to
    !          minimize the residual
    implicit none
    type(t_mesh)::m
    double precision :: lambda(:),r(:,:),delta(:,:),pot_tot(:)
    integer :: nvec,ndelta
    type(t_cvg)::cvg
    
    double precision, external :: ddot
    double precision, parameter::alpha=1.0,beta=0.0
    double precision, allocatable :: normloc
!    double precision, allocatable :: Dinv(:,:)
    integer :: i,j
    double precision :: deltasqr

    deltasqr=m%dx**2
    print *,'Delta > ---------------------'
    print *,'Delta > --- compute_delta ---'
    print *,'Delta > ---------------------'
    !delta(:,:)=0.0
   ! allocate(Dinv(m%N,m%N))
  !  Dinv(:,:)=0.0
    ndelta=0
    do i=1,nvec
       if(cvg%list_cvg(i).eq.0) then
          ndelta=ndelta+1
          do j=1,m%N
             delta(j,ndelta)=r(j,ndelta)/((m%dim/deltasqr+pot_tot(j))-lambda(i))
!             Dinv(j,j)=1.0/((3.0/deltasqr+pot_tot(j))-lambda(i))
          end do
          ! see Victor Eijkhout in "Introduction to scientific and technical computing" edited by Willmore et al
          ! Chap 15 Libraries for Linear Algebra
          ! to get a comprehensive way to use dgemv
 !         call dgemv('N',m%N,m%N,alpha,Dinv,m%N,r(:,ndelta),1,beta,delta(:,ndelta),1)
          !norm=sqrt(ddot(m%N,delta(:,i),1,delta(:,i),1))
          !write(*, '(A10,I4,A2,E12.6)',advance='no') ' Delta > delta(',i,')=',norm
          !delta(:,ndelta)=delta(:,ndelta)+VRitz(:,ndelta)
          normloc=1.0/sqrt(ddot(m%N,delta(:,ndelta),1,delta(:,ndelta),1))
          call dscal(m%N,normloc,delta(:,ndelta),1)
       end if
    end do
    !deallocate(Dinv)
    print *,'Delta > ',ndelta,' new vector(s)'
  end subroutine compute_delta
    
  
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_RESIDUAL()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_residual(r,VRitz,S,nvec,cvg,m,pot_tot)
    implicit none
    type(t_mesh)::m
    integer :: nvec
    double precision :: r(:,:),VRitz(:,:),S(:),pot_tot(:)
    type(t_cvg) :: cvg

    integer :: i,j,k
    double precision :: normloc
    double precision, external :: ddot        
    double precision :: deltasqr

    print *,'Residual > ------------------------'
    print *,'Residual > --- compute residual ---'
    print *,'Residual > ------------------------'
    deltasqr=m%dx**2
    r(:,:)=0.0
    cvg%ncvg=0
    do j=1,nvec
       do i=1,m%N
          r(i,j)=(m%dim/deltasqr+pot_tot(i))*VRitz(i,j)
          do k=1,m%n_neighbors(i)
             r(i,j)=r(i,j)-0.5*VRitz(m%list_neighbors(i,k),j)/deltasqr
          end do
          r(i,j)=r(i,j)-S(j)*VRitz(i,j)
       end do
       normloc=ddot(m%N,r(:,j),1,r(:,j),1)
       write(*,'(A,I4,A,E12.4,A,E12.4)',advance='no') 'Residual > r(',j,')= ',normloc,'/',cvg%ETA
       if (normloc.lt.cvg%ETA) then
          cvg%ncvg=cvg%ncvg+1
          cvg%list_cvg(j)=1
          write(*,*) '--> converged'
       else
          write(*,*)
       end if
       
    end do
  end subroutine compute_residual
  
  ! --------------------------------------------------------------------------------------
  !
  !              RITZ()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Ritz(Vout,Vin,y,nvec)
    implicit none
    double precision :: Vin(:,:),Vout(:,:),y(:,:)
    integer :: nvec

    integer :: i,j
    print *,'Ritz > --------------'
    print *,'Ritz > --- Ritz() ---'
    print *,'Ritz > --------------'
    do i=1,nvec
       Vout(:,i)=y(1,i)*Vin(:,1)
       do j=2,nvec
          Vout(:,i)=Vout(:,i)+y(j,i)*Vin(:,j)
       end do
    end do
  end subroutine Ritz
  
  
  ! --------------------------------------------------------------------------------------
  !
  !              DIAGONALIZATION()
  !
  ! --------------------------------------------------------------------------------------
  subroutine diagonalization(S,H,N)
    implicit none
    integer :: N
    double precision :: H(:,:),S(:)

    integer :: lwork,info
    integer :: lwmax
    double precision,allocatable::work(:)
    parameter(lwmax=100000)
    allocate(work(lwmax))
    lwork=-1
    call  dsyev('vectors','Upper',N,H,N,S,work,lwork,info)
    lwork=min(lwmax,int(work(1)))
    if (lwork.ge.lwmax) then
       write(*,*) 'Diagonalization > WARNING info = ',info
       write(*,*) 'Diagonalization > WARNING lwork=',lwork
       write(*,*) 'Diagonalization > WARNING size of work(1)',int(work(1))
       stop
    end if
    call  dsyev('vectors','Upper',N,H,N,S,work,lwork,info)
    if(info.gt.0) then
       write(*,*) "Diagonalization > WARNING The algorithm computing failed"
       stop
    end if
    deallocate(work)
  end subroutine diagonalization
  
  
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_T()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_T(T,V,nvec,m,pot_tot)
    implicit none
    double precision,allocatable :: V(:,:),T(:,:),pot_tot(:)
    integer :: nvec
    type(t_mesh)::m
    
    integer :: i,j,k,l
    double precision :: deltasqr,acc
    double precision, parameter::alpha=0.0
!    double precision::beta
    deltasqr=m%dx**2
    !$OMP PARALLEL private(acc) 
    !$OMP DO 
    do j=1,nvec
       do i=1,nvec ! Tij
          T(i,j)=0.0
          do k=1,m%N
             acc=(m%dim/deltasqr+pot_tot(k))*V(k,j) ! the potential will be added here
             do l=1,m%n_neighbors(k)
                acc=acc-0.5*V(m%list_neighbors(k,l),j)/deltasqr
!                print *, omp_get_thread_num(),i,j,k,l
             end do
             T(i,j)=T(i,j)+V(k,i)*acc
!             print *,' -> ',omp_get_thread_num(),i,j,T(i,j)
          end do
!          print *,omp_get_thread_num(),i,j,T(i,j)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
!    stop

  end subroutine compute_T
  
  
  ! --------------------------------------------------------------------------------------
  !
  !              INIT_BASIS_SET()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_basis_set(V,nvec,m,wf,tag)
    implicit none
    integer :: nvec
    double precision,allocatable :: V(:,:)
    type(t_mesh)::m
    type(t_wavefunction)::wf
    character(len=*)::tag

    double precision, external :: ddot
    double precision ::normloc
    integer :: i,j
    double precision,allocatable :: Vdump(:,:)
    type(t_GramSchmidt) :: GS
    integer,parameter :: seed = 86456
    
    allocate(Vdump(m%N,nvec))

    if(tag.eq.'rand') then
       print *,'Init_basis_set > Random initialization'
       call srand(seed)
       do i=1,nvec
          do j=1,m%N
             Vdump(j,i)=rand()
          end do
       end do
    else
       print *,'Init_basis_set >  wfc initialization'
       do i=1,nvec
          call dcopy(mesh%N,wf%wfc(:,i),1,Vdump(:,i),1) ! g->d
       end do
    end if

    do i=1,nvec
       normloc=ddot(m%N,Vdump(:,i),1,Vdump(:,i),1)
       normloc=1.0/sqrt(normloc)
       call dscal(m%N,normloc,Vdump(:,i),1)
    end do
    GS%nindep=1
    call GramSchmidt(V,Vdump,nvec,m,GS)

    deallocate(Vdump)
  end subroutine init_basis_set
  ! -----------------------------------------------
  ! Ref.: D. G. Clayton "Gram-Schmidt Orthogonalization", J. Roy. Stat. Soc. C 20, 335 (1971)
  subroutine GramSchmidt(Vout,Vin,nvec,m,GS)
    implicit none
    integer ::  nvec
    double precision,allocatable :: Vin(:,:),Vout(:,:)
    type(t_mesh)::m
    type(t_GramSchmidt) :: GS
    
    integer :: i,k,i0
    double precision, parameter :: ETA=1.0e-6
    double precision,allocatable :: a(:)
    double precision :: normloc
    double precision, external :: ddot

    allocate(a(nvec))
    print *,'GS> ----------------------'
    print *,"GS> Gram-Schmidt algorithm"
    print *,'GS> ----------------------'
    print *,'GS> ',nvec,' vectors to orthogonalize'
    print *,'GS> ',GS%nindep,' vectors are already orthogonalized'
    Vout(:,1:GS%nindep)=Vin(:,1:GS%nindep)
    GS%ndep=0 
    i0=GS%nindep
    do i=i0+1,nvec
       Vout(:,GS%nindep+1)=Vin(:,i)
       do k=1,GS%nindep
          ! We compute the projection of Vini(:,i) on V(:,1-nindep)
          a(k)=ddot(m%N,Vout(:,k),1,Vin(:,i),1)
          !print *,'GS > ',i,k,a(k)
          ! then we remove V(:,k) from V(:,nindep+1)
          call daxpy(m%N,-a(k),Vout(:,k),1,Vout(:,GS%nindep+1),1)
       end do
       ! now wre compute the norm of V(:,nindep+1)
       normloc=sqrt(ddot(m%N,Vout(:,GS%nindep+1),1,Vout(:,GS%nindep+1),1))
       !print *,'GS > norm(',i,')=',norm
       if (normloc.le.ETA) then
          GS%ndep=GS%ndep+1 ! V(:,nindep+1) is not linearly inependent
       else
          normloc=1.0/normloc
          call dscal(m%N,normloc,Vout(:,GS%nindep+1),1)
          !              do k=2,icur
          !                 print *,'<U',k-1,'|U',i,'>=',ddot(N,V(:,k-1),1,V(:,icur),1)
          !              end do
          GS%nindep=GS%nindep+1
       end if
    end do
    print *,'GS> ',GS%ndep,' vectors linearly dependant'
    print *,'GS> ',GS%nindep,' vectors linearly independant'
    print *,'GS> Size of the basis from ',nvec,' to ',GS%nindep
    !call check_ortho(Vout,nvec,m)
    !stop
    nvec=GS%nindep
    deallocate(a)
  end subroutine GramSchmidt
  ! -----------------------------------------------
  subroutine check_ortho(P,nvec,m)
    implicit none
    integer :: nvec
    double precision :: P(:,:)
    type(t_mesh)::m

    double precision, parameter :: ETA=1.0e-6
    integer :: i,j,nfail
    double precision :: pscal
    double precision, external :: ddot
    nfail=-nvec
    print *,"--- check_ortho() ---"      
    do i=1,nvec
       do j=1,nvec
          pscal=ddot(m%N,P(:,i),1,P(:,j),1)
          write(*,'(F10.2)',advance='no') pscal
          if (pscal.gt.ETA) nfail=nfail+1
       end do
       write(*,*)
    end do
    if (nfail.gt.0) then
       print *,nfail,' fail(s)'
       stop
    end if
  end subroutine check_ortho
  

  ! --------------------------------------------------------------------------------------
  !
  !             read_pp()
  !
  ! read pseudpopotential file
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_pp(pp)
    type(t_pseudo)::pp
    integer ::i,j
    open(unit=1,file='/home/bulou/ownCloud/src/Octopus/octopus-7.1/share/pseudopotentials/PSF/H.psf',&
         action="read",form='formatted',status='old')
    read(1, *) 
    read(1, *)
    read(1, *) 
    read(1,*) pp%npot, pp%npotu, pp%n, pp%b, pp%a, pp%zval
    print *, pp%npot, pp%npotu, pp%n, pp%b, pp%a, pp%zval
    ! add extra point for the zero
    pp%n=pp%n+1
    print *,pp%n
    allocate(pp%r(pp%n))
    allocate(pp%pot(pp%npot,pp%n))
    pp%r(1)=0.0
    read(1, *)
    read(1,*) (pp%r(i),i=2,pp%n)
    do i=1,pp%npot
       read(1,*)
       pp%pot(i,1)=0.0
       read(1,*) (pp%pot(i,j),j=2,pp%n)
    end do
    ! do i=2,pp%n
    !    read(1,*) pp%r(i)
    !    print *,i,pp%r(i)
    !    j=j+1
    !    if(j.eq.4) then
    !       read(1,*)
    !       j=0
    !    end if
    ! end do

    close(1)

    do i=1,pp%npot
       pp%pot(i,2:pp%n)=pp%pot(i,2:pp%n)/pp%r(2:pp%n)
    end do
    open(unit=1,file='pot.dat',form='formatted',status='unknown')
    do i=1,pp%n
       write(1,*) pp%r(i),(pp%pot(j,i),j=1,pp%npot)
    end do
    close(1)



  end subroutine read_pp
  ! --------------------------------------------------------------------------------------
  !
  !             Hartree_sd() (steepest descent)
  ! 
  ! --------------------------------------------------------------------------------------
  ! subroutine Hartree2()
  !   implicit none
  !   integer ::n,i,j
  !   double precision::alpha,dr,rmax
  !   double precision,allocatable::source(:),b(:)
  !   double precision,allocatable::U(:),grad(:),y(:)
  !   double precision,allocatable::r(:)
  !   integer,parameter :: seed = 86456
  !   double precision, external :: ddot
  !   call srand(seed)
  !   n=1000
  !   rmax=50.0
  !   dr=rmax/(n-1)
  !   allocate(source(n))
  !   allocate(r(n))
  !   allocate(b(n))
  !   allocate(U(n))
  !   allocate(grad(n))
  !   allocate(y(n))
  !   do i=1,n
  !      r(i)=(i-1)*dr
  !   end do
  !   do i=1,n
  !      source(i)=-4.0*r(i)*exp(-2.0*r(i))
  !   end do

    
  !   U(1)=0
  !   do i=2,n-1
  !      U(i)=rand()
  !   end do
  !   U(n)=1.0

  !   do i=2,n-1
  !      b(i)=source(i)
  !   end do
  !   b(2)=source(2)-U(1)/dr**2
  !   b(n-1)=source(n-1)-U(n)/dr**2

  !   grad(1)=0.0
  !   grad(n)=0.0
  !   y(1)=0.0
  !   y(n)=0.0
  !   do j=1,1000000
  !      grad(2)=(U(3)-2*U(2))/dr**2-b(2)
  !      do i=3,n-2
  !         grad(i)=(U(i+1)-2*U(i)+U(i-1))/dr**2-b(i)
  !      end do
  !      grad(n-1)=(-2*U(n-1)+U(n-2))/dr**2-b(n-1)


  !      y(2)=(grad(3)-2*grad(2))/dr**2
  !      do i=3,n-2
  !         y(i)=(grad(i+1)-2*grad(i)+grad(i-1))/dr**2
  !      end do
  !      y(n-1)=(-2*grad(n-1)+grad(n-2))/dr**2
       
  !      alpha=ddot(n,grad,1,grad,1)/ddot(n,grad,1,y,1)
  !      print *,'alpha=',alpha ,ddot(n,grad,1,y,1),ddot(n,grad,1,grad,1)
       
  !      !    call dscal(mesh%N-2,-alpha,grad,1)
  !      call daxpy(n,-alpha,grad,1,U,1)

       
  !   end do
  !   U(1)=0.0
  !   U(n)=1.0


  !   open(unit=1,file='b.dat',form='formatted',status='unknown')
  !   do i=1,n
  !      write(1,*) r(i),source(i),-(r(i)+1)*exp(-2.0*r(i))+1,U(i)
  !   end do
    
  ! end subroutine Hartree2
  subroutine Hartree_sd(wf,mesh,r)
    implicit none
    double precision::r(:)
    type(t_wavefunction)::wf
    type(t_mesh)::mesh
    double precision,allocatable::source(:),U(:),grad(:),y(:),b(:)
    double precision::alpha,q
    integer::i,j
    double precision, external :: ddot
    integer,parameter :: seed = 86456
    character (len=1024) :: filesave

    allocate(source(wf%N))
    allocate(U(wf%N))
    allocate(b(wf%N))
    allocate(grad(wf%N))
    allocate(y(wf%N))
    call srand(seed)

    q=0.0
    do i=1,mesh%N-1
       q=q+0.5*mesh%dx*(wf%wfc(i,1)**2+wf%wfc(i+1,1)**2)
    end do
    print *,"q=",q 

    do i=2,mesh%N-1
       U(i)=rand()
    end do
    U(1)=0.0
    U(mesh%N)=q


    do i=1,mesh%N
       source(i)=-wf%wfc(i+1,1)**2/r(i+1)
    end do


    do i=1,mesh%N
       b(i)=source(i)
    end do
    b(2)=b(2)-U(1)/mesh%dx**2
    b(mesh%N-1)=b(mesh%N-1)-U(mesh%N)/mesh%dx**2



    
    grad(1)=0.0
    grad(mesh%N)=0.0
    y(1)=0.0
    y(mesh%N)=0.0

    
    do j=1,100000000
       grad(2)=(U(3)-2*U(2))/mesh%dx**2-b(2)
       do i=3,mesh%N-2
          grad(i)=(U(i+1)-2*U(i)+U(i-1))/mesh%dx**2-b(i)
       end do
       grad(mesh%N-1)=(-2*U(mesh%N-1)+U(mesh%N-2))/mesh%dx**2-b(mesh%N-1)
       
       y(2)=(grad(3)-2*grad(2))/mesh%dx**2
       do i=3,mesh%N-2
          y(i)=(grad(i+1)-2*grad(i)+grad(i-1))/mesh%dx**2
       end do
       y(mesh%N-1)=(-2*grad(mesh%N-1)+grad(mesh%N-2))/mesh%dx**2
       
       alpha=ddot(mesh%N,grad,1,grad,1)/ddot(mesh%N,grad,1,y,1)
       print *,'alpha=',alpha,ddot(mesh%N,grad,1,y,1),ddot(mesh%N,grad,1,grad,1)
       
       !    call dscal(mesh%N-2,-alpha,grad,1)
       call daxpy(mesh%N,-alpha,grad,1,U,1)

       if(mod(j,10000).eq.0) then
          write(filesave,'(a,i0,a)') 'b',j,'.dat'
          open(unit=1,file=filesave,form='formatted',status='unknown')
          do i=1,mesh%N
             write(1,*) r(i),source(i-1),U(i)
          end do
       end if
       
    end do


!    open(unit=1,file='b.dat',form='formatted',status='unknown')
!    do i=1,mesh%N
!       write(1,*) r(i),source(i-1),U(i)
!    end do

!    write(1,*) r(mesh%N),source(i-1),x(i)
!    close(1)

    
    deallocate(y)
    deallocate(grad)
    deallocate(source)
    deallocate(U)
  end subroutine Hartree_sd
  
  
end program Hbinitio
