program Hbinitio
  implicit none
  type t_GramSchmidt
     integer :: nindep
     integer :: ndep ! number of linear dependencies discovered
  end type t_GramSchmidt
  type(t_GramSchmidt)::GS
  type t_mesh
     integer :: Nx,Ny,Nz,N
     integer,allocatable :: list_neighbors(:,:),n_neighbors(:)
     double precision :: dx,dy,dz,dv
  end type t_mesh
  type(t_mesh) :: mesh
  type t_cvg
     integer,allocatable:: list_cvg(:)
     integer :: ncvg
     double precision :: ETA
     integer :: nvec_to_cvg 
  end type t_cvg
  type(t_cvg) :: cvg
  integer :: nvecini,nvecmax,nvec
  integer,parameter :: seed = 86456
  double precision,allocatable :: V(:,:) ! wavefunctions
  double precision,allocatable :: VRitz(:,:) ! Ritz's vectors
  double precision,allocatable :: T(:,:) ! reduced matrix T
  double precision,allocatable :: S(:), Sprev(:),dS(:) ! eigenvalues
  double precision,allocatable :: residual(:,:) ! residual
  double precision,allocatable :: delta(:,:) ! delta vectors
  double precision,allocatable :: Vnew(:,:) ! Vnew
  integer :: newnvec
  integer :: iloop
  integer :: loopmax
  integer :: ndelta
  real :: start,inter,end,inter2
  integer :: i
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call cpu_time(start)
  call init_mesh(mesh)  

  nvecini=12
  nvecmax=31
  nvec=nvecini
  allocate(V(mesh%N,nvec))
  call init_basis_set(V,nvec,seed,mesh)

  open(unit=1,file="eigenvalues.dat",form='formatted',status='unknown')
  write(1,*)
  close(1)
  open(unit=1,file="dbg.dat",form='formatted',status='unknown')
  write(1,*)
  close(1)
  loopmax=1000
  iloop=1
  cvg%ncvg=0
  cvg%nvec_to_cvg=7
  cvg%ETA=1.0e-4
  allocate(Sprev(nvecini))
  allocate(dS(nvecini))
  Sprev(:)=0.0
  dS(:)=0.0
  do while((iloop.le.loopmax).and.(cvg%ncvg.lt.cvg%nvec_to_cvg))

     write(*,'(A)') 'Main > #######################################'     
     write(*,'(A,I4,A)') 'Main > ############ scf loop=',iloop,' ############'
     write(*,'(A)') 'Main > #######################################'     

     ! T (reduced matrix) computing
     allocate(T(nvec,nvec))
     call cpu_time(inter)
     call compute_T(T,V,nvec,mesh)
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'compute_T')


     ! Diagonatilzation of T
     allocate(S(nvec))
     call cpu_time(inter)
     call diagonalization(S,T,nvec)
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'Diagonalization')

     dS(:)=S(1:nvecini)-Sprev(:)
     Sprev(:)=S(1:nvecini)
     do i=1,nvecini
        write(*,'(A,F12.6,A,E12.2,A)')"Main > Eigenvalues: ",S(i),'(',dS(i),')'
     end do
     call cpu_time(inter)
     open(unit=1,file="eigenvalues.dat",form='formatted',status='unknown',access='append')
     write(1,*) inter,iloop,S(1:nvecini)
     close(1)
     ! computation of the Ritz's vectors
     allocate(VRitz(mesh%N,nvec))
     call cpu_time(inter)
     call Ritz(VRitz,V,T,nvec)
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'Ritz')

     ! computation of residual
     allocate(residual(mesh%N,nvec))
     allocate(cvg%list_cvg(nvec))
     cvg%list_cvg(:)=0
     cvg%ncvg=0

     call cpu_time(inter)
     call compute_residual(residual,VRitz,S,nvec,cvg,mesh)
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'residual')
     ! computation of delta
     allocate(delta(mesh%N,nvec))
     delta(:,:)=0.0
     call cpu_time(inter)

     call compute_delta(delta,residual,S,nvec,cvg,mesh,ndelta)
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'delta')

     deallocate(V)
     allocate(V(mesh%N,nvec+ndelta))
     V(:,1:nvec)=VRitz(:,:)
     print *,'ndelta=',ndelta
     print *,'nvec=',nvec
     V(:,nvec+1:nvec+ndelta)=delta(:,:ndelta)
     allocate(Vnew(mesh%N,nvec+ndelta))
     newnvec=nvec+ndelta
     GS%nindep=nvec
     call cpu_time(inter)
     call GramSchmidt(Vnew,V,newnvec,mesh,GS)    
     call cpu_time(inter2)
     call dbg(iloop,inter,inter2,'GS')

     print *,'Main > ',GS%ndep,newnvec
     
     deallocate(V)
     if(newnvec.le.nvecmax) then
        nvec=newnvec
     else
        print *,'Main > restart from nvecini'
        nvec=nvecini
     end if
     allocate(V(mesh%N,nvec))
     V(:,:)=Vnew(:,1:nvec)

     !call check_ortho(V,nvec,mesh)
     print *,'Main > New size of the basis ',nvec
     iloop=iloop+1
     deallocate(S)
     deallocate(T)
     deallocate(VRitz)
     deallocate(residual)
     deallocate(delta)
     deallocate(Vnew)
     deallocate(cvg%list_cvg)
  end do
  
  call save_cube(V,1,nvecini,mesh)  
  deallocate(V)
  deallocate(Sprev)
  deallocate(dS)
  call free_mesh(mesh)
  call cpu_time(end)
  if (cvg%ncvg.ge.cvg%nvec_to_cvg) print *,'Main > Job DONE !'
  print '("Main > Total Time = ",e16.6," seconds.")',end-start
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  ! -----------------------------------------------
  subroutine dbg(iloop,inter,inter2,text)
    real :: inter,inter2
    integer :: iloop
    character (len=*) :: text
    open(unit=1,file="dbg.dat",form='formatted',status='unknown',access='append')
    write(1,'(A20,I4,F12.6,F12.6,F12.6)') text,iloop,inter,inter2,inter2-inter
    close(1)
  end subroutine dbg
  ! -----------------------------------------------
  subroutine save_cube(evec,idxmin,idxmax,m)
    implicit none
    double precision :: evec(:,:)
    integer :: idxmin,idxmax
    type(t_mesh) :: m
    
    character(len=*),parameter :: FMT1='(I5,3F12.6)'
    integer :: i,j,k,nn,ifield
    integer :: fileid
    character (len=1024) :: filename
    do fileid=idxmin,idxmax
       write(filename,'(a,i0,a)') 'evec',fileid,'.cube'
       open(unit=1,file=filename,form='formatted',status='unknown')
       write(1,*) ' Cubefile created from 3d.f90 calculation'
       write(1,*) ' H. Bulou, October 2018'
       write(1,FMT1) 1,0.0,0.0,0.0
       write(1,FMT1) m%Nx,m%dx,0.0,0.0
       write(1,FMT1) m%Ny,0.0,m%dy,0.0
       write(1,FMT1) m%Nz,0.0,0.0,m%dz
       write(1,'(I5,4F12.6)') 1,1.0,0.0,0.0,0.0
       do k=1,m%Nz
          ifield=0
          do i=1,m%Nx
             do j=1,m%Ny
                nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx               
                write(1,'(E13.5)',advance='no') evec(nn,fileid)
                ifield=ifield+1
                if (mod(ifield,6).eq.0) then
                   ifield=0
                   write(1,*)
                end if
             end do
          end do
          write(1,*)
       end do
       close(1)
    end do
  end subroutine save_cube

  ! -----------------------------------------------
  subroutine compute_delta(delta,r,lambda,nvec,cvg,m,ndelta)
    ! INPUT: the residual |r>, the Ritz's vectors |VRitz>, the eigenvalues lambda
    ! OUTPUT : the correction |delta> to improve the  Ritz's vectors so that to
    !          minimize the residual
    implicit none
    type(t_mesh)::m
    double precision :: lambda(:),r(:,:),delta(:,:)
    integer :: nvec,ndelta
    type(t_cvg)::cvg
    
    double precision, external :: ddot
    double precision, parameter::alpha=1.0,beta=0.0
    double precision, allocatable :: Dinv(:,:),norm
    integer :: i,j,k,kk
    double precision :: deltasqr

    deltasqr=m%dx**2
    print *,'Delta > ---------------------'
    print *,'Delta > --- compute_delta ---'
    print *,'Delta > ---------------------'
    delta(:,:)=0.0
    allocate(Dinv(m%N,m%N))
    Dinv(:,:)=0.0
    ndelta=0
    do i=1,nvec
       
       if(cvg%list_cvg(i).eq.0) then
          ndelta=ndelta+1
          do j=1,m%N
             Dinv(j,j)=1.0/((3.0/deltasqr)-lambda(i))
          end do
          ! see Victor Eijkhout in "Introduction to scientific and technical computing" edited by Willmore et al
          ! Chap 15 Libraries for Linear Algebra
          ! to get a comprehensive way to use dgemv
          call dgemv('N',m%N,m%N,alpha,Dinv,m%N,r(:,ndelta),1,beta,delta(:,ndelta),1)
          !norm=sqrt(ddot(m%N,delta(:,i),1,delta(:,i),1))
          !write(*, '(A10,I4,A2,E12.6)',advance='no') ' Delta > delta(',i,')=',norm
          !delta(:,ndelta)=delta(:,ndelta)+VRitz(:,ndelta)
          norm=1.0/sqrt(ddot(m%N,delta(:,ndelta),1,delta(:,ndelta),1))

          call dscal(m%N,norm,delta(:,ndelta),1)

       end if
       
    end do
    deallocate(Dinv)
    print *,'Delta > ',ndelta,' new vector(s)'
  end subroutine compute_delta
    
  ! -----------------------------------------------
  subroutine compute_residual(r,VRitz,S,nvec,cvg,m)
    implicit none
    type(t_mesh)::m
    integer :: nvec
    double precision :: r(:,:),VRitz(:,:),S(:)
    type(t_cvg) :: cvg

    integer :: i,j,k
    double precision :: norm
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
          r(i,j)=3.0*VRitz(i,j)/deltasqr
          do k=1,m%n_neighbors(i)
             r(i,j)=r(i,j)-0.5*VRitz(m%list_neighbors(i,k),j)/deltasqr
          end do
          r(i,j)=r(i,j)-S(j)*VRitz(i,j)
       end do
       norm=ddot(m%N,r(:,j),1,r(:,j),1)
       write(*,'(A,I4,A,E12.4,A,E12.4)',advance='no') 'Residual > r(',j,')= ',norm,'/',cvg%ETA
       if (norm.lt.cvg%ETA) then
          cvg%ncvg=cvg%ncvg+1
          cvg%list_cvg(j)=1
          write(*,*) '--> converged'
       else
          write(*,*)
       end if
       
    end do
  end subroutine compute_residual
  
  ! -----------------------------------------------
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
  
  ! -----------------------------------------------
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
  
  ! -----------------------------------------------
  subroutine compute_T(T,V,nvec,m)
    implicit none
    double precision,allocatable :: V(:,:),T(:,:)
    integer :: nvec
    type(t_mesh)::m
    
    integer :: i,j,k,l
    double precision :: deltasqr,acc
    double precision, parameter::alpha=0.0
    double precision::beta
    
    deltasqr=m%dx**2
    do i=1,nvec
       do j=1,nvec ! Tij
          T(i,j)=0.0
          do k=1,m%N
             acc=3.0*V(k,j)/deltasqr ! the potential will be added here
             do l=1,m%n_neighbors(k)
                acc=acc-0.5*V(m%list_neighbors(k,l),j)/deltasqr
             end do
             T(i,j)=T(i,j)+V(k,i)*acc
          end do
       end do
    end do
  end subroutine compute_T
  
  ! -----------------------------------------------
  subroutine init_basis_set(V,nvec,seed,m)
    implicit none
    integer :: nvec,seed
    double precision,allocatable :: V(:,:)
    type(t_mesh)::m

    double precision, external :: ddot
    double precision ::norm
    integer :: i,j
    double precision,allocatable :: Vdump(:,:)
    type(t_GramSchmidt) :: GS
    
    allocate(Vdump(m%N,nvec))
    call srand(seed)
    do i=1,nvec
       do j=1,m%N
          Vdump(j,i)=rand()
       end do
    end do
    do i=1,nvec
       norm=ddot(m%N,Vdump(:,i),1,Vdump(:,i),1)
       norm=1.0/sqrt(norm)
       call dscal(m%N,norm,Vdump(:,i),1)
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
    
    integer :: i,k,j,i0
    double precision, parameter :: ETA=1.0e-6
    double precision,allocatable :: a(:)
    double precision :: norm
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
       norm=sqrt(ddot(m%N,Vout(:,GS%nindep+1),1,Vout(:,GS%nindep+1),1))
       !print *,'GS > norm(',i,')=',norm
       if (norm.le.ETA) then
          GS%ndep=GS%ndep+1 ! V(:,nindep+1) is not linearly inependent
       else
          norm=1.0/norm
          call dscal(m%N,norm,Vout(:,GS%nindep+1),1)
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
    integer :: i,j,k,nfail
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
  

  ! -----------------------------------------------
  subroutine init_mesh(m)
    implicit none
    type(t_mesh)::m
    double precision, parameter :: pi=3.1415927
    double precision,parameter :: Lwidth=pi/sqrt(2.0)
 
    m%Nx=15
    m%Ny=m%Nx
    m%Nz=m%Nx
    m%N=m%Nx*m%Ny*m%Nz

    m%dx=Lwidth/(m%Nx+1)
    m%dy=Lwidth/(m%Ny+1)
    m%dz=Lwidth/(m%Nz+1)
    m%dv=m%dx*m%dy*m%dz
    !m%dv=1.0
    
    allocate(m%n_neighbors(m%N))
    allocate(m%list_neighbors(m%N,6)) !
    m%list_neighbors(:,:)=0
    m%n_neighbors(:)=0
    call compute_list_neighbors(m)
  end subroutine init_mesh
  ! -----------------------------------------------
  subroutine free_mesh(m)
    implicit none
    type(t_mesh) :: m
    deallocate(m%n_neighbors)
    deallocate(m%list_neighbors) !
  end subroutine free_mesh
  ! -----------------------------------------------
  subroutine compute_list_neighbors(m)
    implicit none
    type(t_mesh) :: m
    integer::i,j,k,nn
    !integer,allocatable::n_neighbors(:),list_neighbors(:,:)
    
    do k=1,m%Nz
       do i=1,m%Nx
          do j=1,m%Ny
             nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
             if (k>1) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn-m%Nx*m%Ny
             end if
             if (k<m%Nz) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn+m%Nx*m%Ny
             end if
             if (i>1) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn-m%Ny
             end if
             if (i<m%Nx) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn+m%Ny
             end if
             if (j>1) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn-1
             end if
             if (j<m%Ny) then 
                m%n_neighbors(nn)=m%n_neighbors(nn)+1
                m%list_neighbors(nn,m%n_neighbors(nn))=nn+1
             end if
          end do
       end do
    end do

  end subroutine compute_list_neighbors
  
  
end program Hbinitio
