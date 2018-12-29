module mesh_mod
  use global
  implicit none
contains
  ! -----------------------------------------------
  subroutine new_mesh(m,param)
    implicit none
    type(t_mesh)::m
    type(t_param)::param
    double precision:: Lwidth 
    m%dim=param%dim
    Lwidth=param%box_width
    m%Nx=param%Nx
    m%dx=Lwidth/(m%Nx+1)



    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                3D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(m%dim.eq.3) then
       m%Ny=m%Nx
       m%Nz=m%Nx
       m%dy=Lwidth/(m%Ny+1)
       m%dz=Lwidth/(m%Nz+1)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !                2D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if(m%dim.eq.2) then
       m%Ny=m%Nx
       m%Nz=1
       m%dy=Lwidth/(m%Ny+1)
       m%dz=1.0
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !                1D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else   if(m%dim.eq.1) then
       m%Ny=1
       m%Nz=1
       m%dy=1.0
       m%dz=1.0
    else
       print *,' STOP in new_mesh(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
    m%N=m%Nx*m%Ny*m%Nz
    m%dv=m%dx*m%dy*m%dz
    m%center(1)=Lwidth/2
    m%center(2)=Lwidth/2
    m%center(3)=Lwidth/2

    
    ! max number of neighbors. It depends on m%dim:
    allocate(m%n_neighbors(m%N))
    ! 2*m%dim=2  @1D
    ! 2*m%dim=4  @2D
    ! 2*m%dim=6  @3D
    allocate(m%list_neighbors(m%N,2*m%dim)) !
    m%list_neighbors(:,:)=0
    m%n_neighbors(:)=0

    ! n_bound -> number of inactive neighbors of a point
    ! default = 0
    ! max = 3 (corner)
    allocate(m%n_bound(m%N))
    m%n_bound(:)=0
    ! list_bound -> idx of the inactive neighbors. It corresponds to
    !                       bound(:)
    allocate(m%list_bound(m%N,3)) !
    m%list_bound(:,:)=0
    ! number of element in the boundary surface
    m%nbound=8+4*(m%Nx+m%Ny+m%Nz)+&
         2*(m%Nx*m%Ny+m%Nx*m%Nz+&
         m%Ny*m%Nz)
    print *,'new_mesh > nbound=',m%nbound
    m%nbound=         2*(m%Nx*m%Ny+m%Nx*m%Nz+&
         m%Ny*m%Nz)
    print *,'new_mesh > nbound=',m%nbound
    allocate(m%bound(m%nbound))


    call compute_list_neighbors(m)
  end subroutine new_mesh
  ! -----------------------------------------------
  subroutine free_mesh(m)
    implicit none
    type(t_mesh) :: m
    deallocate(m%n_neighbors)
    deallocate(m%list_neighbors) !
    deallocate(m%n_bound)
    deallocate(m%list_bound) !
  end subroutine free_mesh
  ! -----------------------------------------------
  subroutine compute_list_neighbors(m)
    implicit none
    type(t_mesh) :: m
    integer::i,j,k,nn,idx
    !integer,allocatable::n_neighbors(:),list_neighbors(:,:)

    if(m%dim.eq.3) then   ! 3D
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

       idx=1
       call update_bound(idx,1        ,1        ,1,-1,0,0,m)
       call update_bound(idx,1        ,1        ,1,0,-1,0,m)
       call update_bound(idx,1        ,1        ,1,0,0,-1,m)
       call update_bound(idx,m%Nx,1        ,1,1,0,0,m)
       call update_bound(idx,m%Nx,1        ,1,0,-1,0,m)
       call update_bound(idx,m%Nx,1        ,1,0,0,-1,m)
       call update_bound(idx,1        ,m%Ny,1,-1,0,0,m)
       call update_bound(idx,1        ,m%Ny,1,0,1,0,m)
       call update_bound(idx,1        ,m%Ny,1,0,0,-1,m)
       call update_bound(idx,m%Nx,m%Ny,1,1,0,0,m)
       call update_bound(idx,m%Nx,m%Ny,1,0,1,0,m)
       call update_bound(idx,m%Nx,m%Ny,1,0,0,-1,m)
       call update_bound(idx,1        ,1        ,m%Nz,-1,0,0,m)
       call update_bound(idx,1        ,1        ,m%Nz,0,-1,0,m)
       call update_bound(idx,1        ,1        ,m%Nz,0,0,1,m)
       call update_bound(idx,m%Nx,1        ,m%Nz,1,0,0,m)
       call update_bound(idx,m%Nx,1        ,m%Nz,0,-1,0,m)
       call update_bound(idx,m%Nx,1        ,m%Nz,0,0,1,m)
       call update_bound(idx,1        ,m%Ny,m%Nz,-1,0,0,m)
       call update_bound(idx,1        ,m%Ny,m%Nz,0,1,0,m)
       call update_bound(idx,1        ,m%Ny,m%Nz,0,0,1,m)
       call update_bound(idx,m%Nx,m%Ny,m%Nz,1,0,0,m)
       call update_bound(idx,m%Nx,m%Ny,m%Nz,0,1,0,m)
       call update_bound(idx,m%Nx,m%Ny,m%Nz,0,0,1,m)
       do i=2,m%Nx-1
          call update_bound(idx,i,1,1,          0,-1,0,m)
          call update_bound(idx,i,1,1,          0,0,-1,m)
          call update_bound(idx,i,m%Ny,1,        0,1,0,m)
          call update_bound(idx,i,m%Ny,1,         0,0,-1,m)
          call update_bound(idx,i,1,m%Nz,         0,-1,0,m)
          call update_bound(idx,i,1,m%Nz,          0,0,1,m)
          call update_bound(idx,i,m%Ny,m%Nz,  0,1,0,m)
          call update_bound(idx,i,m%Ny,m%Nz,  0,0,1,m)
       end do
       do j=2,m%Ny-1
          call update_bound(idx,1,j,1,-1,0,0,m)
          call update_bound(idx,1,j,1,0,0,-1,m)
          call update_bound(idx,m%Nx,j,1,1,0,0,m)
          call update_bound(idx,m%Nx,j,1,0,0,-1,m)
          call update_bound(idx,1,j,m%Nz,-1,0,0,m)
          call update_bound(idx,1,j,m%Nz,0,0,1,m)
          call update_bound(idx,m%Nx,j,m%Nz,1,0,0,m)
          call update_bound(idx,m%Nx,j,m%Nz,0,0,1,m)
       end do
       do k=2,m%Nz-1
          call update_bound(idx,1,1,k,-1,0,0,m)
          call update_bound(idx,1,1,k,0,-1,0,m)
          call update_bound(idx,m%Nx,1,k,1,0,0,m)
          call update_bound(idx,m%Nx,1,k,0,-1,0,m)
          call update_bound(idx,1,m%Ny,k,-1,0,0,m)
          call update_bound(idx,1,m%Ny,k,0,1,0,m)
          call update_bound(idx,m%Nx,m%Ny,k,1,0,0,m)
          call update_bound(idx,m%Nx,m%Ny,k,0,1,0,m)
       end do
       do i=2,m%Nx-1
          do j=2,m%Ny-1
             call update_bound(idx,i,j,1,0,0,-1,m)
             call update_bound(idx,i,j,m%Nz,0,0,1,m)
          end do
       end do
       do i=2,m%Nx-1
          do k=2,m%Nz-1
             call update_bound(idx,i,1,k,0,-1,0,m)
             call update_bound(idx,i,m%Ny,k,0,1,0,m)
          end do
       end do
       do j=2,m%Ny-1
          do k=2,m%Nz-1
             call update_bound(idx,1,j,k,-1,0,0,m)
             call update_bound(idx,m%Nx,j,k,1,0,0,m)
          end do
       end do



    else    if(m%dim.eq.2) then       ! 2D
       do i=1,m%Nx
          do j=1,m%Ny
             nn=j+(i-1)*m%Ny
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
    else if(m%dim.eq.1) then      ! 1D
       do i=1,m%Nx
          if (i>1) then 
             m%n_neighbors(i)=m%n_neighbors(i)+1
             m%list_neighbors(i,m%n_neighbors(i))=i-1
          end if
          if (i<m%Nx) then 
             m%n_neighbors(i)=m%n_neighbors(i)+1
             m%list_neighbors(i,m%n_neighbors(i))=i+1
          end if
       end do
    else
       print *,' STOP in compute_list_neighbors(): dimension=',m%dim,' not yet implemented!'
       stop
    end if


    

    
  end subroutine compute_list_neighbors
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! update_bound(idx,i,j,k,di,dj,dk,m)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_bound(idx,i,j,k,di,dj,dk,m)
         implicit none
         integer :: idx,i,j,k,nn,di,dj,dk
         type(t_mesh)::m
         nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
         m%n_bound(nn)=m%n_bound(nn)+1
         m%list_bound(nn,m%n_bound(nn))=idx
         m%bound(idx)%q(1)=m%dx*(i+di)
         m%bound(idx)%q(2)=m%dy*(j+dj)
         m%bound(idx)%q(3)=m%dz*(k+dk)
         m%bound(idx)%d=sqrt((m%bound(idx)%q(1)-m%center(1))**2+&
              (m%bound(idx)%q(2)-m%center(2))**2+&
              (m%bound(idx)%q(3)-m%center(3))**2)
         idx=idx+1
       end subroutine update_bound

end module mesh_mod
