module IO
  use global
  use tools
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              SAVE_CUBE()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_cube_3D(data,filename,m)
    implicit none
    double precision :: data(:)
!    integer :: idxmin,idxmax
    type(t_mesh) :: m
    character (len=1024) :: filename
    
    character(len=*),parameter :: FMT1='(I5,3F12.6)'
    integer :: i,j,k,nn,ifield
    
    open(unit=1,file=filename,form='formatted',status='unknown')
    write(1,*) ' Cubefile created from Hbinitio.f90 calculation'
    write(1,*) ' H. Bulou, November 2018'
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
             write(1,'(E13.5)',advance='no') data(nn)
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
  end subroutine save_cube_3D
  ! --------------------------------------------------------------------------------------
  !
  !              save_wavefunction(param,mesh,V,wf)
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_wavefunction(param,mesh,V,wf)
    implicit none
    type(t_param)::param
    type(t_mesh)::mesh
    double precision::V(:,:)
    type(t_wavefunction)::wf
    
    integer :: i,j,k
    character (len=1024) :: filecube
    
    do i=1,param%nvecmin
       call norm(mesh,V(:,i))
       call dcopy(wf%N,V(:,i),1,wf%wfc(:,i),1)
       if(mesh%dim.eq.3) then    ! 3D
          write(filecube,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/evec',i,'.cube'
          call save_cube_3D(V(:,i),filecube,mesh)
       else if(mesh%dim.eq.2) then   ! 2D
          write(filecube,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/evec',i,'.dat'
          open(unit=1,file=filecube,form='formatted',status='unknown')
          do j=1,mesh%Nx
             do k=1,mesh%Ny
                write(1,*) j*mesh%dx,k*mesh%dy,V(j+(k-1)*mesh%Nx,i)
             end do
          end do
          close(1)
       else if(mesh%dim.eq.1) then  ! 1D
          write(filecube,'(a,i0,a)') 'evec',i,'.dat'
          open(unit=1,file=filecube,form='formatted',status='unknown')
          do j=1,mesh%N
             write(1,*) j*mesh%dx,V(j,i)
          end do
          close(1)
       else
          print *,' STOP in main(): dimension=',mesh%dim,' not yet implemented!'
          stop
       end if
    end do
  end subroutine save_wavefunction

end module IO
