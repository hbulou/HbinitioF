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
    double precision :: data(:),val
!    integer :: idxmin,idxmax
    type(t_mesh) :: m
    character (len=1024) :: filename
    
    character(len=*),parameter :: FMT1='(I5,3F12.6)'
    integer :: i,j,k,nn,ifield

    print *,"# save_cube_3D > saving ",trim(filename)
    
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
             nn=m%ijk_to_idx(i,j,k)%n
             val=0.0
             if(m%node(nn)%active) val=data(nn)
             write(1,'(E13.5)',advance='no') val
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
!    print *,"# save_cube_3D > ok"
  end subroutine save_cube_3D
  ! --------------------------------------------------------------------------------------
  !
  !              save_wavefunction(param,mesh,V,molecule)
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_wavefunction(param,mesh,V,molecule)
    implicit none
    type(t_param)::param
    type(t_mesh)::mesh
    double precision::V(:,:)
    type(t_molecule)::molecule
    
    integer :: i,j,k
    character (len=1024) :: filecube
    
    do i=1,param%nvec_to_cvg
       call norm(mesh,V(:,i))
       call dcopy(molecule%mesh%nactive,V(:,i),1,molecule%wf%wfc(:,i),1)
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
          do j=1,mesh%nactive
             write(1,*) j*mesh%dx,V(j,i)
          end do
          close(1)
       else
          print *,' STOP in main(): dimension=',mesh%dim,' not yet implemented!'
          stop
       end if
    end do

  end subroutine save_wavefunction
  ! --------------------------------------------------------------------------------------
  !
  !                             save_config()
  !
  ! subroutine to save the configuration of the calculation in order to restart it
  ! later if necessary
  ! --------------------------------------------------------------------------------------
  subroutine save_config(V,m,nvecmin,param)
    implicit none
    type(t_mesh)::m
    type(t_param)::param
    double precision :: V(:,:)
    integer::nvecmin,i,j
    character (len=1024)::filename

    write(filename,'(a,a)') param%prefix(:len_trim(param%prefix)),'/evectors.dat'
    print *,"# save_config > saving ",trim(filename)
    open(unit=1,file=filename,form='formatted',status='unknown')
    do i=1,m%nactive
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
       do i=1,m%nactive
          read(1,*) (V(i,j),j=1,nvecmin)
       end do
       close(1)
    else
       print *,"### ERROR ### evectors.dat doesn't exist"
       stop
    end if
  end subroutine read_config

end module IO
