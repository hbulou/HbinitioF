module param_mod
  use global
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              read_param()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_param(param)
    implicit none
    character (len=1024)::line

    type(t_param)::param
    integer::lline,eqidx
    double precision, parameter :: pi=3.1415927

    param%ieof=0
    param%loopmax=1000
    param%restart=.FALSE.
    param%init_wf=.TRUE.
    param%nvecmin=20
    param%nvecmax=41
    param%Nx=30
    param%noccstate=1
    param%nvec_to_cvg=20
    param%ETA=1.0e-3
    param%dim=1
    param%box_width=pi/sqrt(2.0)
    param%Iperturb=1.0
    param%sigma=1.0
    open(unit=1,file='inp',form='formatted')
    do while(.not.(is_iostat_end(param%ieof)))
       read(1,*,iostat=param%ieof) line
       lline=len_trim(line)
       eqidx=index(line,"=")
       print *,'###',eqidx,lline
       if(line(1:eqidx-1).eq."restart") then
          if(line(eqidx+1:lline).eq.'.TRUE.') then
             param%restart=.TRUE.
          else
             param%restart=.FALSE.
          end if
       end if
       if(line(1:eqidx-1).eq."init_wf") then
          if(line(eqidx+1:lline).eq.'.TRUE.') then
             param%init_wf=.TRUE.
          else
             param%init_wf=.FALSE.
          end if
       end if
       if(line(1:eqidx-1).eq."loopmax") then
          read(line(eqidx+1:lline),*) param%loopmax
       end if
       if(line(1:eqidx-1).eq."nvecmin") then
          read(line(eqidx+1:lline),*) param%nvecmin
       end if
       if(line(1:eqidx-1).eq."nvecmax") then
          read(line(eqidx+1:lline),*) param%nvecmax
       end if
       if(line(1:eqidx-1).eq."Nx") then
          read(line(eqidx+1:lline),*) param%nx
       end if
       if(line(1:eqidx-1).eq."noccstate") then
          read(line(eqidx+1:lline),*) param%noccstate
       end if
       if(line(1:eqidx-1).eq."ETA") then
          read(line(eqidx+1:lline),*) param%ETA
       end if
       if(line(1:eqidx-1).eq."nvec_to_cvg") then
          read(line(eqidx+1:lline),*) param%nvec_to_cvg
       end if
       if(line(1:eqidx-1).eq."box_width") then
          read(line(eqidx+1:lline),*) param%box_width
       end if
       if(line(1:eqidx-1).eq."dimension") then
          read(line(eqidx+1:lline),*) param%dim
       end if
       if(line(1:eqidx-1).eq."sigma") then
          read(line(eqidx+1:lline),*) param%sigma
       end if
       if(line(1:eqidx-1).eq."Iperturb") then
          read(line(eqidx+1:lline),*) param%Iperturb
       end if
       line=''
    end do
    close(1)


    print *,'#restart=',param%restart
    print *,'#init_wf=',param%init_wf
    print *,'#loopmax=',param%loopmax
    print *,'#nvecmin=',param%nvecmin
    print *,'#nvecmax=',param%nvecmax
    print *,'#ETA=',param%ETA
    print *,'#nvec_to_cvg=',param%nvec_to_cvg
    print *,'#box_width=',param%box_width
    print *,'#Nx=',param%nx
    print *,'#noccstate=',param%noccstate
    print *,'#dh=',param%box_width/(param%Nx+1)
    print *,'#Dimension of the mesh=',param%dim
    print *,'#Magnitude of the perturbation=',param%Iperturb
    print *,'#Spread of the perturbation=',param%sigma
    
  end subroutine read_param
end module param_mod
