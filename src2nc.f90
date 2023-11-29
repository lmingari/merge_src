program src2nc
    use netcdf
    implicit none
    !
    integer, parameter  :: nens      = 128
    integer, parameter  :: OUT_NT    = 800
    integer, parameter  :: OUT_DT    = 200 !seconds
    integer, parameter  :: OUT_NZ    = 800
    real(8), parameter  :: OUT_DZ    = 60  !meters
    !
    character(len=128)  :: fname_src = "millennium.src"
    character(len=128)  :: fname_nc  = "millennium.src.nc"
!    character(len=128)  :: fname_wf  = "weights.dat"
    !
    character(len=*), parameter :: form1 = '(i7,1x,i7)'
    character(len=*), parameter :: form2 = '(i7,1x,i7)'
    character(len=*), parameter :: form3 = '(e16.9)'
    character(len=*), parameter :: form4 = '(2(1x,f11.6),2x,f9.0,2x,100(e16.9,1x))'
    !
    character(len=5)    :: path
    integer             :: ncid,istat,i
    integer,allocatable :: t(:)
    real(8),allocatable :: z(:)
    real(8),allocatable :: mfr(:)    ! Mass flow rate (kg/s)
    real(8),allocatable :: src(:,:)  ! Linear source (kg/s/m)
    real(8),allocatable :: w(:)      ! Weight factors
    !
    write(*,*) "Program to convert an ensemble of src files into netCDF"
    !
!    if(command_argument_count().ge.1) then
!        call get_command_argument(1,fname_wf)
!    else
!        write(*,*) "Using default file name for file of weight factors"
!    endif
    !
    allocate(src(OUT_NZ,OUT_NT))
    allocate(mfr(       OUT_NT))
    allocate(t  (       OUT_NT))
    allocate(z  (OUT_NZ       ))
    !
    t  = [(OUT_DT*i,i=0,OUT_NT-1)]
    z  = [(OUT_DZ*i,i=0,OUT_NZ-1)]
    !
    ! Reading weight factors
!    write(*,*) "Reading file: ", trim(fname_wf)
!    istat = get_weights(fname_wf,w)
    !
    write(*,*) "Creating file: ", trim(fname_nc)
    ncid = create_nc(fname_nc)
    !
    do i=1,nens
        write(path,'(I0.4,a1)') i,'/'
        !
        src(:,:) = 0.0
        mfr(:)   = 0.0
        !
        write(*,*) "Reading member: ", i
        istat = get_src(path//fname_src,t,mfr,src)
        if(istat.ne.0) exit
        call save_nc(ncid)
        if(istat.ne.0) exit
    end do
    !
    istat = free_memory(ncid)
    if(istat.ne.0) stop "Stopping with errors"
    write(*,*) "Succesfully run"
    !
    contains
        !
        function create_nc(fname) result(ncid)
            implicit none
            !
            character(len=*), intent(in) :: fname
            integer                      :: ncid
            !
            integer :: istat
            integer, dimension(3) :: dimids,varids
            integer               :: varid
            !
            istat = nf90_create(fname, NF90_CLOBBER, ncid)
            !
            istat = nf90_def_dim(ncid, "lev",  OUT_NZ, dimids(1))
            istat = nf90_def_dim(ncid, "time", OUT_NT, dimids(2))
            istat = nf90_def_dim(ncid, "ens",  nens,   dimids(3))
            !
            istat = nf90_def_var(ncid, "lev",  NF90_FLOAT, dimids(1), varids(1))
            istat = nf90_def_var(ncid, "time", NF90_INT,   dimids(2), varids(2))
!            istat = nf90_def_var(ncid, "w",    NF90_FLOAT, dimids(3), varids(3))
            istat = nf90_def_var(ncid, "src",  NF90_FLOAT, dimids,    varid)
            istat = nf90_def_var(ncid, "mfr",  NF90_FLOAT, [dimids(2),dimids(3)], varid)
            !
            istat = nf90_enddef(ncid)
            !
            istat = nf90_put_var(ncid, varids(1), z)
            istat = nf90_put_var(ncid, varids(2), t)
!            istat = nf90_put_var(ncid, varids(3), w)
            !
            return
            !
        end function create_nc
        !
        subroutine save_nc(ncid)
            implicit none
            !
            integer, intent(in) :: ncid
            integer :: istat,varid
            integer :: iens = 0
            !
            iens = iens + 1
            !
            istat = nf90_inq_varid(ncid, "src", varid)
            istat = nf90_put_var(ncid, varid, src, &
                start = [1,1,iens],                &
                count = [OUT_NZ,OUT_NT,1]          )
            !
            istat = nf90_inq_varid(ncid, "mfr", varid)
            istat = nf90_put_var(ncid, varid, mfr, &
                start = [1,iens],                  &
                count = [OUT_NT,1]                 )
            !
        end subroutine save_nc
        !
        function free_memory(ncid) result(istat)
            implicit none
            !
            integer, intent(in) :: ncid
            integer             :: istat
            !
            if(allocated(t))   deallocate(t)
            if(allocated(z))   deallocate(z)
            if(allocated(src)) deallocate(src)
            if(allocated(mfr)) deallocate(mfr)
            !
            istat = nf90_close(ncid)
            !
        end function free_memory
        !
        function get_src(fname,t,mfr,src) result(istat)
            implicit none
            !
            character(len=*), intent(in)        :: fname
            integer, allocatable, intent(in)    :: t(:)
            real(8), allocatable, intent(inout) :: mfr(:)
            real(8), allocatable, intent(inout) :: src(:,:)
            integer                             :: istat
            !
            integer :: nin
            integer :: ip,it,iz
            integer :: t1,t2
            integer :: np,nb
            real(8) :: m
            real(8) :: x,y
            real(8) :: dt,dz
            integer :: it1,it2,tmin,tmax
            integer :: iz1,iz2,zmin,zmax
            real(8) :: z1,z2
            logical :: first
            real(8),allocatable :: lev(:)
            real(8),allocatable :: work2d(:,:)
            !
            !
            ! Open files
            !
            open(newunit=nin,      &
                 file=fname,       &
                 iostat=istat,     &
                 status='old',     &
                 err=300           )
            if(istat.ne.0) then
                write(*,*) "Error reading: ", trim(fname)
                return
            end if
            !
            dt  = real(OUT_DT,8)
            dz  = real(OUT_DZ,8)
            !
            ! Read files
            !
            first = .True.
            do
                read(nin,form1,iostat=istat) t1,t2
                if(istat.ne.0) then
                    istat=0
                    exit
                endif
                read(nin,form2) np,nb
                read(nin,form3) m
                !
                it1 = t1/OUT_DT+1
                it2 = t2/OUT_DT+1
                !
                do it=it1,it2
                    tmin = max(t(it)  ,t1)
                    tmax = min(t(it+1),t2)
                    mfr(it) = mfr(it) + (tmax-tmin)/dt*m
                end do
                !
                if(first) then
                    allocate(lev(np))
                    allocate(work2d(nb,np))
                    first = .False.
                end if
                !
                if(size(lev).ne.np) then
                    write(*,*) "Inconsistent number of levels"
                    istat = -1
                    return
                end if
                !
                do ip=1,np
                    read(nin,form4) x,y,lev(ip),work2d(:,ip)
                end do 
                !
                if(it2.ge.OUT_NT) then
                    write(*,*) "Max. time steps was exceeded"
                    istat = -1
                    return
                end if
                !
                do ip=1,np-1 
                    z1  = lev(ip)
                    z2  = lev(ip+1)
                    m   = sum(work2d(:,ip))/(z2-z1)
                    iz1 = int(z1/OUT_DZ)+1
                    iz2 = int(z2/OUT_DZ)+1
                    do it=it1,it2
                        tmin = max(t(it)  ,t1)
                        tmax = min(t(it+1),t2)
                        do iz=iz1,iz2
                            zmin = max(z(iz),  z1)
                            zmax = min(z(iz+1),z2)
                            src(iz,it) = src(iz,it) + (tmax-tmin)/dt*(zmax-zmin)/dz*m
                        end do
                    end do
                end do
            end do
            !
            close(nin)

            return
300         write(*,*) "Error fatal"
            istat=-1
            return
            !
        end function get_src
        !
        function get_weights(fname,w) result(istat)
            implicit none
            !
            character(len=*), intent(in)        :: fname
            real(8), allocatable, intent(inout) :: w(:)
            integer                             :: istat
            !
            integer :: nin,n,i
            !
            if(allocated(w)) deallocate(w)
            !
            open(newunit=nin,      &
                 file=fname,       &
                 iostat=istat,     &
                 status='old',     &
                 err=300           )

            n=0
            do
                read(nin,*,iostat=i)
                if(i.ne.0) exit
                n=n+1
            end do
            !
            if(n.gt.0) then
                write(*,*) "Reading ",n," weight factors. Expected: ", nens
                rewind(nin)
                allocate(w(n))
                do i=1,n
                    read(nin,*) w(i)
                end do
            end if
            !
200         close(nin)
            return
            !
300         write(*,*) "Error fatal"
            istat=-1
            return
        end function get_weights
        !
end program src2nc
