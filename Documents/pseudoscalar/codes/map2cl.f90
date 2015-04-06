PROGRAM map2Cl
    USE healpix_types
    USE alm_tools,      ONLY: map2alm, alm2map, alm2cl
    USE fitstools,      ONLY: getsize_fits,input_map,write_bintab,getnumext_fits
    USE head_fits,      ONLY: add_card
    USE misc_utils,     ONLY: assert,assert_alloc,fatal_error,wall_clock_time,string
    USE pix_tools,      ONLY: nside2npix,npix2nside,convert_nest2ring
    USE extension,      ONLY: getEnvironment
    USE extension,      ONLY:  getArgument, nArguments
    USE paramfile_io,   ONLY: paramfile_handle, parse_init,parse_int,parse_string,parse_real,parse_double,concatnl
    USE udgrade_nr
    
    IMPLICIT NONE
     
    INTEGER(I4B) :: nlmax,nmmax,n_rings,n_rings_read,nmw8
    INTEGER(KIND=I4B) :: ordering,IQU,read1,read2,ncl,part
    COMPLEX(KIND=DPC), DIMENSION(:,:,:), ALLOCATABLE :: alm1,alm2
    INTEGER(I4B) :: l_loop,ipnest,ipring,imap,extno,n_ext,column,l_conv1,l_conv2 
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: dummy,map_in,map_out,cl_auto
    REAL(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: w8ring_TQU,dw8
    REAL(KIND=DP), DIMENSION(1:2) :: cos_theta_cut
    INTEGER(KIND=I4B) :: polar_fits,nmaps,nside_in,npix_in,npix_out,status
    INTEGER(KIND=I4B) :: type, nside_out
    REAL(KIND=DP) :: fmissval=-1.6375e-30
    REAL(DP),parameter::PII = (2.0_dp)*asin(1.0_dp),d2r = PII/180.0, r2d = 180.0/PII

    CHARACTER(LEN=80), DIMENSION(1:120) :: header, header_in,header1
    CHARACTER(LEN=filenamelen)          :: TMap,QUMap
    CHARACTER(LEN=FILENAMELEN)          :: infile_w8
    LOGICAL(lgt)                        :: do_double
    CHARACTER(LEN=filenamelen)          :: parafile = ''
    CHARACTER(LEN=*), PARAMETER :: code = "Cl calculator"
    type(paramfile_handle) :: handle 
!=============================================================================
    if (nArguments() == 0) then
        parafile=''
    else 
        if (nArguments() /= 1) then
            write(*,*) "(' Usage: Cl calculator [parameter file name]')"
            stop 1
        end if
        call getArgument(1,parafile)
    end if

    handle = parse_init(parafile)

    nside_out = parse_int(handle,'nside_out')
    
    IQU = parse_int(handle,'IQU')
  
    TMap = parse_string(handle,'TMap',filestatus = 'old')
    
    QUMap = parse_string(handle,'QUMap',filestatus = 'old')

    infile_w8=parse_string(handle,'infile_w8')
    
    nlmax = 2*nside_out

    nmmax = nlmax

    npix_out = nside2npix(nside_out)
    
    if (IQU .eq. 1) then
        ncl=1
    else if (IQU .eq. 3) then
        ncl=4
    end if
  
!=========================================================================
    allocate(map_out(0:npix_out-1,1:IQU),stat = status)
    call assert_alloc(status, code, "map_out")
    
    do part=1,3
        if ((part .eq. 1) .or. (part .eq. 2)) then
            extno = 0
        else if (part .eq. 3) then
            extno = 1
        end if 
        if (part .eq. 1) then 
            n_ext = getnumext_fits(TMap)
            npix_in = getsize_fits(TMap, nmaps=nmaps, ordering=ordering, nside=nside_in, type=type,&
            polarisation = polar_fits)
            print*, extno
        else if ((part .eq. 2) .or. (part .eq. 3)) then
            n_ext = getnumext_fits(QUMap)
            npix_in = getsize_fits(QUMap, nmaps=nmaps, ordering=ordering, nside=nside_in, type=type,&
            polarisation = polar_fits)
            print*, extno
        end if
    
        if (npix_in .ne. nside2npix(nside_in)) then
            write(*,*)"      "//code//"> ERROR: Npix. Aborting."
            stop
        end if
    
        if ((nside_in .lt. 0) .or. (ordering .eq. 0)) then
            write(*,*)"      "//code//"> ERROR: Nside and/or pixelization scheme unknown. Aborting."
            stop
        end if

        write(*,*) 'Extensions   =', n_ext
        write(*,*) 'nmaps        =', nmaps
        write(*,*) 'Pixelization =', ordering
        write(*,*) 'nside        =', nside_in
        write(*,*) 'file type    =', type
        write(*,*) 'polarisation =', polar_fits 

        column = 1
        
        allocate(dummy(0:npix_in-1,1:column),stat = status)
        call assert_alloc(status, code, "dummy")
   
        if (part .eq. 1) then 
            call input_map(TMap, dummy, npix_in, column, &
            fmissval=fmissval, header=header_in,extno=extno)
        else if ((part .eq. 2) .or. (part .eq. 3))  then
            call input_map(QUMap, dummy, npix_in, column, &
            fmissval=fmissval, header=header_in,extno=extno)
        end if
    
        do l_conv1=0,npix_in-1
            do l_conv2=1,column
                dummy(l_conv1,l_conv2)=dummy(l_conv1,l_conv2)*(1.0E3_dp)
            end do
        end do
        
        allocate(map_in(0:npix_in-1,1:column),stat = status)
        call assert_alloc(status, code, "map_in")
    
        if (nside_in .ne. nside_out) then
            if (ordering .eq. 1) then
                call udgrade_ring(dummy,nside_in,map_in,nside_out,fmissval)
            else if(ordering .eq. 2) then
                call udgrade_nest(dummy,nside_in,map_in,nside_out,fmissval)
            end if
        else if (nside_in .eq. nside_out) then
            map_in=dummy
        end if
    
        if (part .eq. 1) then
            do read1=0,npix_out-1
                map_out(read1,1)=map_in(read1,1)
            end do
        else if ((IQU .eq. 3) .and. (part .eq. 2)) then
            do read1=0,npix_out-1
                map_out(read1,2)=map_in(read1,1)
            end do
        else if ((IQU .eq. 3) .and. (part .eq. 3)) then
            do read1=0,npix_out-1
                map_out(read1,3)=map_in(read1,1)
            end do
        end if
        deallocate (map_in,dummy)
    end do

    if (ordering .eq. 2) then 
        write(*,*)"      "//code//"> Converting Nest to Ring "
        call convert_nest2ring(nside_out,map_out)
    end if

       
!---------------------------get alm's of the input map--------------------------            
   
    allocate(alm1(1:IQU, 0:nlmax, 0:nmmax),stat = status)
    call assert_alloc(status, code, "alm1")
     
    allocate(alm2(1:IQU, 0:nlmax, 0:nmmax),stat = status)
    call assert_alloc(status, code, "alm2")
  
    alm1 = CMPLX(0.0e0_dp,0.0e0_dp)
    alm2 = CMPLX(0.0e0_dp,0.0e0_dp)
    
    allocate(w8ring_TQU(1:2*nside_out,1:IQU),stat = status)
    call assert_alloc(status,code,"w8ring_TQU")

    allocate(dw8(1:2*nside_out,1:IQU),stat = status)
    call assert_alloc(status,code,"dw8")

    allocate(cl_auto(0:nlmax,1:ncl))
    call assert_alloc(status, code,"cl_auto")
    dw8=0.0_dp

    if (trim(infile_w8)/="") then
        n_rings = 2*nside_out
        n_rings_read = getsize_fits(infile_w8, nmaps=nmw8)

        if (n_rings_read /= n_rings) then
            write(*,*) " "
            write(*,*) "wrong ring weight file:"//trim(infile_w8)
            call fatal_error(code)
        endif
        
        nmw8 = min(nmw8, IQU)

        write(*,*) "      "//code//"> Inputting Quadrature ring weights "
        call input_map(infile_w8, dw8, n_rings, nmw8, fmissval=0.0_dp)
    end if
    
    w8ring_TQU = 1.0_dp + dw8

    deallocate(dw8)

    cos_theta_cut= 0.0e0_dp
    
    write(*,*)"      "//code//"> Calculating alm from map "
    call map2alm (nside_out,nlmax,nmmax,map_out,alm1,cos_theta_cut,w8ring_TQU)
    write(*,*)"      "//code//"> Calculating Cls from alm "
    call alm2cl(nlmax, nmmax, alm1, cl_auto)

    open(unit=2,file='Cls_Planck2.dat' ,status='unknown')
    
    do l_loop=5,500
        write(2,*) l_loop, (real(l_loop,dp)*(real(l_loop,dp)+1.0_dp)/(2*PII))*cl_auto(l_loop,1),&
        (real(l_loop,dp)*(real(l_loop,dp)+1.0_dp)/(2*PII))*cl_auto(l_loop,2),&
        (real(l_loop,dp)*(real(l_loop,dp)+1.0_dp)/(2*PII))*cl_auto(l_loop,3)
    end do
    
    deallocate(map_out,alm1,cl_auto)
    close(2)
END PROGRAM map2Cl
