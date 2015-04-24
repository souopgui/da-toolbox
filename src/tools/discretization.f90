!> \file discretization.f90
!! Discretization tools for rectangular domains
!! @author Innocent Souopgui
!!
!<

!> Discretization module, define data structures and subroutines for discretization.
!!
!<
module discretization
    use general_constant
    use general_tools
    use debug_tools
implicit none

    !> \brief User defined type for discretization parameters
    !! \details
    !<
    type discretization_param
        integer :: ndim !< number of physical dimension of the problem, 0-3
        !x direction
        real(kind=cp) :: xmin!< minimum coordinate in the x direction of the physical domain
        real(kind=cp) :: xsize!< size of physical domain the in the x direction
        integer :: nx!< number of discretization points in the x direction of the physical domain
        real(kind=cp) :: dx !< space step in the x direction of the physical domain, uniform discretization
        real(kind=cp), dimension(:), pointer :: x => null()!< Coordinates of discretization points in the x direction of the physical domain
        logical :: x_allocated = .false. !< says is the array x has been allocated

        !y direction
        real(kind=cp) :: ymin!< minimum coordinate in the y direction of the physical domain
        real(kind=cp) :: ysize!< size of physical domain the in the y direction
        integer :: ny!< number of discretization points in the y direction of the physical domain
        real(kind=cp) :: dy !< space step in the y direction of the physical domain, uniform discretization
        real(kind=cp), dimension(:), pointer :: y => null()!< Coordinates of discretization points in the y direction of the physical domain
        logical :: y_allocated = .false. !< say is the array y has been allocated

        !z direction
        real(kind=cp) :: zmin!< minimum coordinate in the z direction of the physical domain
        real(kind=cp) :: zsize!< size of physical domain the in the z direction
        integer :: nz!< number of discretization points in the z direction of the physical domain
        real(kind=cp) :: dz !< space step in the z direction of the physical domain, uniform discretization
        real(kind=cp), dimension(:), pointer :: z => null()!< Coordinates of discretization points in the z direction of the physical domain
        logical :: z_allocated = .false. !< says is the array z has been allocated

        !time
        real(kind=cp) :: tmin!< Minimum time
        real(kind=cp) :: tsize!< Size of the simulation window (in s)
        integer :: nt!< Total number of time steps in the simulation
        real(kind=cp) :: dt!< time step, uniform discretization

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                  User defined parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! The subdomain A was introduced as the response region in
        ! the sensitivity analysis problem
        ! we expect the user to provide min and size, we will keep min and max
        real(cp) :: A_xmin, A_xmax, A_xsize !< lower boundary and size of A subdomain in x
        real(cp) :: A_ymin, A_ymax, A_ysize !< lower boundary and size of A subdomain in y
        real(cp) :: A_zmin, A_zmax, A_zsize !< lower boundary and size of A subdomain in z

        ! ENO / WENO approximation parameters
        !< @brief left cells count in the stencil for interpolation in
        !< the ENO approximation
        integer :: eno_r
        !< @brief cells count in the stencil for interpolation in ENO approximation
        integer :: eno_k
        integer :: lb, ub !< @brief bounds of the part of array used in computation
        !< @brief extended bounds taken into account extention for
        !< ENO/WENO approximations.
        !>
        integer :: ext_lb, ext_ub
        !< @brief reduced bounds; define the part of the state variable that
        !< goes to the control vector, useful when the boundaries are not
        !< considered as part of the control vector, otherwise, identical
        !< to lb and ub .
        !>
        integer :: red_lb, red_ub
        real(cp), dimension(:), pointer :: &
            eno_ce =>null(),&
            eno_cw =>null(),&
            weno_de=>null(),&
            weno_dw=>null()
    end type discretization_param

    !> discretization parameters, global variable
    type(discretization_param), save :: tg_dp

contains

    !> @brief Load discretization parameters from namelist
    !! @param[in,out] td_dp discretization parameter, the [in] attribute is used to free dynamic memory if necessary
    !! @param[in] pFName name of the namelist file
    !! @details
    !! Most of the time, boundary conditions are provided by and external
    !! process and do not need to be controlled. State variable always included
    !! the boundaries, however, they must be excluded from the control variable
    !! unless they are controlled.
    !! we will possibly add this two more parameters:
    !! param[in] bred reduction in the boundaries
    !! param[in] bext extension in the boundaries for ENO/WENO schemes
    !! @a bred gives the number of grid points
    !! to exclude from the state variable:
    !! @a bred is a two column array that defines the number of
    !! grid points/cells used for the boundary conditions.
    !!
    !! When using schemes like ENO/WENO, the state variable used in computation
    !! is extended by many more grid/cells to match the stencil at the boundaries.
    !! In this case too, the computing variable is an extended state variable.
    !! The extension cells need to be removed to get the actual state variable.
    !! @a bext is is a two column array that defines the number of extension
    !! grid points/cells.
    !!
    !! The number of rows of of @a bred and @a bext is equal to the number
    !! of dimensions of the problem. Those arrays are set up as follows:
    !!   - the row i corresponds to the dimension i (1 for x, 2 for y and 3 for z)
    !!   - the column 1 corresponds to lower bounds
    !!   - the column 2 corresponds to upper bounds
    !! For example,
    !! bred(1,1) is the number of grid points used for the lower boundary in x
    !! bred(1,2) is the number of grid points used for the upper boundary in x
    !! bred(2,1) is the number of grid points used for the lower boundary in y
    !! bred(2,2) is the number of grid points used for the upper boundary in y
    !! bred(3,1) is the number of grid points used for the lower boundary in z
    !! bred(3,2) is the number of grid points used for the upper boundary in z
    !! The same logic applies to @a bext
    !!
    !! There is a strong difference between the two: @a bred is used to get the
    !! control vector from the state vector; @a bext is used to get the actual
    !! state vector from the extended state vector.
    !!
    !<
    subroutine load_discretization(td_dp, pFName)
        integer, parameter:: ip_numnam = 68
        type(discretization_param), intent(in out) :: td_dp
        character(len=*), intent(in) :: pFName
        !integer, dimension(:,:), optional, intent(in) :: bred, bext
        !local variables
        character(len=ip_fnl) :: domainType, timeType, tmp_str
        real(kind=cp), dimension(3) :: coordMin, delta_space, domainSize
        real(kind=cp) :: tMin, delta_t, simulTime
        integer, dimension(3) :: nPoint
        integer :: nDimension, nTime
        logical :: positive_size

        namelist/NAM_discretization/&
        nDimension ,&
        domainType ,&
        timeType   ,&
        coordMin   ,&
        tMin       ,&
        nPoint     ,&
        nTime      ,&
        delta_space,&
        delta_t    ,&
        domainSize ,&
        simulTime

        coordMin = 0.0
        delta_space = 0.0
        domainSize = 0.0
        nPoint = -1
        OPEN(ip_numnam,FILE=pFName,FORM='FORMATTED',STATUS='OLD')
        !If one is going to read many blocs in the same file,
        !it is recommended to rewind
        REWIND(ip_numnam)
        READ(ip_numnam, NAM_discretization)!reading the block
        CLOSE(ip_numnam)

        select case (toLower(timeType))
        case('d')!discrete information
            simulTime = nTime*delta_t
        case('c')!continuous information
            delta_t = simulTime/nTime
        case default
            call debug('domain typ should be:')
            call debug('  - "c" for continuous or "d" for discrete, check the namelist')
            call stop_program('In load_discretization, bad type of domain type, see above ')
        end select
        select case (toLower(domainType))
        case('d')!discrete information
            domainSize = nPoint*delta_space
        case('c')!continuous information
            delta_space = domainSize/(nPoint-1)
        case default
            call stop_program('Bad type of domain type, should be "c" for continuous or "d" for discrete, check the namelist')
        end select
        positive_size = (count(domainSize>0.0)==nDimension)

        if( (.not.positive_size).or.(count(nPoint>0)/=nDimension) )then
            tmp_str = 'The size of delta_space/domainSize and nPoint should correspond to nDimension.'
            tmp_str = trim(tmp_str)//' Negative values or zeroes are not allowed! check the namelist'
            call stop_program(tmp_str)
        end if

        call set_discretization(td_dp, coordMin(1:nDimension)&
        , delta_space(1:nDimension), nPoint(1:nDimension)&
        , tMin, delta_t, nTime)
    end subroutine load_discretization

    !> Set the discretization parameters
    !! @param [in,out] td_dp discretization data structure
    !! @param [in] cMin array giving the minimum coordinate in each dimension of the physical space.
    !! @param [in] delta  array giving the discretization step in each dimension of the physical space.
    !! @param [in] np  array giving the number of discretization points in each dimension of the physical space.
    !! @param [in] tMin starting time of the simulation
    !! @param [in] delta_t time step
    !! @param [in] nTime total number of time steps in the simulation
    !!
    !! \details \a coordmin(i), \a delta(i), \a np(i) correspond to the ith dimension of the physical space.
    !<
    subroutine set_discretization(td_dp, cMin, delta, np, tMin, delta_t, nTime)
        TYPE(discretization_param), INTENT(INOUT) :: td_dp
        REAL(KIND=cp), dimension(:), intent(in) :: cMin, delta
        REAL(KIND=cp), optional, intent(in) :: tMin, delta_t
        INTEGER, dimension(:), intent(in) :: np
        INTEGER, optional, intent(in) :: nTime
        !local variables
        INTEGER :: nDimension, t_nParam

        nDimension = size(cMin)
        if(&
            (nDimension/=size(delta))&
            .or.(nDimension/=size(np)) &
            .or.(nDimension>3)&
            )then
            call stop_program('cMin, delta and np should have the same size, that should be <= 3')
        end if

        td_dp%ndim = nDimension

        if(nDimension>=1)then
            call setCoordinates(td_dp, 'x', np(1), delta(1), cMin(1))
        end if
        if(nDimension>=2)then
            call setCoordinates(td_dp, 'y', np(2), delta(2), cMin(2))
        end if
        if(nDimension>=3)then
            call setCoordinates(td_dp, 'z', np(3), delta(3), cMin(3))
        end if

        t_nParam = count( (/present(tmin),present(delta_t),present(nTime)/) )
        select case(t_nParam)
        case(3)
            call setCoordinates(td_dp, 't', nTime, delta_t, tmin)
        case(0)
        case default
            !call stop_program('set_discretization must be call with exactly none or all of the 3 times parameters')
        end select

    end subroutine set_discretization

    !> @brief load a subdomain from the namelist
    !! @param[in,out] td_dp discretization parameter, the [in] attribute
    !!   is used to free dynamic memory if necessary
    !! @param[in] pFName name of the namelist file
    !! @param[in] dName name of the domain
    !! @details the coordinate of the subdomain are stored in the fields
    !! A_[x,y,z]min, A_[x,y,z]max and A_[x,y,z]size of td_dp
    !<
    subroutine load_subdomain( td_dp, pFName, dName )
        integer, parameter:: fId = 68
        type(discretization_param), intent(in out) :: td_dp
        character(len=*), intent(in) :: pFName, dName
        !local variables
        real(kind=cp), dimension(3) :: xmin, xsize
        character(len=ip_snl) :: domainName, ala_dName
        namelist/NAM_DOMAIN/&
            domainName  ,&
            xmin,&
            xsize
        !
        open(fId,file=pFName,form='FORMATTED',status='OLD')
        domainName = "{#}@"
        ala_dName = TRIM(dName)
        call uppercase(ala_dName)
        do while(domainName.ne.ala_dName)
            read(fId, NAM_DOMAIN)!reading the block
            call uppercase(domainName)
            call debug(domainName, "In load_subdomain, found varName", tag=dALLWAYS)
        end do
        close(fId)
        if(td_dp%ndim>=1)then
            td_dp%A_xmin = xmin(1)
            td_dp%A_xmax = xmin(1)+xsize(1)
            td_dp%A_xsize = xsize(1)
        end if
        !
        if(td_dp%ndim>=2)then
            td_dp%A_ymin = xmin(2)
            td_dp%A_ymax = xmin(2)+xsize(2)
            td_dp%A_ysize = xsize(2)
        end if
        !
        if(td_dp%ndim>=2)then
            td_dp%A_zmin = xmin(2)
            td_dp%A_zmax = xmin(2)+xsize(2)
            td_dp%A_zsize = xsize(2)
        end if
    end subroutine load_subdomain

    !> @brief load the eno parameters from the namelist
    !! @param[in,out] td_dp discretization parameter, the [in] attribute
    !!   is used to free dynamic memory if necessary
    !! @param[in] pFName name of the namelist file
    !<
    subroutine load_enoWeno(td_dp, pFName)
        integer, parameter:: ip_fid = 68
        type(discretization_param), intent(in out) :: td_dp
        character(len=*), intent(in) :: pFName
        !local variables
        integer :: il_eno_r, il_eno_k
        namelist/NAM_ENO/&
        il_eno_r, &
        il_eno_k

        open(ip_fid,file=pFName,form='formatted',status='old')
        !If one is going to read many blocs in the same file, it is recommended to rewind
        rewind(ip_fid)
        read(ip_fid, NAM_ENO)!reading the block
        close(ip_fid)

        td_dp%eno_r  = il_eno_r
        td_dp%eno_k  = il_eno_k
        td_dp%ext_lb = td_dp%lb - td_dp%eno_k!+1
        td_dp%ext_ub = td_dp%ub + td_dp%eno_k

        !eno multipliers for east boundaries
        if (associated(td_dp%eno_ce)) deallocate(td_dp%eno_ce)
        allocate( td_dp%eno_ce(0:td_dp%eno_k-1) )
        td_dp%eno_ce = eno_multipliers(td_dp%eno_r, td_dp%eno_k)

        !eno multipliers for west boundaries
        if (associated(td_dp%eno_cw)) deallocate(td_dp%eno_cw)
        allocate( td_dp%eno_cw(0:td_dp%eno_k-1) )
        td_dp%eno_cw = eno_multipliers(td_dp%eno_r-1, td_dp%eno_k)

        !weno parameter d
        call weno_compute_d(td_dp)
    end subroutine load_enoWeno

    !> Set the coordinates of the discretization points in the \a dir direction of the physical domain or in the time domain
    !! @param [in,out] td_dp discretization data structure, the [in] attribute is used to free dynamic memory if necessary
    !! @param [in] axis dimension where the coordinates are needed
    !! @param [in] n number of discrete coordinates along the axis @a axis
    !! @param [in] delta discretization step along the axis @a axis
    !! @param [in] minCoord minimum coordinate along the axis @a axis
    !!
    !! \details dir should be :
    !! - 'x' or 'X' for the x axis of the physical domain
    !! - 'y' or 'Y' for the Y axis of the physical domain
    !! - 'z' or 'Z' for the Z axis of the physical domain
    !! - 't' or 'T' for the T axis
    !<
    subroutine setCoordinates(td_dp, axis, n, delta, minCoord)
        TYPE(discretization_param), INTENT(INOUT) :: td_dp
        character(len=1), intent(in) :: axis
        INTEGER, INTENT(IN) :: n
        REAL(cp), INTENT(IN) :: delta
        REAL(cp), OPTIONAL, INTENT(IN) :: minCoord
        !local variables
        REAL(cp), DIMENSION(:), POINTER :: coords
        integer :: i

        if((n<0).or.(delta<0))then
        call stop_program('Bad discretization parameters in the dimension '//axis//'. Check the namelist')
        end if
        if(toLower(axis)/='t')then
        allocate(coords(n))
        do i = 1,n
            coords(i) = minCoord + real(i-1,cp)*delta
        end do
        end if
        select case(toLower(axis))
        case('x')
            if (td_dp%x_allocated) deallocate(td_dp%x)
            td_dp%x => coords
            td_dp%x_allocated = .TRUE.
            td_dp%xmin = minCoord
            td_dp%xSize= (n-1)*delta
            td_dp%nx = n
            td_dp%dx = delta
        case('y')
            if (td_dp%y_allocated) deallocate(td_dp%y)
            td_dp%y => coords
            td_dp%y_allocated = .TRUE.
            td_dp%ymin = minCoord
            td_dp%ySize= (n-1)*delta
            td_dp%ny = n
            td_dp%dy = delta
        case('z')
            if (td_dp%z_allocated) deallocate(td_dp%z)
            td_dp%z => coords
            td_dp%z_allocated = .TRUE.
            td_dp%zmin = minCoord
            td_dp%zSize= (n-1)*delta
            td_dp%nz = n
            td_dp%dz = delta
        case('t')
            td_dp%tmin = minCoord
            td_dp%tsize = n*delta
            td_dp%nt = n
            td_dp%dt = delta
        case default
            deallocate(coords)
            call stop_program("Bad dimension; can set coordinates only in the physical dimension x, y and z ")
        end select
    end subroutine setCoordinates

    !> @brief compute multipliers for approximation of boundaries values
    !> in ENO schemes
    !! @param [in] id_r left cells count in the stencil, relative to
    !!  the cell for which approximation is needed
    !! @param [in] id_k cells count in the stencil
    !! @details
    !! When calling this procedure, you know what you are doing.
    !! It is not harmful, but it was design only for 1D problem
    !! When solving problems in more than 1 dimension, make sure that it adapts
    !! to your case
    !<
    function eno_multipliers(id_r, id_k) result(rla_ce)
        integer, intent(in) :: id_r, id_k
        real(kind=cp), dimension(0:id_k-1) :: rla_ce

        select case (id_k)
            case(1)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/ 1. /), cp)
                    case(0)
                        rla_ce = real((/ 1. /), cp)
                end select
            case(2)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/  3./2, -1./2 /), cp)
                    case(0)
                        rla_ce = real((/  1./2,  1./2 /), cp)
                    case(1)
                        rla_ce = real((/ -1./2,  3./2 /), cp)
                end select
            case(3)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/ 11./6, -7./6,  1./3 /), cp)
                    case(0)
                        rla_ce = real((/  1./3,  5./6, -1./6 /), cp)
                    case(1)
                        rla_ce = real((/ -1./6,  5./6,  1./3 /), cp)
                    case(2)
                        rla_ce = real((/  1./3, -7./6, 11./6 /), cp)
                end select
            case(4)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/ 25./12, -23./12,  13./12,  -1./4 /), cp)
                    case(0)
                        rla_ce = real((/   1./4,  13./12,  -5./12,  1./12 /), cp)
                    case(1)
                        rla_ce = real((/ -1./12,   7./12,   7./12, -1./12 /), cp)
                    case(2)
                        rla_ce = real((/  1./12,  -5./12,  13./12,   1./4 /), cp)
                    case(3)
                        rla_ce = real((/  -1./4,  13./12, -23./12, 25./12 /), cp)
                end select
            case(5)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/ 137./60, -163./60, 137./60,  -21./20,    1./5 /), cp)
                    case(0)
                        rla_ce = real((/    1./5,   77./60, -43./60,   17./60,  -1./20 /), cp)
                    case(1)
                        rla_ce = real((/  -1./20,    9./20,  47./60,  -13./60,   1./30 /), cp)
                    case(2)
                        rla_ce = real((/   1./30,  -13./60,  47./60,    9./20,  -1./20 /), cp)
                    case(3)
                        rla_ce = real((/  -1./20,   17./60, -43./60,   77./60,    1./5 /), cp)
                    case(4)
                        rla_ce = real((/    1./5,  -21./20, 137./60, -163./60, 137./60 /), cp)
                end select
            case(6)
                select case (id_r)
                    case(-1)
                        rla_ce = real((/ 49./20, -71./20,   79./20, -163./60,  31./30,  -1./6 /), cp)
                    case(0)
                        rla_ce = real((/   1./6,  29./20,  -21./20,   37./60, -13./60,  1./30 /), cp)
                    case(1)
                        rla_ce = real((/ -1./30,  11./30,   19./20,  -23./60,   7./60, -1./60 /), cp)
                    case(2)
                        rla_ce = real((/  1./60,  -2./15,   37./60,   37./60,  -2./15,  1./60 /), cp)
                    case(3)
                        rla_ce = real((/ -1./60,   7./60,  -23./60,   19./20,  11./30, -1./30 /), cp)
                    case(4)
                        rla_ce = real((/  1./30, -13./60,   37./60,  -21./20,  29./20,   1./6 /), cp)
                    case(5)
                        rla_ce = real((/  -1./6,  31./30, -163./60,   79./20, -71./20, 49./20 /), cp)
                end select
            case(7)
                select case (id_r)
                    case(-1)
                    rla_ce = real((/ 363./140 -617./140, 853./140, -2341./420, 667./210, -43./42, 1./7 /), cp)
                    case(0)
                    rla_ce = real((/ 1./7, 223./140, -197./140, 153./140, -241./420, 37./210, -1./42 /), cp)
                    case(1)
                    rla_ce = real((/ -1./42, 13./42, 153./140, -241./420, 109./420, -31./420, 1./105 /), cp)
                    case(2)
                    rla_ce = real((/ 1./105, -19./210, 107./210, 319./420, -101./420, 5./84, -1./140 /), cp)
                    case(3)
                    rla_ce = real((/ -1./140, 5./84, -101./420, 319./420, 107./210, -19./210, 1./105 /), cp)
                    case(4)
                    rla_ce = real((/ 1./105, -31./420, 109./420, -241./420, 153./140, 13./42, -1./42 /), cp)
                    case(5)
                    rla_ce = real((/ -1./42, 37./210, -241./420, 153./140, -197./140, 223./140, 1./7 /), cp)
                    case(6)
                    rla_ce = real((/ 1./7, -43./42, 667./210, -2341./420, 853./140, -617./140, 363./140 /), cp)
                end select
        end select
    end function eno_multipliers


    !> @brief compute WENO paramters for 1D problems
    !! @param [in] td_tp discretization data structure
    !! @details
    !! When calling this procedure, you know what you are doing.
    !! It is not harmful, but it was design only for 1D problem
    !! When solving problems in more than 1 dimension, make sure that it adapts
    !! to your case
    !<
    subroutine weno_compute_d(td_tp)
        !see 2.54
        type(discretization_param), intent(inout) :: td_tp
        integer :: ibr

        if (associated(td_tp%weno_de)) deallocate(td_tp%weno_de)
        allocate( td_tp%weno_de(0:td_tp%eno_k-1) )
        select case(td_tp%eno_k)
            case (1)
                td_tp%weno_de = real( (/ 1. /), cp)
            case (2)
                td_tp%weno_de = real( (/ 2./3, 1./3 /), cp)
            case (3)
                td_tp%weno_de = real( (/ 3./10, 3./5, 1./10 /), cp)
        end select
        if (associated(td_tp%weno_dw)) deallocate(td_tp%weno_dw)
        allocate( td_tp%weno_dw(0:td_tp%eno_k-1) )
        !symmetric computation for weno_dw
        do ibr = 0, td_tp%eno_k-1
            td_tp%weno_dw(ibr) = td_tp%weno_dw(td_tp%eno_k - 1 - ibr)
        end do
    end subroutine weno_compute_d

    !> Print the discretization parameters
    !! @param [in] td_dp discretization data structure
    !<
    subroutine print_discretization_param(td_dp)
        type(discretization_param), intent(in) :: td_dp

        call debug('Printing discretization parameters ', tag=dALLWAYS)
        if(td_dp%ndim>=1)then
            call debug(td_dp%xmin, 'xmin', tag=dALLWAYS)
            call debug(td_dp%nx, 'nx', tag=dALLWAYS)
            call debug(td_dp%xSize, 'xSize', tag=dALLWAYS)
            call debug(td_dp%dx, 'delta_x', tag=dALLWAYS)
        end if
        if(td_dp%ndim>=2)then
            call debug(td_dp%ymin, 'ymin', tag=dALLWAYS)
            call debug(td_dp%ny, 'ny', tag=dALLWAYS)
            call debug(td_dp%ySize, 'ySize', tag=dALLWAYS)
            call debug(td_dp%dy, 'delta_y', tag=dALLWAYS)
        end if
        if(td_dp%ndim>=3)then
            call debug(td_dp%zmin, 'zmin', tag=dALLWAYS)
            call debug(td_dp%nz, 'nz', tag=dALLWAYS)
            call debug(td_dp%zSize, 'zSize', tag=dALLWAYS)
            call debug(td_dp%dz, 'delta_z', tag=dALLWAYS)
        end if

        call debug(td_dp%tmin, 'tmin', tag=dALLWAYS)
        call debug(td_dp%nt, 'nt', tag=dALLWAYS)
        call debug(td_dp%tSize, 'tSize', tag=dALLWAYS)
        call debug(td_dp%dt, 'delta_t', tag=dALLWAYS)
        call debug('-----------------------------------', tag=dALLWAYS)
    end subroutine print_discretization_param

end module discretization

