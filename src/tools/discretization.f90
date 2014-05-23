!> \file discretization.f90
!! Discretization tools for rectangular domains
!! @author Innocent Souopgui
!!
!<

!> @brief Discretization module, define data structures and subroutines for finite difference-type discretization.
!! @todo
!! Check the behavior of modules that uses the discretization module. A correction has been made on May 23rd 2014 by Innocent.
!! The modification should affect function that computes something as a function of the physical coordinates given by the discretization parameter
!! In the normal case, it should not affects computations that does not used physical coordinates directly.
!! Now, there is a discretization point at each end of the domain in each direction. The number of discretization points if 1 unit greater than the number of cells.
!! Previously, there was no discretization point at the min coordinate, The number of discretization points corresponded to the number of cells.
!!
!<
MODULE discretization
    USE general_constant
    use general_tools
    USE debug_tools
IMPLICIT NONE

    !> @brief User defined type for discretization parameters
    !! @details
    !<
    TYPE discretization_param
      integer :: ndim !< number of physical dimension of the problem, 0-3
      !x direction
      REAL(KIND=cp) :: xmin!< Minimum coordinate in the x direction of the physical domain
      REAL(KIND=cp) :: xsize!< Size of physical domain the in the x direction
      INTEGER :: nx!< number of discretization points in the x direction of the physical domain
      REAL(KIND=cp) :: dx !< space step in the x direction of the physical domain, uniform discretization
      REAL(KIND=cp), DIMENSION(:), POINTER :: x => NULL()!< Coordinates of discretization points in the x direction of the physical domain
      logical :: x_allocated = .FALSE. !< says is the array x has been allocated

      !y direction
      REAL(KIND=cp) :: ymin!< Minimum coordinate in the y direction of the physical domain
      REAL(KIND=cp) :: ysize!< Size of physical domain the in the y direction
      INTEGER :: ny!< number of discretization points in the y direction of the physical domain
      REAL(KIND=cp) :: dy !< space step in the y direction of the physical domain, uniform discretization
      REAL(KIND=cp), DIMENSION(:), POINTER :: y => NULL()!< Coordinates of discretization points in the y direction of the physical domain
      logical :: y_allocated = .FALSE. !< say is the array y has been allocated

      !z direction
      REAL(KIND=cp) :: zmin!< Minimum coordinate in the z direction of the physical domain
      REAL(KIND=cp) :: zsize!< Size of physical domain the in the z direction
      INTEGER :: nz!< number of discretization points in the z direction of the physical domain
      REAL(KIND=cp) :: dz !< space step in the z direction of the physical domain, uniform discretization
      REAL(KIND=cp), DIMENSION(:), POINTER :: z => NULL()!< Coordinates of discretization points in the z direction of the physical domain
      logical :: z_allocated = .FALSE. !< says is the array z has been allocated

      !time
      REAL(KIND=cp) :: tmin!< Minimum time
      REAL(KIND=cp) :: tsize!< Size of the simulation window (in s)
      INTEGER :: nt!< Total number of time steps in the simulation
      REAL(KIND=cp) :: dt!< time step, uniform discretization
    END TYPE discretization_param

    !> discretization parameters, global variable
    type(discretization_param), save :: tg_dp

CONTAINS

  !> \brief Load discretization parameters from namelist
  !! \param[in,out] td_dp discretization parameter, the [in] attribute is used to free dynamic memory if necessary
  !! \param[in] pFName name of the namelist file
  !<
  SUBROUTINE load_discretization(td_dp, pFName)
    INTEGER, PARAMETER:: ip_numnam = 68
    TYPE(discretization_param), INTENT(IN OUT) :: td_dp
    CHARACTER(len=*), INTENT(IN) :: pFName
    CHARACTER(len=ip_fnl) :: domainType, timeType, tmp_str
    REAL(KIND=cp), dimension(3) :: coordMin, delta_space, domainSize
    REAL(KIND=cp) :: tMin, delta_t, simulTime
    INTEGER, dimension(3) :: nPoint
    INTEGER :: nDimension, nTime
    logical :: positive_size

    NAMELIST/NAM_discretization/&
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
    nPoint = -1
    OPEN(ip_numnam,FILE=pFName,FORM='FORMATTED',STATUS='OLD')
    !If one is going to read many blocs in the same file, it is recommended to rewind
    REWIND(ip_numnam)
    READ(ip_numnam, NAM_discretization)!reading the block
    CLOSE(ip_numnam)

    select case (toLower(timeType))
      case('d')!discrete information
        simulTime = nTime*delta_t
      case('c')!continuous information
        delta_t = simulTime/nTime
      case default
        call stop_program('Bad type of domain type, should be "c" for continuous or "d" for discrete, check the namelist')
    end select
    select case (toLower(domainType))
      case('d')!discrete information
        domainSize = (nPoint-1)*delta_space
        positive_size = (count(delta_space>epsilon(1.0))==nDimension)
      case('c')!continuous information
        positive_size = (count(domainSize>epsilon(1.0))==nDimension)
        delta_space = domainSize/(nPoint-1)
      case default
        call stop_program('Bad type of domain type, should be "c" for continuous or "d" for discrete, check the namelist')
    end select


    if( (.not.positive_size).or.(count(nPoint>0)/=nDimension)&
    )then
      tmp_str = 'The size of delta_space/domainSize and nPoint should correspond to nDimension.'
      tmp_str = trim(tmp_str)//' Negative values or zeroes are not allowed! check the namelist'
      call stop_program(tmp_str)
    end if

    call set_discretization(td_dp, coordMin(1:nDimension)&
      , delta_space(1:nDimension), nPoint(1:nDimension)&
      , tMin, delta_t, nTime)
  END SUBROUTINE load_discretization

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

  !> Set the coordinates of the discretization points in the \a dir direction of the physical domain or in the time domain
  !! @param [in,out] td_dp discretization data structure, the [in] attribute is used to free dynamic memory if necessary
  !! @param [in] dir dimension where the coordinates are needed
  !! @param [in] n number of discrete coordinates in the dimension dir
  !! @param [in] delta discretization step in the dimension dir
  !! @param [in] minCoord minimum coordinate in the dimension dir
  !!
  !! \details dir should be :
  !! - 'x' or 'X' for the x dimension of the physical domain
  !! - 'y' or 'Y' for the Y dimension of the physical domain
  !! - 'z' or 'Z' for the Z dimension of the physical domain
  !! - 't' or 'T' for the Z dimension of the physical domain
  !<
  subroutine setCoordinates(td_dp, dir, n, delta, minCoord)
    TYPE(discretization_param), INTENT(INOUT) :: td_dp
    character(len=1), intent(in) :: dir
    INTEGER, INTENT(IN) :: n
    REAL(cp), INTENT(IN) :: delta
    REAL(cp), OPTIONAL, INTENT(IN) :: minCoord
    !local variables
    REAL(cp), DIMENSION(:), POINTER :: coords
    integer :: i

    if((n<0).or.(delta<0))then
      call stop_program('Bad discretization parameters in the dimension '//dir//'. Check the namelist')
    end if
    if(toLower(dir)/='t')then
      allocate(coords(n))
      do i = 1,n
        coords(i) = minCoord + real(i-1,cp)*delta
      end do
    end if
    select case(toLower(dir))
      case('x')
        if (td_dp%x_allocated) deallocate(td_dp%x)
        td_dp%x => coords
        td_dp%x_allocated = .TRUE.
        td_dp%xmin = minCoord
        td_dp%xSize= n*delta
        td_dp%nx = n
        td_dp%dx = delta
      case('y')
        if (td_dp%y_allocated) deallocate(td_dp%y)
        td_dp%y => coords
        td_dp%y_allocated = .TRUE.
        td_dp%ymin = minCoord
        td_dp%ySize= n*delta
        td_dp%ny = n
        td_dp%dy = delta
      case('z')
        if (td_dp%z_allocated) deallocate(td_dp%z)
        td_dp%z => coords
        td_dp%z_allocated = .TRUE.
        td_dp%zmin = minCoord
        td_dp%zSize= n*delta
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

  !> Print the discretization parameters
  !! @param [in] td_dp discretization data structure
  !<
  SUBROUTINE print_discretization_param(td_dp)
      TYPE(discretization_param), INTENT(IN) :: td_dp

      call debug('Printing discretization parameters ', tag=dAllways)
      if(td_dp%ndim>=1)then
        call debug(td_dp%xmin, 'xmin', tag=dAllways)
        call debug(td_dp%nx, 'nx', tag=dAllways)
        call debug(td_dp%xSize, 'xSize', tag=dAllways)
        call debug(td_dp%dx, 'delta_x', tag=dAllways)
      end if
      if(td_dp%ndim>=2)then
        call debug(td_dp%ymin, 'ymin', tag=dAllways)
        call debug(td_dp%ny, 'ny', tag=dAllways)
        call debug(td_dp%ySize, 'ySize', tag=dAllways)
        call debug(td_dp%dy, 'delta_y', tag=dAllways)
      end if
      if(td_dp%ndim>=3)then
        call debug(td_dp%zmin, 'zmin', tag=dAllways)
        call debug(td_dp%nz, 'nz', tag=dAllways)
        call debug(td_dp%zSize, 'zSize', tag=dAllways)
        call debug(td_dp%dz, 'delta_z', tag=dAllways)
      end if

      call debug(td_dp%tmin, 'tmin', tag=dAllways)
      call debug(td_dp%nt, 'nt', tag=dAllways)
      call debug(td_dp%tSize, 'tSize', tag=dAllways)
      call debug(td_dp%dt, 'delta_t', tag=dAllways)
      call debug('-----------------------------------', tag=dAllways)
    END SUBROUTINE print_discretization_param

END MODULE discretization
