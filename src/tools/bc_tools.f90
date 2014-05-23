!!Special module for boundary conditions
!!@author Innocent Souopgui
MODULE bc_tools
  USE general_constant
  USE general_tools
  USE debug_tools
  USE general_discretization
  IMPLICIT NONE


  !> \brief User defined type for boundary conditions parameters
  !! Defines parameters that are used to define or manage the boundary conditions
  !<
  TYPE bc_param
    CHARACTER(LEN=ip_snl) ::&
      aa_title,& !> title of the output file
      aa_hist ,& !> history string for the output file
      aa_prefix  !> prefix of the BC filename
    INTEGER :: i_ndim = -1 !> number of physical dimensions for the problem
    REAL(dp):: r_dt   = -1.0_cp
    INTEGER :: i_ndata = -1 !> \brief Number of boundary data when there is no grouping by axis
    !> Number of data points in each axis direction (used when there is a grouping by axis)
    !! In the case of a grouping by axis, the number of data on a boundary face orthogonal to the x axis is given by ny*nz; nx*nz for a face orthogonal to y axis and nx*ny for a face orthogonal to z axis. With nx=ia_nxyz(1), ny=ia_nxyz(2) and nz=ia_nxyz(3)
    !<
    INTEGER, DIMENSION(:), POINTER :: ia_nxyz
!     INTEGER ::&
!       i_ndata = -1 !> \brief Number of boundary data when there is no grouping by axis
!       i_nx = 1 !> Number of boundary data on the boundaries in the x direction (used when there is a grouping by axis)
!       i_ny = 1 !> Number of boundary data on the boundaries in the y direction (used when there is a grouping by axis)
!       i_nz = 1 !> Number of boundary data on the boundaries in the x direction (used when there is a grouping by axis)
!     !In the case of a grouping by axis, the number of data on a boundary face orthogonal to the x axis is given by i_ny*i_nz; i_nx*i_nz for a face orthogonal to y axis and i_nx*i_ny for a face orthogonal to z axis.

    !> \brief grouping status of BC data.
    !! Are data grouped by axis limits? For example, bc for x_max(east), xmin (west), ymax(north), ymin(south), ...
    !! Grouping by axis limit facilitates the management for rectangular domains and regular discretization.
    !! It can be impractical for non rectangular domain or irregular discretization.
    !<
    LOGICAL :: l_by_axis
    !> \brief integer coordinates status
    !! Integer coordinates refer to the indices in a multi-D array, it is useful for rectangular domains and regular discretization.
    !! used only when there is no grouping by axis
    !<
    LOGICAL :: l_icoord
    !> \brief real coordinates status
    !! The real coordinates are the coordinates in the physical space
    !! used only when there is no grouping by axis
    !<
    LOGICAL :: l_rcoord
    !> \brief ordinal coordinates status
    !! ordinal coordinates refer to the indices in a 1-D array, it is useful when the system state is store in a 1D array
    !! used only when there is no grouping by axis
    !<
    LOGICAL :: l_ocoord
  END TYPE bc_param

  !> \brief User defined type for boundary data at runtime
  !<
  TYPE bc_structure
    INTEGER :: i_ts
    REAL(dp):: r_date !> date associated with the boundary data
    !> \brief grouping status of BC data.
    !! Are data grouped by axis limits? For example, bc for x_max(east), xmin (west), ymax(north), ymin(south), ...
    !! Grouping by axis limit facilitates the management for rectangular domains and regular discretization.
    !! It can be impractical for non rectangular domain or irregular discretization.
    !<
    LOGICAL :: l_by_axis
    !> \brief boundary data when there is no grouping by axis
    !! useful for unstructured grid or non rectangular domain
    !<
    REAL(dp), DIMENSION(:), POINTER :: ra_data => NULL() !> boundary data for vectorized state vector
    !If the boundary data are grouped by axis, the variables below are used
    !for the corners, prefered storage is the axis that comes first. any boundary data belonging to xmin or xmax face is saved as xmin data, respectively xmax data, even if the same data belongs to ymin or ymax face.
    REAL(dp), DIMENSION(:,:), POINTER ::&
      ra_xmin_bc => NULL(),&!> data for xmin(West face) bounday
      ra_xmax_bc => NULL(),&!> data for xmax(East face) bounday
      ra_ymin_bc => NULL(),&!> data for ymin(South face) bounday
      ra_ymax_bc => NULL(),&!> data for ymax(North face) bounday
      ra_zmin_bc => NULL(),&!> data for zmin(Bottom face) bounday
      ra_zmax_bc => NULL()  !> data for zmax(Top face) bounday

!     !> \brief real coordinates status
!     !! this logical field says if real-type coordinates are presents.  Real coordinates are the coordinates in the physical domain
!     !<
!     LOGICAL :: l_rcoord = .FALSE.
!     !> \brief integer coordinates status
!     !! this logical field says if integer-type coordinates are presents. Integer coordinates are the indices in a multi-D array
!     !<
!     LOGICAL :: l_icoord = .FALSE.
!     !> \brief ordinal index coordinates status
!     !! this logical field says if ordinal index coordinates are presents. ordinal index coordinates are the indices in a 1D array
!     !<
!     LOGICAL :: l_ocoord = .FALSE.
!     !> \brief Allocation status of the integer coordinates array.
!     !! this logical field says if there is a proper array allocated for the given structure variable. For a time dependent problem, the location of boundary data can be the same for every time step. In such a case, only one array of coordinates needs to be allocated for all time step in order to save memory space. The variable l_icoord_allocated is use to manage memory deallocation in this case.
!     !<
!     LOGICAL :: l_icoord_allocated = .FALSE.
!     !> \brief Allocation status of the real coordinates array. See l_icoord_allocated for details
!     LOGICAL :: l_rcoord_allocated = .FALSE.
!     !> \brief Allocation status of the ordinal index coordinates array. See l_icoord_allocated for details
!     LOGICAL :: l_ocoord_allocated = .FALSE.
!     !> \brief boundary data coordinates, indices in the discretization grid
!     !! This array gives the indices of the boundary data in the computation grid.
!     !! used only when there is no grouping by axis
!     !<
!     INTEGER, DIMENSION(:,:), POINTER ::  ia_icoord => NULL()
!     !> \brief boundary data coordinates, real coordinates in the computation domain
!     !! used only when there is no grouping by axis
!     !<
!     REAL(dp), DIMENSION(:,:), POINTER :: ra_rcoord => NULL()
!     !> \brief boundary data coordinates, ordinal index in the vectorized system state
!     !! used only when there is no grouping by axis
!     !<
!     INTEGER, DIMENSION(:), POINTER ::  ia_ocoord => NULL()
  END TYPE bc_structure

  INTERFACE extract_bc
    !MODULE PROCEDURE extract_bc_1d
    MODULE PROCEDURE extract_bc_2d
    !MODULE PROCEDURE extract_bc_3d
  END INTERFACE extract_bc

  INTERFACE set_bc
    !MODULE PROCEDURE set_bc_1d
    MODULE PROCEDURE set_bc_2d
    !MODULE PROCEDURE set_bc_3d
  END INTERFACE set_bc

  INTERFACE set_bcAdj
    !MODULE PROCEDURE set_bcAdj_1d
    MODULE PROCEDURE set_bcAdj_2d
    !MODULE PROCEDURE set_bcAdj_3d
  END INTERFACE set_bcAdj

CONTAINS

  !> \brief load ctl parameters from a namelist file
  !! \param[in,out] td_bcp data structure to be loaded
  !! \param[in] ada_namelist name of the namelist file
  !! \param[in] ada_varName name of the variable associated with td_bcp
  !! \param[in] ndim number of dimensions of the physical space
  !! \param[in] nxyz (optional) number of boundary data in each axis direction, the default value 1 is used for the direction where no information is provided. The number of boundary data includes the corner; Basically, it is equal to 2 plus the number of data in the corresponding dimension of the computational domain; Boundary data correspond to ghost points.
  !! \param[in] dt time step of the model evolution
  !! The program automatically stops in IO error if the variable is not found in the namelist file
  !<
  SUBROUTINE load_bcp(td_bcp, ada_namelist, ada_varName, ndim, nxyz, dt)
    INTEGER, PARAMETER:: ip_numnam = 68
    TYPE(bc_param), INTENT(IN OUT) :: td_bcp
    CHARACTER(len=*), INTENT(IN) :: ada_namelist, ada_varName
    INTEGER, INTENT(IN) :: ndim
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: nxyz
    REAL(cp), OPTIONAL, INTENT(IN) :: dt
    !local variables
    !REAL(cp) :: rl_sigmaB
    CHARACTER(len=ip_snl) :: varName, title, history, ala_varName
    LOGICAL :: by_axis, icoord, rcoord, ocoord
    INTEGER :: ibi
    NAMELIST/NAM_BC_PARAM/&
      varName,&
      title  ,&
      history,&
      by_axis,&
      icoord ,&
      rcoord ,&
      ocoord

    OPEN(ip_numnam,FILE=ada_namelist,FORM='FORMATTED',STATUS='OLD')
    varName = "{#}@"
    ala_varName = TRIM(ada_varName)
    CALL uppercase(ala_varName)
    DO WHILE(varName.NE.ala_varName)
      READ(ip_numnam, NAM_BC_PARAM)!reading the block
      CALL uppercase(varName)
      CALL debug(varName, "In load_bcp, found varName = ", tag=dALLWAYS)
    END DO
    CLOSE(ip_numnam)
    td_bcp%aa_title  = title
    td_bcp%aa_hist   = history
    td_bcp%l_by_axis = by_axis
    td_bcp%l_icoord  = icoord
    td_bcp%l_rcoord  = rcoord
    td_bcp%l_ocoord  = ocoord
    td_bcp%i_ndim    = ndim
    ALLOCATE( td_bcp%ia_nxyz(ndim) )
    IF( PRESENT(nxyz) )THEN
      DO ibi=1,SIZE(nxyz)
        td_bcp%ia_nxyz(ibi) = nxyz(ibi)
      END DO
      DO ibi=SIZE(nxyz)+1, ndim
        td_bcp%ia_nxyz(ibi) = 1
      END DO
    ELSE
      td_bcp%ia_nxyz = 1
    END IF
    IF( PRESENT(dt) )THEN
      td_bcp%r_dt = dt
    ELSE
      td_bcp%r_dt = 1.0_cp
    END IF

    CALL lowercase(varName)
    td_bcp%aa_prefix  = TRIM( varName )
  END SUBROUTINE load_bcp

  !> \brief print ctl parameters
  !! \param[in] td_bcp data structure for ctl parameters
  !!
  !<
  SUBROUTINE print_bcp(td_bcp)
    TYPE(bc_param), INTENT(IN) :: td_bcp

    CALL debug('', 'printing bc_param---------------------------------', tag=dALLWAYS)
    CALL debug(td_bcp%aa_title , '  aa_title  = ', tag=dALLWAYS)
    CALL debug(td_bcp%aa_hist  , '  aa_hist   = ', tag=dALLWAYS)
    CALL debug(td_bcp%aa_prefix, '  aa_prefix = ', tag=dALLWAYS)
    CALL debug(td_bcp%l_by_axis, '  l_by_axis = ', tag=dALLWAYS)
    CALL debug(td_bcp%l_icoord , '  l_icoord  = ', tag=dALLWAYS)
    CALL debug(td_bcp%l_rcoord , '  l_rcoord  = ', tag=dALLWAYS)
    CALL debug(td_bcp%l_ocoord , '  l_ocoord  = ', tag=dALLWAYS)
    CALL debug(td_bcp%i_ndim   , '  i_ndim    = ', tag=dALLWAYS)
    IF(td_bcp%l_by_axis)THEN
      CALL debug(td_bcp%ia_nxyz  , '  ia_nxyz   = ', tag=dALLWAYS)
    ELSE
      CALL debug(td_bcp%i_ndata  , '  i_ndata   = ', tag=dALLWAYS)
    END IF
    CALL debug('', '.......................................................', tag=dALLWAYS)
  END SUBROUTINE print_bcp

  !> \brief Initializes boundary condition structure
  !! \param[in out] td_bc boundary condition structure to be initialized
  !! \param[in] td_bcp parameters to use for the initialization
  !<
  SUBROUTINE init_bcs(td_bc, td_bcp)
    TYPE(bc_structure), INTENT(IN OUT) :: td_bc
    TYPE(bc_param), INTENT(IN) :: td_bcp
    !local variables
    INTEGER :: nx, ny, nz

    IF( .NOT.td_bcp%l_by_axis ) CALL stop_program('In init_bcs, only grouping by axis is supported for now')
    td_bc%l_by_axis= td_bcp%l_by_axis
    IF( td_bcp%i_ndim>=1 ) THEN
      nx = td_bcp%ia_nxyz(1)
    ELSE
      nx = 1
    END IF
    IF( td_bcp%i_ndim>=2 ) THEN
      ny = td_bcp%ia_nxyz(2)
    ELSE
      ny = 1
    END IF
    IF( td_bcp%i_ndim>=3 ) THEN
      nz = td_bcp%ia_nxyz(3)
    ELSE
      nz = 1
    END IF

    IF( td_bc%l_by_axis )THEN
      IF( td_bcp%i_ndim>3 )CALL stop_program("In init_bcs, can not manage more than 3 dimensional problems")
      ALLOCATE( td_bc%ra_xmin_bc( ny, nz ), td_bc%ra_xmax_bc( ny, nz ) )
      IF( td_bcp%i_ndim>=2 )ALLOCATE( td_bc%ra_ymin_bc( nx-2, nz   ), td_bc%ra_ymax_bc( nx-2, nz   ) )
      IF( td_bcp%i_ndim==3 )ALLOCATE( td_bc%ra_zmin_bc( nx-2, ny-2 ), td_bc%ra_zmax_bc( nx-2, ny-2 ) )
    ELSE
      CALL stop_program("In init_bcs, only the grouping by axis is implemented for now")
      !ALLOCATE ( td_bc%ra_data(td_bcp%i_ndata) )
    END IF

  END SUBROUTINE init_bcs

  SUBROUTINE print_bcs(td_bc)
    TYPE(bc_structure), INTENT(IN OUT) :: td_bc

    CALL debug('', 'Printing bc_structure---------------------------------', tag=dALLWAYS)
    IF(td_bc%l_by_axis)THEN
     CALL debug(SHAPE(td_bc%ra_xmin_bc)  , '  SHAPE(ra_xmin_bc)   = ', tag=dALLWAYS)
     CALL debug(SHAPE(td_bc%ra_xmax_bc)  , '  SHAPE(ra_xmax_bc)   = ', tag=dALLWAYS)
     CALL debug(SHAPE(td_bc%ra_ymin_bc)  , '  SHAPE(ra_xmin_bc)   = ', tag=dALLWAYS)
     CALL debug(SHAPE(td_bc%ra_ymax_bc)  , '  SHAPE(ra_xmax_bc)   = ', tag=dALLWAYS)
    ELSE
      CALL debug(SHAPE(td_bc%ra_data)  , '  SHAPE(ra_data)   = ', tag=dALLWAYS)
    END IF
    CALL debug('', '.......................................................', tag=dALLWAYS)
  END SUBROUTINE print_bcs

  SUBROUTINE reset_bcs(td_bc)
    TYPE(bc_structure), INTENT(IN OUT) :: td_bc
    !local variables

    IF( td_bc%l_by_axis )THEN
      IF( ASSOCIATED(td_bc%ra_xmin_bc) )THEN
        DEALLOCATE ( td_bc%ra_xmin_bc, td_bc%ra_xmax_bc )
        td_bc%ra_xmin_bc => NULL()
        td_bc%ra_xmax_bc => NULL()
      END IF
      IF( ASSOCIATED(td_bc%ra_ymin_bc) )THEN
        DEALLOCATE ( td_bc%ra_ymin_bc , td_bc%ra_ymax_bc )
        td_bc%ra_ymin_bc => NULL()
        td_bc%ra_ymax_bc => NULL()
      END IF
      IF( ASSOCIATED(td_bc%ra_zmin_bc) )THEN
        DEALLOCATE ( td_bc%ra_zmin_bc , td_bc%ra_zmax_bc )
        td_bc%ra_zmin_bc => NULL()
        td_bc%ra_zmax_bc => NULL()
      END IF
    ELSE
      CALL stop_program("In init_bcs, only the grouping by axis is implemented for now")
      !ALLOCATE ( td_bc%ra_data(td_bcp%i_ndata) )
    END IF

  END SUBROUTINE reset_bcs


  !> \brief extract boundary condition from a state variable
  !! \param[in out] td_bc data structure for the boundary condition
  !! \param[in] rda_state 2D array (state variable from which the extraction is performed)
  !! \param[in] td_rs (optional) data structure describing the subdomain
  !! \param[in] ts (optional) time state associated with the model state
  !! \param[in] nghost (optional) number of ghost data at the boundaries (each boundary) of rda_state
  !! td_rs can described either the full domain (if BC needed for the full domain) or the subdomain of interest
  !! if td_rs is not provided, the full domain is considered
  !! if nghost is not provided, it is assumed that there is no ghost cell
  !! This subroutine assumes that data structures are initialized with correct data and size.
  !<
  SUBROUTINE extract_bc_2d(td_bc, rda_state, td_rs, ts, nghost)
    TYPE(bc_structure), INTENT(IN OUT)      :: td_bc
    REAL(cp), DIMENSION(:,:), INTENT(IN)    :: rda_state
    TYPE(rectangular_subdomain), OPTIONAL, INTENT(IN) :: td_rs
    INTEGER, OPTIONAL, INTENT(IN) :: ts, nghost
    !local variables
    INTEGER :: imin, imax, jmin, jmax, i_idx, j_idx, il_nghost

    IF( PRESENT(nghost) )THEN
      il_nghost = nghost
    ELSE
      il_nghost = 0
    END IF
    i_idx = 1
    j_idx = 2
    IF( PRESENT(td_rs) )THEN
      imin = 1+td_rs%ia_shift(i_idx)+il_nghost!base 1 index in fortran
      imax = imin + td_rs%ia_size(i_idx) - 1
      jmin = 1+td_rs%ia_shift(j_idx)+il_nghost!base 1 index in fortran
      jmax = jmin + td_rs%ia_size(j_idx) - 1
      IF( (imax>UBOUND(rda_state,i_idx)-il_nghost).OR.(jmax>UBOUND(rda_state,j_idx)-il_nghost) )&
        CALL stop_program("In extract_bc_2d, subdomain exeeds the size of the entire domain")
    ELSE
      imin = il_nghost + 1
      imax = UBOUND(rda_state,i_idx)-il_nghost
      jmin = il_nghost + 1
      jmax = UBOUND(rda_state,j_idx)-il_nghost
    END IF
    !imin:imax, jmin:jmax give the range of the computational domain, boundary are set on the ghost points

    IF( td_bc%l_by_axis )THEN
      td_bc%ra_xmin_bc(:,1) = rda_state(imin-1, jmin-1:jmax+1)
      td_bc%ra_xmax_bc(:,1) = rda_state(imax+1, jmin-1:jmax+1)
      td_bc%ra_ymin_bc(:,1) = rda_state(imin:imax, jmin-1)
      td_bc%ra_ymax_bc(:,1) = rda_state(imin:imax, jmax+1)
    ELSE
      CALL stop_program("In extract_bc_2d, only the grouping by axis is implemented for now")
    END IF

    IF( PRESENT(ts) )THEN
      td_bc%i_ts = ts
    ELSE
      td_bc%i_ts = -1
    END IF
  END SUBROUTINE extract_bc_2d


  !> \brief set boundary condition on a state variable
  !! \param[in out] td_bc data structure for the boundary condition
  !! \param[in] rda_state data structure of the state
  !! \param[in] nghost (optional) number of ghost data at the boundaries (each boundary) of rda_state
  !! This subroutine assumes that data structures are initialized with correct data and size
  !! if nghost is not provided, it is assumed that there is no ghost cell
  !<
  SUBROUTINE set_bc_2d(td_bc, rda_state, nghost)
    TYPE(bc_structure), INTENT(IN)      :: td_bc
    REAL(cp), DIMENSION(:,:), INTENT(IN OUT)    :: rda_state
    INTEGER, OPTIONAL, INTENT(IN) :: nghost
    !local variables
    INTEGER :: imin, imax, jmin, jmax, i_idx, j_idx, il_nghost

    IF( PRESENT(nghost) )THEN
      il_nghost = nghost
    ELSE
      il_nghost = 0
    END IF
    i_idx = 1
    j_idx = 2
    imin = il_nghost + 1
    imax = UBOUND(rda_state,i_idx)- il_nghost
    jmin = il_nghost + 1
    jmax = UBOUND(rda_state,j_idx)- il_nghost
    !imin:imax, jmin:jmax give the range of the computational domain, boundary are set on the internal ghost points

    IF( td_bc%l_by_axis )THEN
      rda_state(imin-1, jmin-1:jmax+1) = td_bc%ra_xmin_bc(:,1)
      rda_state(imax+1, jmin-1:jmax+1) = td_bc%ra_xmax_bc(:,1)
      rda_state(imin:imax, jmin-1)     = td_bc%ra_ymin_bc(:,1)
      rda_state(imin:imax, jmax+1)     = td_bc%ra_ymax_bc(:,1)
    ELSE
      CALL stop_program("In set_bc_2d, only the grouping by axis is implemented for now")
    END IF
  END SUBROUTINE set_bc_2d


  !> \brief Adjoint of set_bc_2d
  !! \param[in] rda_statead data structure of the state
  !! This subroutine assumes that data structures are initialized with correct data and size
  !<
  SUBROUTINE set_bcAdj_2d(td_bc, rda_statead, nghost)
    TYPE(bc_structure), INTENT(IN)      :: td_bc
    REAL(cp), DIMENSION(:,:), INTENT(IN OUT) :: rda_statead
    INTEGER, OPTIONAL, INTENT(IN) :: nghost
    !local variables
    INTEGER :: imin, imax, jmin, jmax, i_idx, j_idx, il_nghost

    IF( PRESENT(nghost) )THEN
      il_nghost = nghost
    ELSE
      il_nghost = 0
    END IF
    i_idx = 1
    j_idx = 2
    imin = il_nghost + 1
    imax = UBOUND(rda_statead,i_idx)- il_nghost
    jmin = il_nghost + 1
    jmax = UBOUND(rda_statead,j_idx)- il_nghost
    !imin:imax, jmin:jmax give the range of the computational domain, boundary are set on the ghost points

    IF( td_bc%l_by_axis )THEN
      rda_statead(imin:imax, jmax+1)     = 0.0_cp
      rda_statead(imin:imax, jmin-1)     = 0.0_cp
      rda_statead(imax+1, jmin-1:jmax+1) = 0.0_cp
      rda_statead(imin-1, jmin-1:jmax+1) = 0.0_cp
    ELSE
      CALL stop_program("In set_bcAdj_2d, only the grouping by axis is implemented for now")
    END IF
  END SUBROUTINE set_bcAdj_2d

END MODULE bc_tools
