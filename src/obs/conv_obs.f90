!> \file conv_obs.f90
!! User defined types and subroutines/functions used sample and process observations
!!
!!
!<


!> Module defining data structures and general subroutines for observations in a data assimilation problem.
!!
!! \todo rewrite \a conv_obsgap_elt to make use of subsample
!<
module conv_obs
    !use balaise_constant
    use debug_tools
    use general_tools
    use nctools
implicit none
    !> @brief User defined type for observations at runtime
    !! PRIVATE (removed for compatibility with old compilers): i_max_nobs,
    !!  i_ndim, i_nobs, l_icoord_allocated, l_rcoord_allocated
    !! @a ra_sigma and @a ra_obsgap will be remove in the future.
    !!  They have been replaced
    !!  by ra_std (standard deviation) and ra_inn (for innovation).
    !<
    type obs_structure
        character(len=ip_snl) :: obs_fName, ogap_fName
        character(len=ip_ccl) :: varName !< Name of the variable associated with the obs data structure
        !> @brief maximum number of observations, usefull when the totalnumber of observations
        !! to be used is less than the total number of observation in the file
        !<
        integer :: i_max_nobs = 0
        integer :: i_ndim = 0 !< number of dimensions for the problem
        integer :: i_nobs   = 0 !< \brief total number of observations
        integer  :: i_ts
        real(dp) :: r_date,& !< date of observation
        r_cost,&!< cost associated with this obs
        r_costb !< adjoint variable associated with r_cost
        !> @brief real coordinates
        !! this logical field says if real-type coordinates are presents.
        !!  Real coordinates are the coordinates in the physical domain
        !<
        logical :: l_rcoord = .false.
        !> @brief integer coordinates
        !! this logical field says if integer-type coordinates are presents.
        !! Integer coordinates are the indices in a vector. They are usefull
        !! when observations correspond to a subsampling of the state variable
        !<
        logical :: l_icoord = .false.
        !> @brief Allocation status of the integer coordinates array.
        !! this logical field says if there is a proper array allocated for the given
        !!   obs_structure variable. For a time dependent problem, the location of
        !!   observation can be the same for every observation time. In such case,
        !!   only one array of coordinates can be allocated for all observation time
        !!   in order to save memory space. The variable l_icoord_allocated is use
        !!   to manage memory deallocation in this case.
        !<
        logical :: l_icoord_allocated = .false.
        !> @brief Allocation status of the real coordinates array. See l_icoord_allocated for details
        logical :: l_rcoord_allocated = .false.
        real(dp), dimension(:), pointer :: ra_obs    => null()!< observation data
        real(dp), dimension(:), pointer :: ra_std    => null()!< standard deviation of observation error
        real(dp), dimension(:), pointer :: ra_Rm1    => null()!< inverse covariance (diag matrix)
        real(dp), dimension(:), pointer :: ra_inn    => null()!< innovation, difference between the observation and the model output
        real(dp), dimension(:), pointer :: ra_sigma  => null()!< duplicate of @a ra_std, used for compatibility
        real(dp), dimension(:), pointer :: ra_obsgap => null()!< obsgap, duplicate of @a ra_inn for compatibility

        !wavelet related variables
        integer :: i_obs_level = 0!< observation level, this variable changes the interpretation of integer coordinates
        integer, dimension(3) :: ia_M_vector = (/0, 0, 0/)!<wavelet grid at the coarsest scale (level =1)
        integer, dimension(3) :: ia_nobs = (/1, 1, 1/)!<wavelet grid at the coarsest scale (level =1)
        real(dp), dimension(3) :: ra_dx !<Step of regularly spaced observations
        !> @brief Observation indices in the discretization grid
        !! This array gives the indices of the observations in the full grid at
        !!  observation level. It is to be differentiated from the indices in the
        !!  adapted grid that must be computed at each time
        !<
        integer, dimension(:,:), pointer :: ia_icoord => NULL()
        !> @brief Observation fractional indices. This is the substitute of integer
        !>  coordinates when Observations are not located exactly at discretization
        !>  grid points. They are mostly used for interpolation
        !<
        real(dp), dimension(:,:), pointer :: ra_fidx => NULL()
        !> @brief Observation coordinates, it represents
        !>  - the x-y-z in a Cartesian coordinate system;
        !>  - the radius and angle in a polar coordinate system;
        !>  - the radius, azimuthal angle and polar angle in a spherical coordinate system
        !>  - the longitude, latitude and elevation in a geographic coordinate system
        !>  - etc.
        !<
        real(dp), dimension(:,:), pointer :: ra_rcoord => NULL()
    end type obs_structure

    !> \brief user define type for conventional observation parameters
    !! For now, we suppose that observation are regularly spaced
    !! @todo Changes:
    !! - ia_shifts to ia_lshifts (left shifts)
    !! - ia_max to ia_rshifts (right shifts)
    !! as ia_shifts, ia_lshifts(i) gives the minimum index in the ith physical direction, counting from the right.
    !! ia_rshifts(i) gives the maximum index in the ith physical direction, counting from the right.
    !! Meaning for the physical direction i: the first observation point is at
    !! index ia_lshifts(i), and the last one does not exceed ia_rshifts(i). No
    !! observation point before ia_lshifts(i), and no observation point after
    !! ia_rshifts(i). All observation points are in the range [ia_lshifts(i), ia_rshifts(i)]
    !! Those changes induces some other changes in the subroutines that generates the indexes.
    !<
    TYPE convObs_param
        CHARACTER(LEN=ip_snl) ::&
        aa_title,& !< title of the output file
        aa_hist    !< history string for the output file
        INTEGER :: i_ndim = 0 !< number of physical dimensions of the problem
        INTEGER, DIMENSION(:), POINTER::&
        ia_shifts => NULL(),& !< shifts along the physical space, ia_shifts(i) gives the min index of the observation grid point in the ith dimension
        ia_steps  => NULL(),& !< steps between obs along the physical space, ia_steps(i) gives the number of grid points between two successive observations along the ith dimension
        ia_max    => NULL()   !< max along the physical space, ia_max(i) gives the max index of the observation grid point in the ith dimension
        INTEGER :: &
        i_tshift,& !< shift along the time coordinate, this is the index of the time step associated with the first observation
        i_tstep ,& !< step along the time coordinate, number of time steps between two successive observations
        i_tmax  !< time step of the last observation
        REAL(cp) :: r_sigmaR !< Standard deviation of observation error
        LOGICAL  :: l_sampleObs = .FALSE. !< if True, observation are sample.

        CHARACTER(LEN=ip_snl) :: aa_prefix !< prefix for the file name, this can be the name of the observed variable
        REAL(cp) :: r_dt
        INTEGER, DIMENSION(:), POINTER :: ia_samplingTS=> NULL() !< sampling time steps, arays of ordered time steps at which observations are sampled.
    END TYPE convObs_param

    !> \brief routine to set the observation data
    !<
    interface set_obs
        module procedure set_obs_with_icoord
        module procedure set_obs_with_fidx
        module procedure set_obs_data
    end interface

    !> brief interface for obsgap subroutines for a single variable
    !! \details For now, we have routines available only for regular grid and
    !! rectangular physical domain. We also assume that the set of observation
    !! points is a subset of the grid point so that each observation can be
    !! identified by the indices of the associate grid point. These indices are
    !! designated by integer coordinate of the observation.
    !! However, the data structure for observation can be used for any shape of
    !! the physical space domain, any type of discretization grid, providing the
    !! indexes of observation location on the grid. It can also be easily used
    !! with observation points different from the actual grid points.
    !<
    INTERFACE conv_obsgap
        MODULE PROCEDURE conv_obsgap_elt_1D !< obs gap for 1D array of the state vector in 1D physical space
        MODULE PROCEDURE conv_obsgap_elt_2D !< obs gap for 2D array of the state vector in 2D physical space
        !MODULE PROCEDURE conv_obsgap_3D
    END INTERFACE

    INTERFACE ts_cost
        MODULE PROCEDURE ts_cost_1D
        MODULE PROCEDURE ts_cost_1D_single_obs
        MODULE PROCEDURE ts_cost_2D
    END INTERFACE

    INTERFACE ts_costb
        MODULE PROCEDURE ts_costb_1D
        !MODULE PROCEDURE ts_costb_1D_single_obs
        MODULE PROCEDURE ts_costb_2D
    END INTERFACE
CONTAINS

    !> \brief get the size of the control vector
    !! \param[in] td_obs observation data structure
    !<
    FUNCTION get_obsSize(td_obs) RESULT(il_nobs)
        TYPE(obs_structure), INTENT(IN) :: td_obs
        INTEGER :: il_nobs

        il_nobs = td_obs%i_nobs
    END FUNCTION get_obsSize

    !> \brief set observation data with integer coordinates
    !! \param[in,out] td_obs observation data structure
    !! \param[in] ada_fName observation file name
    !<
    SUBROUTINE set_obs_fName(td_obs, ada_fName)
        TYPE(obs_structure), INTENT(INOUT) :: td_obs
        CHARACTER(len=*), INTENT(IN) :: ada_fName
        !local variables

        td_obs%obs_fName = ada_fName
    END SUBROUTINE set_obs_fName

    !> \brief set observation data assuming that the coordinates are already set
    !! \param[in,out] td_obs observation data structure
    !! \param[in] rda_data observation data
    !<
    SUBROUTINE set_obs_data(td_obs, rda_data)
        TYPE(obs_structure), INTENT(INOUT) :: td_obs
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_data

        td_obs%ra_obs = rda_data
    END SUBROUTINE set_obs_data

    !> \brief set observation data with integer coordinates
    !! \param[in,out] td_obs observation data structure
    !! \param[in] rda_data observation data
    !! \param[in] ida_coord coordinates
    !<
    SUBROUTINE set_obs_with_icoord(td_obs, rda_data, ida_coord)
        TYPE(obs_structure), INTENT(INOUT) :: td_obs
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_data
        INTEGER , DIMENSION(:, :), INTENT(IN) :: ida_coord
        !local variables
        td_obs%ra_obs = rda_data
        td_obs%ia_icoord = ida_coord
    END SUBROUTINE set_obs_with_icoord

    !> \brief set observation data with fractional indexes
    !! \param[in,out] td_obs observation data structure
    !! \param[in] rda_data observation data
    !! \param[in] fidx fractional indices
    !<
    subroutine set_obs_with_fidx(td_obs, rda_data, fidx)
        type(obs_structure), intent(inout) :: td_obs
        real(dp), dimension(:), intent(in) :: rda_data
        real(dp), dimension(:, :), intent(in) :: fidx
        !local variables

        td_obs%ra_obs = rda_data
        td_obs%ra_fidx = fidx
    end subroutine set_obs_with_fidx

    !> \brief default initialization of obs_structure
    !! \param[in,out] td_os observation data structure
    !<
    subroutine set_default_obs(td_os)
        type(obs_structure), intent(inout) :: td_os

        td_os%i_ts   = -1
        td_os%r_date = -1.0_dp
        if( (td_os%i_nobs>0).and.(td_os%i_ndim>0) )then
            td_os%ra_obs = 0.0_dp
            td_os%ra_std = 1.0_dp
            td_os%ra_rm1 = 1.0_dp
            td_os%ra_inn = 0.0_dp

            if(td_os%l_icoord) td_os%ia_icoord  = -1
            if(td_os%l_rcoord ) td_os%ra_rcoord = -999.0_dp
            td_os%ra_fidx = -999.0_dp
        end if
    end subroutine set_default_obs

    !> @brief Copy observation data structure
    !! @param[in] td_src data structure to copy from (source)
    !! @param[in] td_dest data structure to copy to (destination)
    !<
    subroutine copy_obs( td_src, td_dest )
        type(obs_structure), intent(in) :: td_src
        type(obs_structure), intent(inout) :: td_dest

        call set_obsSize( td_dest, td_src%i_nobs, td_src%i_ndim, td_src%l_icoord&
        , td_src%l_rcoord, td_src%varName )

        td_dest%i_ts   = td_src%i_ts
        td_dest%r_date = td_src%r_date
        if( (td_dest%i_nobs>0).and.(td_dest%i_ndim>0) )then
            td_dest%ra_obs = td_src%ra_obs
            td_dest%ra_std = td_src%ra_std
            td_dest%ra_Rm1 = td_src%ra_Rm1
            td_dest%ra_inn = td_src%ra_inn

            if(td_dest%l_icoord) td_dest%ia_icoord  = td_src%ia_icoord
            if(td_dest%l_rcoord ) td_dest%ra_rcoord = td_src%ra_rcoord
            td_dest%ra_fidx = td_src%ra_fidx
        end if
    end subroutine copy_obs

  !> \brief Allocates space for observation
  !! \param[in] td_os observation data structure
  !! \param[in] id_nobs total number of observations
  !! \param[in] ndim number of dimension (for coordinates)
  !! \param[in] icoord logical flag saying if integer coordinates are presents
  !! \param[in] rcoord logical flag saying if real coordinates are presents
  !! \param[in] vName variable name
  !! @details For now the space is automatically allocated for fractional coord
  !! @todo add a parameter for fractional coordinates
  !!   Proposition: use a single parameter for all the coordinates flags.
  !!   for example an enumerated variable that take all the possible combinations
  !!   - ALL
  !!   - ICOORD
  !!   - RCOORD
  !!   - FCOORD
  !!   - IRCOORD
  !!   - IFCOORD
  !!   - FRCOORD
  !!
  !!   and uses a structure like:
  !!   if (coordFlag in [ALL, ICOORD, IRCOORD, IFCOORD]) allocate icoord
  !!   if (coordFlag in [ALL, FCOORD, FRCOORD, IFCOORD]) allocate fcoord
  !!   if (coordFlag in [ALL, RCOORD, IRCOORD, FRCOORD]) allocate rcoord
  !<
  SUBROUTINE set_obsSize( td_os, id_nobs, ndim, icoord, rcoord, vName )
    TYPE(obs_structure), INTENT(INOUT) :: td_os
    INTEGER, INTENT(IN) :: id_nobs, ndim
    LOGICAL, INTENT(IN) :: icoord, rcoord
    character(len=*), optional, intent(in) :: vName

    !if the current size or dimensionality are different from the required ones, reallocate
    !call debug([td_os%i_nobs, id_nobs] , "set_obsSize, [td_os%i_nobs, id_nobs]  = ")
    !call debug([td_os%i_ndim, ndim]    , "set_obsSize, [td_os%i_ndim, ndim]     = ")
    !call debug(int([td_os%l_icoord, icoord]), "set_obsSize, [td_os%l_icoord, icoord] = ")
    !call debug(int([td_os%l_rcoord, rcoord]), "set_obsSize, [td_os%l_rcoord, rcoord] = ")
    IF( (td_os%i_nobs /= id_nobs).OR.(td_os%i_ndim /= ndim).OR.&
        (td_os%l_icoord .NEQV. icoord).OR.(td_os%l_rcoord .NEQV. rcoord) )THEN
      !if the current size and dimensionality are nonzero, then deallocate
      IF( (td_os%i_nobs>0).AND.(td_os%i_ndim>0) )THEN
        DEALLOCATE(td_os%ra_obs, td_os%ra_std, td_os%ra_Rm1, td_os%ra_inn,&
                   td_os%ra_fidx)

        IF(td_os%l_rcoord) DEALLOCATE(td_os%ra_rcoord)
        IF(td_os%l_icoord ) DEALLOCATE(td_os%ia_icoord)
        nullify(td_os%ra_obs, td_os%ra_std, td_os%ra_Rm1, td_os%ra_inn&
               , td_os%ra_obsgap, td_os%ra_sigma &
               , td_os%ra_rcoord, td_os%ia_icoord, td_os%ra_fidx)
      END IF
      ! if the required size and dimensionality are nonzero, then allocate
      IF( (id_nobs>0).AND.(ndim>0) )THEN
        ALLOCATE( td_os%ra_obs(id_nobs), td_os%ra_inn(id_nobs),&
                  td_os%ra_std(id_nobs), td_os%ra_Rm1(id_nobs),&
                  td_os%ra_fidx(ndim, id_nobs) &
        )
        td_os%ra_obsgap => td_os%ra_inn
        td_os%ra_sigma  => td_os%ra_std

        IF(icoord) ALLOCATE( td_os%ia_icoord(ndim, id_nobs) )
        IF(rcoord) ALLOCATE( td_os%ra_rcoord(ndim, id_nobs) )
        CALL set_default_obs(td_os)
      END IF
        td_os%i_nobs   = id_nobs
        td_os%i_ndim   = ndim
        td_os%l_rcoord = rcoord
        td_os%l_icoord = icoord
    END IF
    td_os%varName = vName
  END SUBROUTINE set_obsSize

    !> @brief reset the obs data structure, free memory space
    !! @param[in,out] obs
    !! @param[in,out] o1, o2, o3, o4, o5, o6, o7, o8, o9, o10 additional observation to reset
    !! @details
    !!   this procedure can reset a maximum of 11 obs data structure
    !<
    subroutine reset_obs(obs, o1, o2, o3, o4, o5, o6, o7, o8, o9, o10)
        type(obs_structure), intent(in out) :: obs
        type(obs_structure), optional, intent(in out) :: &
            o1, o2, o3, o4, o5, o6, o7, o8, o9, o10

        call set_obsSize(obs, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o1))call set_obsSize(o1, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o2))call set_obsSize(o2, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o3))call set_obsSize(o3, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o4))call set_obsSize(o4, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o5))call set_obsSize(o5, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o6))call set_obsSize(o6, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o7))call set_obsSize(o7, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o8))call set_obsSize(o8, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o9))call set_obsSize(o9, 0, 0, icoord=.false. , rcoord=.false.)
        if(present(o10))call set_obsSize(o10, 0, 0, icoord=.false. , rcoord=.false.)
    end subroutine reset_obs

    !> @brief remove rejected observations from the data structure
    !! @param[in,out] obs
    !! @param[in] rej rejection status, .true. if rejected
    !! @details
    !<
    subroutine remove_rejected(obs, rej)
        type(obs_structure), intent(in out) :: obs
        logical, dimension(:), intent(in) :: rej
        !local variables
        integer :: nKept, k, pos
        type(obs_structure) :: tl_obs

        !call debug('Entering remove_rejected', tag=dALLWAYS)
        if(size(rej)/=obs%i_nobs)then
            call stop_program('In conv_obs::remove_rejected, incompatible obs and rej size')
        end if

        nKept = count( .not.rej )
        if((obs%i_nobs>nKept).and.(nKept>0))then
            call set_obsSize(tl_obs, nKept, obs%i_ndim, icoord=obs%l_icoord&
                , rcoord=obs%l_rcoord, vName=obs%varName)

            pos=1
            do k = 1, obs%i_nobs
                if(.not.rej(k))then
                    tl_obs%ra_obs(pos) = obs%ra_obs(k)
                    tl_obs%ra_inn(pos) = obs%ra_inn(k)
                    tl_obs%ra_std(pos) = obs%ra_std(k)
                    tl_obs%ra_fidx(:, pos) = obs%ra_fidx(:, k)
                    if(obs%l_icoord)then
                        tl_obs%ia_icoord(:, pos) = obs%ia_icoord(:, k)
                    end if
                    if(obs%l_rcoord)then
                        tl_obs%ra_rcoord(:, pos) = obs%ra_rcoord(:, k)
                    end if
                    pos = pos+1
                end if
            end do

            call set_obsSize(obs, nKept, tl_obs%i_ndim, icoord=tl_obs%l_icoord&
            , rcoord=tl_obs%l_rcoord, vName=tl_obs%varName)

            obs%ra_obs = tl_obs%ra_obs
            obs%ra_inn = tl_obs%ra_inn
            obs%ra_std = tl_obs%ra_std
            obs%ra_fidx =tl_obs%ra_fidx
            if(obs%l_icoord) obs%ia_icoord = tl_obs%ia_icoord
            if(obs%l_rcoord) obs%ra_rcoord = tl_obs%ra_rcoord

        else if(nKept==0)then
            call reset_obs(obs)
        end if
        !call debug('Exiting remove_rejected', tag=dALLWAYS)
    end subroutine remove_rejected

  subroutine print_os_info(td_os)
    type(obs_structure), intent(in) :: td_os

    call print_var(td_os%i_nobs      , '  total number of observations    = ')
    call print_var(td_os%i_obs_level , '  Observation level for wavelet   = ')
    if( associated(td_os%ra_rcoord) )&
         call print_var(td_os%ra_rcoord, '  real coordinates               = ')
    if( associated(td_os%ia_icoord) )&
         call print_var(td_os%ia_icoord, '  int coordinates            = ')
  end subroutine print_os_info

    !>\brief print observation for diagnostics
    SUBROUTINE print_os(td_os)
        TYPE(obs_structure), INTENT(IN) :: td_os

        WRITE(*,                 *) 'com_tools::print_obs : obs_structure-----------------------'
        CALL print_var(td_os%obs_fName, '  Observation file name.......... = ')
        CALL print_os_info(td_os)
        CALL print_var(td_os%ra_obs      , '  observation data............... = ')
        CALL print_var(td_os%ra_std    , '  standard deviation in obs error = ')
    END SUBROUTINE print_os

    !>\brief print obsgap for diagnostics
    SUBROUTINE print_os_gap(td_os)
        TYPE(obs_structure), INTENT(IN) :: td_os

        WRITE(*,                 *) 'com_tools::print_obs : obs_structure-----------------------'
        CALL print_var(td_os%ogap_fName, '  Obs gap file name...... = ')
        CALL print_os_info(td_os)
        CALL print_var(td_os%ra_inn   , '  innovation data .......... = ')
        CALL print_var(td_os%ra_Rm1      , '  diagonal of the cov matrix = ')
    END SUBROUTINE print_os_gap

    !> \brief obsgap routine for a 1D single variable
    !! \param [in,out] tda_obs array of obs_structure
    !! \param [in] rda_state state variable
    !! \param [in] id_ts time step associated with the state variable
    !! \details This routine is designed for a problem in 1 physical dimension
    !! with a state variable represented by a 1D array at each time step. For now,
    !! it is assumed that observation are taken at computation time and
    !! computation grid points so that there is no need for interpolation. tda_obs
    !! is a vector of obs_structure, its size is given by the number of time step
    !! where observation are sampled. Each element of this vector contains the
    !! observation for a given observation time
    !<
    SUBROUTINE conv_obsgap_elt_1D(tda_obs, rda_state, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT)    :: tda_obs
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_state
        INTEGER, INTENT(IN) :: id_ts
        INTEGER :: il_idx, ibi, il_ix

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix = tda_obs(il_idx)%ia_icoord(1, ibi)
                tda_obs(il_idx)%ra_inn(ibi) = rda_state( il_ix ) - tda_obs(il_idx)%ra_obs(ibi)
            END DO
        END IF
    END SUBROUTINE conv_obsgap_elt_1D

    !> \brief compute the cost for a given time step
    !! \param[in,out] tda_obs array of observation data structures
    !! \param[in] rda_state state variable
    !! \param[in] id_ts time step associated with the adjoint state variable
    !! \details This routine is designed for a problem in 1 physical dimension with
    !! a state variable represented by a 1D array at each time step. For now, it is
    !! assumed that observation are taken at computation time and computation grid
    !! points so that there is no need for interpolation. tda_obs is a vector of
    !! obs_structure, its size is given by the number of time step where
    !! observation are sampled. Each element of this vector contains the
    !! observation for a given observation time. This routine saves the obsgap so
    !! that it is not necessary to keep
    !<
    SUBROUTINE ts_cost_1D(tda_obs, rda_state, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT) :: tda_obs
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_state
        INTEGER, INTENT(IN) :: id_ts
        INTEGER :: il_idx, ibi, il_ix

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix = tda_obs(il_idx)%ia_icoord(1, ibi)
                !the next line of code helps to not passed the state to the associated adjoint routine
                tda_obs(il_idx)%ra_inn(ibi) = rda_state( il_ix ) - tda_obs(il_idx)%ra_obs(ibi)
            END DO
            tda_obs(il_idx)%r_cost = 0.5_cp*SUM( tda_obs(il_idx)%ra_Rm1 * tda_obs(il_idx)%ra_inn**2 )
        END IF
    END SUBROUTINE ts_cost_1D

    !> \brief compute the cost for a given time step
    !! \param[in,out] td_obs array of observation data structures
    !! \param[in] rda_state state variable
    !! \details This routine is designed for a problem in 1 physical dimension
    !! with a state variable represented by a 1D array at each time step. For now,
    !! it is assumed that observation are taken at computation time and
    !! computation grid points so that there is no need for interpolation. tda_obs
    !! is a vector of obs_structure, its size is given by the number of time step
    !! where observation are sampled. Each element of this vector contains the
    !! observation for a given observation time. This routine saves the obsgap so
    !! that it is not necessary to keep
    !<
    SUBROUTINE ts_cost_1D_single_obs(td_obs, rda_state)
        TYPE(obs_structure), INTENT(INOUT) :: td_obs
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_state
        INTEGER :: ibi, il_ix

        DO ibi = 1, SIZE(td_obs%ra_obs)
        il_ix = td_obs%ia_icoord(1, ibi)
        !the next line of code helps to not pass the state to the associated adjoint routine
        td_obs%ra_inn(ibi) = rda_state( il_ix ) - td_obs%ra_obs(ibi)
        END DO
        td_obs%r_cost = 0.5_cp*SUM( td_obs%ra_Rm1 * td_obs%ra_inn**2 )
    END SUBROUTINE ts_cost_1D_single_obs

    !> \brief compute the cost for a given time step, 2 physical dimensions
    !! \param[in,out] tda_obs array of observation data structures
    !! \param[in] rda_state state variable
    !! \param[in] id_ts time step associated with the adjoint state variable
    !! \details This routine is designed for a problem in 2 physical dimension
    !! with a state variable represented by a 2D array at each time step. For now,
    !! it is assumed that observation are taken at computation time and
    !! computation grid points so that there is no need for interpolation. tda_obs
    !! is a vector of obs_structure, its size is given by the number of time step
    !! where observation are sampled. Each element of this vector contains the
    !! observation for a given observation time. This routine saves the obsgap so
    !! that it is not necessary to keep
    !<
    SUBROUTINE ts_cost_2D(tda_obs, rda_state, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT) :: tda_obs
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: rda_state
        INTEGER, INTENT(IN) :: id_ts
        INTEGER :: il_idx, ibi, il_ix1, il_ix2

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix1 = tda_obs(il_idx)%ia_icoord(1, ibi)
                il_ix2 = tda_obs(il_idx)%ia_icoord(2, ibi)
                !the next line of code helps to not passed the state to the associated adjoint routine
                tda_obs(il_idx)%ra_inn(ibi) = rda_state( il_ix1, il_ix2 ) - tda_obs(il_idx)%ra_obs(ibi)
            END DO
            tda_obs(il_idx)%r_cost = 0.5_cp*SUM( tda_obs(il_idx)%ra_Rm1 * tda_obs(il_idx)%ra_inn**2 )
            !CALL debug(tda_obs(il_idx)%ra_Rm1, 'In ts_cost_2D, tda_obs(il_idx)%ra_Rm1 = ',tag=dALLWAYS)
        END IF
    END SUBROUTINE ts_cost_2D

    !> \brief obsgap routine for a 2D single variable
    !! \param[in,out] tda_obs array of obs_structure
    !! \param[in] rda_state state variable
    !! \param[in] id_ts time step associated with the state variable
    !! \details This routine is design for a problem in 2 physical dimensions with
    !! a state variable represented by a 2D array at each time step. For now, it
    !! is assumed that observation are taken at computation time and computation
    !! grid points so that there is no need for interpolation. tda_obs is a vector
    !! of obs_structure, its size is given by the number of time step where
    !! observation are sampled. Each element of this vector contains the observation
    !! for a given observation time
    !<
    SUBROUTINE conv_obsgap_elt_2D(tda_obs, rda_state, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT)    :: tda_obs
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_state
        INTEGER, INTENT(IN) :: id_ts
        INTEGER :: il_idx, ibi, il_ix, il_iy

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix = tda_obs(il_idx)%ia_icoord(1, ibi)
                il_iy = tda_obs(il_idx)%ia_icoord(2, ibi)
                tda_obs(il_idx)%ra_inn(ibi) = rda_state( il_ix, il_iy ) - tda_obs(il_idx)%ra_obs(ibi)
            END DO
        END IF
    END SUBROUTINE conv_obsgap_elt_2D

    !> Checks if observation are supposed to be sampled at time step \a id_ts
    !! @param [in] tda_obs vector of observation (data structure)
    !! @param [in] id_ts time step of interest
    !! @param [in] id_idx index in \a tda_obs of the observation of interest if it exists
    !!
    !!
    !<
    FUNCTION find_conv_obs( tda_obs, id_ts, id_idx )RESULT(ll_found)
        TYPE(obs_structure), DIMENSION(:), intent(in) :: tda_obs
        INTEGER, INTENT(IN)  :: id_ts!index of the time step of interest
        INTEGER, INTENT(OUT) :: id_idx! Index of the particular obs relative to the numbers of obs
        LOGICAL              :: ll_found
        INTEGER              :: ibi

        id_idx  = -1 !
        ll_found = .FALSE.
        DO ibi = SIZE(tda_obs),1,-1
            IF( id_ts == tda_obs(ibi)%i_ts ) THEN
                id_idx = ibi
                ll_found = .TRUE.
            END IF
        END DO

    END FUNCTION find_conv_obs

    !> \brief set the number of dimension of the problem
    !! \param[in,out] td_cop data structure containing the conventional observation parameters
    !! \param[in] id_ndim number of dimension of the problem
    !! \details this routine also allocate space for dynamic fields of the data structure
    !<
    SUBROUTINE set_cop_ndim(td_cop, id_ndim)
        TYPE(convObs_param), INTENT(INOUT) :: td_cop
        INTEGER, INTENT(IN) :: id_ndim

        IF(td_cop%i_ndim>0)DEALLOCATE( td_cop%ia_shifts, td_cop%ia_steps, td_cop%ia_max )

        IF( id_ndim>0 )THEN
        td_cop%i_ndim = id_ndim
        ALLOCATE(&
            td_cop%ia_shifts( id_ndim ),&
            td_cop%ia_steps ( id_ndim ),&
            td_cop%ia_max   ( id_ndim ) &
        )
        END IF
    END SUBROUTINE set_cop_ndim

    !> \brief computes the sampling time steps for twin observation
    !! \param[in,out] td_cop data structure containing the conventional observation parameters
    !<
    SUBROUTINE compute_cop_samplingTS(td_cop)
        TYPE(convObs_param), INTENT(INOUT) :: td_cop
        !local variables
        INTEGER :: il_tmp
        INTEGER, DIMENSION(1) :: ila_ncoordTS
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ila_tmpTS

        il_tmp = regular_coord_count(&
            ida_max=(/td_cop%i_tmax/), shifts=(/td_cop%i_tshift/),&
            steps=(/td_cop%i_tstep/), ncoords=ila_ncoordTS&
        )
        IF( ASSOCIATED(td_cop%ia_samplingTS) )DEALLOCATE( td_cop%ia_samplingTS )
        ALLOCATE(&
        ila_tmpTS( 1, ila_ncoordTS(1) ),&
        td_cop%ia_samplingTS( ila_ncoordTS(1) )&
        )
        CALL build_regular_int_coord(&
            ila_tmpTS, ncoords=ila_ncoordTS, &
            shifts=(/td_cop%i_tshift/), steps=(/td_cop%i_tstep/) &
        )
        td_cop%ia_samplingTS(:) = ila_tmpTS(1, :)
    !     CALL debug(ila_ncoordTS, 'In simulator_init_assim, ila_ncoordTS = ', tag=dALLWAYS)
    !     CALL debug(ila_tmpTS, 'In simulator_init_assim, ila_tmpTS = ', tag=dALLWAYS)
    !     CALL debug(td_cop%ia_samplingTS, 'In simulator_init_assim, td_cop%ia_samplingTS = ', tag=dALLWAYS)
    !     CALL debug(SHAPE(td_cop%ia_samplingTS), 'In simulator_init_assim, SHAPE(td_cop%ia_samplingTS) = ', tag=dALLWAYS)
    END SUBROUTINE compute_cop_samplingTS

    !> \brief print conventional observation parameters
    !! \param[in] td_cop data structure containing the conventional observation parameters
    !!
    !<
    SUBROUTINE print_cop(td_cop)
        TYPE(convObs_param), INTENT(IN) :: td_cop

        CALL debug('', 'printing convObs_param---------------------------------', tag=dALLWAYS)
        CALL debug(td_cop%aa_title   , '  aa_title...', tag=dALLWAYS)
        CALL debug(td_cop%aa_hist    , '  aa_hist....', tag=dALLWAYS)
        CALL debug(td_cop%ia_shifts  , '  ia_shifts..', tag=dALLWAYS)
        CALL debug(td_cop%ia_steps   , '  ia_steps...', tag=dALLWAYS)
        CALL debug(td_cop%ia_max     , '  ia_max.....', tag=dALLWAYS)
        CALL debug(td_cop%i_tshift   , '  i_tshift...', tag=dALLWAYS)
        CALL debug(td_cop%i_tstep    , '  i_tstep....', tag=dALLWAYS)
        CALL debug(td_cop%i_tmax     , '  i_tmax.....', tag=dALLWAYS)
        CALL debug(td_cop%r_sigmaR   , '  r_sigmaR...', tag=dALLWAYS)
        CALL debug(td_cop%l_sampleObs, '  l_sampleObs', tag=dALLWAYS)
        CALL debug(td_cop%aa_prefix  , '  aa_prefix..', tag=dALLWAYS)
        CALL debug(td_cop%r_dt       , '  r_dt.......', tag=dALLWAYS)
        CALL debug(td_cop%ia_samplingTS, '  ia_samplingTS= ', tag=dALLWAYS)
        CALL debug('', '.......................................................', tag=dALLWAYS)
    END SUBROUTINE print_cop

!   !> Copies the data structure src to dest
!   !! \param[in] src data structure to be copied
!   !! \param[in] dest destination data structure
!   !! \param[in] destVarName name of the destination variable
!   !!
!   !<
!   subroutine copy_cop(src, dest, destVarName)
!     TYPE(convObs_param), INTENT(IN) :: td_cop
!     TYPE(convObs_param), INTENT(INOUT) :: td_cop
!     CHARACTER(len=*), INTENT(IN) :: ada_namelist, ada_varName
!     !
!
!   end subroutine copy_cop

  !> \brief Loads conventional observation parameters from a namelist file
  !! \param[in,out] td_cop data structure to be loaded
  !! \param[in] ada_namelist name of the namelist file
  !! \param[in] ada_varName name of the variable associated with td_cop
  !! \param[in] ndim number of dimensions of the physical space
  !! \param[in] tshift shift in the time dimension
  !! \param[in] nt (optional) maximum number of time steps at which the last observation can be sampled.
  !! \param[in] shifts (optional) shifts along the physical space, ia_shifts(i) gives the min index of the observation grid point in the ith dimension
  !! \param[in] maxi (optional) !> max along the physical space, ia_max(i) gives the max index of the observation grid point in the ith dimension
  !! \param[in] dt (optional) time steps used for model integration
  !! \param[in] compute (optional) says if the unloaded parameters (samplingTS) should be computed, by default, those parameters are computed
  !! \param[in] newVarName (optional) var name to be used in the data structure. This parameter is useful when the observation are sampled at the same location for many variable. In that case, the input file describes only one that is loaded for all
  !!
  !! \details This subroutine will cause the program to stop if there is no NAM_OBS bloc in ada_namelist with varName equal to ada_varName.
  !! if present, \a nt does not give the number of observation, it gives the
  !! number of time steps at wich observation are authorized. The number of
  !! observation takes into account the number of time steps between successive
  !! observations. Some values from the input file are ovewritten if the
  !! corresponding argument if present: tshift -> il_tshift, nt -> il_tmax,
  !! shifts -> ila_shifts, max -> ila_max. This can be useful in a test case
  !! when discretization parameters changed frequently, and that one can forget
  !! to update observation parameters in the input file.
  !<
  SUBROUTINE load_cop(td_cop, ada_namelist, ada_varName, ndim, tshift, nt, shifts, maxi, compute, dt, newVarName)
    INTEGER, PARAMETER:: ip_numnam = 68
    TYPE(convObs_param), INTENT(INOUT) :: td_cop
    CHARACTER(len=*), INTENT(IN) :: ada_namelist, ada_varName
    INTEGER, INTENT(IN) :: ndim
    INTEGER, OPTIONAL, INTENT(IN) :: nt, tshift
    REAL(cp), OPTIONAL, INTENT(IN) :: dt
    INTEGER, DIMENSION(ndim), OPTIONAL, INTENT(IN) :: shifts, maxi
    LOGICAL, OPTIONAL, INTENT(IN) :: compute
    CHARACTER(len=*), INTENT(IN),OPTIONAL :: newVarName
    !local var
    CHARACTER(len=ip_snl) :: varName, title, history, ala_varName
    REAL(cp) :: rl_sigmaR
    LOGICAL :: sampleObs, ll_compute
    INTEGER :: il_tshift, il_tstep, il_tmax
    INTEGER, DIMENSION(3) :: ila_shift, ila_steps, ila_max !a physical space is supposed to have at most 3 dimensions
    NAMELIST/NAM_OBS/&
      varName  ,&
      sampleObs,&
      title    ,&
      history  ,&
      rl_sigmaR,&
      ila_shift,&
      ila_steps,&
      ila_max  ,&
      il_tshift,&
      il_tstep ,&
      il_tmax

    OPEN(ip_numnam,FILE=ada_namelist,FORM='FORMATTED',STATUS='OLD')
    varName = "{#}@"
    ala_varName = TRIM(ada_varName)
    CALL uppercase(ala_varName)
    DO WHILE(varName.NE.ala_varName)
      READ(ip_numnam, NAM_OBS)!reading the block
      CALL uppercase(varName)
      CALL debug(varName, "In load_cop, found varName", tag=dALLWAYS)
    END DO
    CLOSE(ip_numnam)

    IF( PRESENT(nt) )THEN
      il_tmax = nt
    END IF
    IF( PRESENT(tshift ) )THEN
      il_tshift = tshift
    END IF
    IF( PRESENT(shifts) )THEN
      ila_shift(1:ndim) = shifts
    END IF
    IF( PRESENT(maxi) )THEN
      ila_max(1:ndim) = maxi
    END IF
    IF( PRESENT(compute) )THEN
      ll_compute = compute
    ELSE
      ll_compute = .TRUE.
    END IF

    CALL set_cop_ndim( td_cop, ndim )
    td_cop%aa_title   = title
    td_cop%aa_hist    = history
    IF( ANY(ila_shift(1:ndim)<0) )THEN
      CALL stop_program('In load_cop, the shift(s) in space should not be negative for variable '//ada_varName)
    ELSE
      td_cop%ia_shifts  = ila_shift(1:ndim)
    END IF
    IF( ANY(ila_steps(1:ndim)<0) )THEN
      CALL stop_program('In load_cop, the step(s) in space should not be negative for variable '//ada_varName)
    ELSE
      td_cop%ia_steps   = ila_steps(1:ndim)
    END IF
    IF( ANY(ila_max(1:ndim)<0) )THEN
      CALL stop_program('In load_cop, the max(s) in space should not be negative for variable '//ada_varName)
    ELSE
      td_cop%ia_max     = ila_max(1:ndim)
    END IF
    IF(il_tshift<0)THEN
      CALL stop_program('In load_cop, the shift in time should not be negative for variable '//ada_varName)
    ELSE
      td_cop%i_tshift   = il_tshift
    END IF
    IF(il_tstep<0)THEN
      CALL stop_program('In load_cop, the step in time should not be negative for variable '//ada_varName)
    ELSE
    td_cop%i_tstep    = il_tstep
    END IF
    IF(il_tmax<0)THEN
      CALL stop_program('In load_cop, the max in time should not be negative for variable '//ada_varName)
    ELSE
    td_cop%i_tmax     = il_tmax
    END IF
    td_cop%r_sigmaR   = rl_sigmaR
    td_cop%l_sampleObs= sampleObs
    CALL lowercase(varName)
    if(present(newVarName))then
      td_cop%aa_prefix  = TRIM( toLower(newVarName) )
      td_cop%aa_title   = trim(title)//': '// trim(td_cop%aa_prefix)
    else
      td_cop%aa_prefix  = TRIM( varName )
    end if
    IF( PRESENT( dt ) )THEN
      td_cop%r_dt = dt
    ELSE
      td_cop%r_dt = -1.0_cp
    END IF
    IF(ll_compute) CALL compute_cop_samplingTS(td_cop)
  END SUBROUTINE load_cop

!   !> \brief initializes data structures for observation sampling in twin experiments
!   !! \param[in,out] td_obs data structure for observations
!   !! \param[in,out] td_cop data structure containing observation parameters
!   !!
!   !<
!   SUBROUTINE init_conv_obs_sampling(td_obs, td_cop)
!     TYPE(obs_structure), INTENT(INOUT) :: td_obs
!     TYPE(convObs_param), INTENT(IN) :: td_cop
!
!     IF(td_cop%l_make_obs)THEN
!       CALL initObsOutput(&
!           tm_uobs_out, td_cop%aa_fName, id_nobs=il_nobs, &
!           id_ndim=il_ndim, dt=td_cop%r_dt &
!       )
!
!       tm_uobs_out%aa_title   = "Observations of the u component of the velocity"
!       tm_uobs_out%aa_history = make_ncHistory("BALAISE direct model for twin experiments")
!       CALL ncObsCreate(tm_uobs_out)
!       CALL set_obsSize(tm_uobs, il_nobs, id_coord_ndim=il_ndim, icoord=.TRUE., rcoord=.FALSE.)
!       CALL set_default_obs(tm_uobs)
!       CALL build_regular_int_coord(&
!           tm_uobs%ia_icoord,  ncoords=ila_ncoords(1:il_ndim), &
!           shifts=td_cop%ia_shifts(1:il_ndim),   steps=td_cop%ia_steps(1:il_ndim) &
!       )
!       !CALL debug(tm_uobs%ia_icoord, '  tm_uobs%ia_icoord = ')
!       tm_uobs%ra_std = td_cop%r_u_sigma
!     END IF
! END SUBROUTINE init_conv_obs_sampling

!---------------------------- Adjoint code section ------------------------------
    !> Adjoint routine associated with ts_cost_1D
    !! \details this routine is an optimuized adjoint and do not have all the details of a mechanical adjoint code
    !! @param tda_obs, rda_stateb, id_ts for details, @see ts_cost_2D
    !<
    SUBROUTINE ts_costb_1D(tda_obs, rda_stateb, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT) :: tda_obs
        REAL(dp), DIMENSION(:), INTENT(IN OUT) :: rda_stateb
        INTEGER, INTENT(IN) :: id_ts
        !local variables
        REAL(dp), DIMENSION(:), ALLOCATABLE :: rla_ogapb
        INTEGER :: il_idx, ibi, il_ix

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            ALLOCATE( rla_ogapb( SIZE(tda_obs(il_idx)%ra_obs) ) )
            rla_ogapb = 0.0_cp
            rla_ogapb = rla_ogapb + tda_obs(il_idx)%r_costb * tda_obs(il_idx)%ra_Rm1 * tda_obs(il_idx)%ra_inn
            tda_obs(il_idx)%r_costb = 0.0_cp
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix = tda_obs(il_idx)%ia_icoord(1, ibi)
                rda_stateb( il_ix ) = rda_stateb( il_ix ) + rla_ogapb(ibi)
                rla_ogapb(ibi) = 0.0_cp
            END DO
        END IF
    END SUBROUTINE ts_costb_1D

    !> Adjoint routine associated with ts_cost_2D
    !! @details this routine is an optimuized adjoint and do not have all the details of a mechanical adjoint code
    !! @param tda_obs, rda_stateb, id_ts for details, @see ts_cost_2D
    !<
    SUBROUTINE ts_costb_2D(tda_obs, rda_stateb, id_ts)
        TYPE(obs_structure), DIMENSION(:), INTENT(INOUT) :: tda_obs
        REAL(dp), DIMENSION(:,:), INTENT(IN OUT) :: rda_stateb
        INTEGER, INTENT(IN) :: id_ts
        !local variables
        REAL(dp), DIMENSION(:), ALLOCATABLE :: rla_ogapb!there is no field obsgapb in the observation data structure
        INTEGER :: il_idx, ibi, il_ix1, il_ix2

        IF( find_conv_obs( tda_obs, id_ts, il_idx ) )THEN
            ALLOCATE( rla_ogapb( SIZE(tda_obs(il_idx)%ra_obs) ) )
            rla_ogapb = 0.0_cp
            rla_ogapb = rla_ogapb + tda_obs(il_idx)%r_costb * tda_obs(il_idx)%ra_Rm1 * tda_obs(il_idx)%ra_inn
            tda_obs(il_idx)%r_costb = 0.0_cp
            DO ibi = 1, SIZE(tda_obs(il_idx)%ra_obs)
                il_ix1 = tda_obs(il_idx)%ia_icoord(1, ibi)
                il_ix2 = tda_obs(il_idx)%ia_icoord(2, ibi)
                rda_stateb( il_ix1, il_ix2 ) = rda_stateb( il_ix1, il_ix2 ) + rla_ogapb(ibi)
                rla_ogapb(ibi) = 0.0_cp
            END DO
        END IF
    END SUBROUTINE ts_costb_2D
end module conv_obs
