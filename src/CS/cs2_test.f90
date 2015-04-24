PROGRAM cs2_test
  USE debug_tools
  USE com_tools
  USE fftw_tools
  USE cs_tools
  USE general_tools
IMPLICIT NONE
  REAL(cp), DIMENSION(:), ALLOCATABLE :: rga_state , rga_obsgap
  REAL(cp), DIMENSION(:), ALLOCATABLE :: rga_stateb, rga_obsgapb !adjoint variables
  TYPE(exchange_param) :: tg_ep
  TYPE(obs_structure)  :: tg_obs
  REAL(cp) :: rg_costb

  CALL init_all()
  !main program
  !CALL debug('Entering model program*******************')
  CALL read_ep_data(tg_ep, EP_DATA, RTIME_FILE)
  !CALL debug('Entering model program*******************')
  !tg_obs%i_nb_dim  = 1
  tg_obs%obs_fName  = make_fileName(tg_ep, OBS_DATA, INPUT_FILE)
  tg_obs%ogap_fName = make_fileName(tg_ep, OGAP_DATA, RTIME_FILE)

  !CALL print_ep(tg_ep)
  CALL debug(tg_ep%aa_solver_action, 'Running cs2_test for ')
  SELECT CASE(tg_ep%aa_solver_action)
    CASE (RUN_COST) !cost function recquired
      CALL init_cost(tg_ep, tg_obs)
      CALL cost(tg_ep%ra_ctl, tg_ep%r_cost)
      CALL write_ep_data(tg_ep, COST_DATA, RTIME_FILE, 'RUN_COST')
    CASE (RUN_GRADIENT) !gradient of the cost function recquired
      CALL init_grad(tg_ep, tg_obs)
      CALL adjoint_zeroing(tg_ep)
      CALL costb(tg_ep%ra_ctl, tg_ep%ra_grad, tg_ep%r_cost, rg_costb)
      CALL write_ep_data(tg_ep, GRAD_DATA, RTIME_FILE, 'RUN_GRADIENT')
      !CALL print_ep(tg_ep)
    CASE (MAKE_OBS, MAKE_AOBS, MAKE_TOBS) !make observations for twin experiments
      !CALL print_ep(tg_ep)
      CALL make_twin_obs(tg_ep)
    CASE (MAKE_OBS_FROM_COORD) !make observations from given coordinates for twin experiments
      CALL make_twin_obs_from_coord(tg_ep)
    CASE (MAKE_CTL)! generates a control vector, the true trajectory and the true obs
      CALL make_random_ctl(tg_ep)
    CASE (MAKE_BG)! generates and saves a background control vector
      CALL make_zero_bg(tg_ep)
    CASE (RUN_DIRECT, MAKE_DMT, MAKE_ADMT, MAKE_TDMT) ! run direct model and save trajectory
      CALL make_model_trajectory(tg_ep)
      !WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_solver_action)//'>'
    CASE (RUN_ADJOINT)! adjoint model run recquired
      !CALL init_adjoint()
      !CALL adjoint(...)
      WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_solver_action)//'>'
    CASE DEFAULT
      CALL stop_program('Unknow action <'//TRIM(tg_ep%aa_solver_action)//'>')
  END SELECT
  CALL debug('Ending model program  *******************')
  !end of main program
  CALL finalize_all()
CONTAINS

  SUBROUTINE init_all()
    CALL init_fftw()
    CALL init_cs_tools()
  END SUBROUTINE init_all

  SUBROUTINE finalize_all()
    CALL finalize_fftw()
    CALL finalize_cs_tools()
  END SUBROUTINE finalize_all

  !> initializes direct model
  SUBROUTINE init_direct(td_ep)
    TYPE(exchange_param), INTENT(IN) :: td_ep

    ALLOCATE( rga_state( get_ctlSize(td_ep) ) )
  END SUBROUTINE init_direct

  !> initializes the environment for the computation of the cost function
  SUBROUTINE init_cost(td_ep, td_obs)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    TYPE(obs_structure), INTENT(INOUT) :: td_obs

    CALL init_direct(td_ep)
    CALL read_obs(td_obs)
    ALLOCATE( rga_obsgap ( get_obsSize(td_obs) ) )
  END SUBROUTINE init_cost

  !> initializes the environment for the computation of the gradient of the cost function
  !<
  SUBROUTINE init_grad(td_ep, td_obs)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    TYPE(obs_structure), INTENT(INOUT) :: td_obs

    CALL read_obsgap(td_obs)
    !CALL print_os(td_obs)
    ALLOCATE( rga_state  ( get_ctlSize(td_ep ) ),&
              rga_stateb ( get_ctlSize(td_ep ) ),&
              rga_obsgap ( get_obsSize(td_obs) ),&
              rga_obsgapb( get_obsSize(td_obs) ) &
    )
  END SUBROUTINE init_grad

  !> zeroing adjoint variables
  !<
  SUBROUTINE adjoint_zeroing(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    rga_stateb    = 0.0_cp
    rga_obsgapb   = 0.0_cp
    td_ep%ra_grad = 0.0_cp
    rg_costb      = 1.0_cp
  END SUBROUTINE adjoint_zeroing

  !> direct model for the compressive sensing problem
  SUBROUTINE direct(rda_ctl, rda_state)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_state
    !local variable
    REAL(cp), DIMENSION( SIZE(rda_ctl) ) :: rla_ctl
    !CALL debug('','In direct ')
    rla_ctl = rda_ctl + tg_ep%ra_b_ctl
    CALL ifftw_trans(rla_ctl, rda_state)
  END SUBROUTINE direct

  !> adjoint model for the compressive sensing problem
  SUBROUTINE adjoint(rda_ctl, rda_ctlb, rda_state, rda_stateb)
    REAL(cp), DIMENSION(:), INTENT(IN)    :: rda_ctl
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_ctlb
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_state
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_stateb
    !local variable
    REAL(cp), DIMENSION( SIZE(rda_ctl) ) ::rla_ctlb
    rla_ctlb     = 0.0_cp
    !
    CALL nothing(rda_ctl)
    CALL nothing(rda_state)
    !
    !recomputation
    !adjoint code
    !CALL debug(SUM(rda_stateb), 'In adjoint: SUM(rda_stateb) = ')
    !CALL debug(SUM(rla_ctlb_tmp), 'In adjoint: SUM(rla_ctlb_tmp) = ')
    CALL fftw_trans(rda_stateb, rla_ctlb)
    !CALL debug(SUM(rda_stateb), 'In adjoint: SUM(rda_stateb) = ')
    !CALL debug(SUM(rla_ctlb_tmp), 'In adjoint: SUM(rla_ctlb_tmp) = ')
    rda_ctlb = rda_ctlb + rla_ctlb
    rla_ctlb = 0.0_cp
    rda_stateb = 0.0_cp
  END SUBROUTINE adjoint

  !> computes the cost function
  SUBROUTINE cost(rda_ctl, rd_cost)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl
    REAL(cp), INTENT(OUT) :: rd_cost
    !local variables
    REAL(cp) :: rl_ocost, rl_rcost

    CALL debug('', 'In cost: starting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    CALL direct(rda_ctl, rga_state)
    rga_obsgap = rga_state( tg_obs%ia_icoord(1,:) ) - tg_obs%ra_obs

    rl_ocost = SUM(tg_obs%ra_Rm1*rga_obsgap**2)
    rl_rcost = SUM( ABS(rda_ctl) )

    CALL debug( (/rl_ocost, rl_rcost/), '(/rl_ocost, rl_rcost/) = ' )

    rd_cost = rl_ocost/get_obsSize(tg_obs) + tg_ep%r_cs_reg*rl_rcost/SIZE(rda_ctl)
    !saving obsgap
    tg_obs%ra_obsgap = rga_obsgap
    CALL debug(rd_cost, 'In cost: rd_cost = ')
    !CALL debug(tg_obs%ra_obsgap, 'In cost: writing obsgap, tg_obs%ra_obsgap = ')
    CALL write_obsgap(tg_obs)
    !saving exchange parameters
    CALL debug('', 'In cost: ending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
  END SUBROUTINE cost

  !> computes the gradient of the cost function
  SUBROUTINE costb(rda_ctl, rda_ctlb, rd_cost, rd_costb)
    REAL(cp), DIMENSION(:), INTENT(IN)    :: rda_ctl
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_ctlb
    REAL(cp), INTENT(IN)    :: rd_cost
    REAL(cp), INTENT(INOUT) :: rd_costb
    !local variables
    INTEGER  :: ibi
    REAL(cp) :: rl_ocostb, rl_rcostb
    !zeroing local adjoint variables
    rl_ocostb = 0.0_cp
    rl_rcostb = 0.0_cp
    !
    CALL nothing(rd_cost)
    !
    CALL debug('', 'In costb: starting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    !CALL debug(rga_stateb, 'In costb: rga_stateb = ')
    !restauring obsgap
    rga_obsgap = tg_obs%ra_obsgap
    !adjoint
    rl_ocostb = rl_ocostb + rd_costb/get_obsSize(tg_obs)
    rl_rcostb = rl_rcostb + tg_ep%r_cs_reg*rd_costb/SIZE(rda_ctl)
    rd_costb  = 0.0_cp
    DO ibi = 1, SIZE(rda_ctl)
      rda_ctlb(ibi) = rda_ctlb(ibi) + rl_rcostb*SIGN( 1.0_cp, rda_ctl(ibi) )
      !CALL debug( (/rda_ctl(ibi), SIGN( 1.0_cp, rda_ctl(ibi) )/), '(/rda_ctl(ibi), SIGN( 1.0_cp, rda_ctl(ibi) )/) = ' )
    END DO
    rl_rcostb = 0.0_cp
    !DO ibi = 1, SIZE(rga_obsgap)
    !  rga_obsgapb(ibi) = rga_obsgapb(ibi) + tg_obs%ra_Rm1(ibi)*rga_obsgap(ibi)*rl_ocostb
    !END DO
    rga_obsgapb   = rga_obsgapb + tg_obs%ra_Rm1*rga_obsgap*rl_ocostb
    rl_ocostb = 0.0_cp
    !CALL debug(rga_obsgapb, 'In costb: rga_obsgapb = ')
    rga_stateb( tg_obs%ia_icoord(1,:) ) = rga_stateb( tg_obs%ia_icoord(1,:) ) + rga_obsgapb
    rga_obsgapb = 0.0_cp
    !CALL debug(rga_stateb, 'In costb: rga_stateb = ')
    CALL adjoint(rda_ctl, rda_ctlb, rga_state, rga_stateb)
    !CALL debug(rda_ctlb, 'In costb: rda_ctlb = ')
    CALL debug('', 'In costb: ending <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
  END SUBROUTINE costb

  !> @brief Generates and saves a control vector
  !! @param [in,out] td_ep exchange parameter
  !<
  SUBROUTINE make_random_ctl(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    CHARACTER(LEN=IP_ACTION_LEN) :: ala_tmpAction

    ala_tmpAction = td_ep%aa_solver_action!saving the action
    CALL set_default_ctl_param(td_ep)
    CALL debug('In make_random_ctl, after set_default_ctl_param(td_ep)' )
    CALL randval_randidx(td_ep%ra_ctl, td_ep%r_nz_ratio, td_ep%r_max_coef)
    CALL write_ep_data(td_ep, TCTL_DATA, OUTPUT_FILE, 'make_random_ctl' )
    td_ep%aa_solver_action = MAKE_OBS
    CALL make_twin_obs(td_ep)!make and save twin obs and associated model trajectory
    td_ep%aa_solver_action = ala_tmpAction!restauring the action
  END SUBROUTINE make_random_ctl

  !> @brief Generates a zero background control vector
  !! @param [in,out] td_ep exchange parameter
  !<
  SUBROUTINE make_zero_bg(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep

    CALL set_default_ctl_param(td_ep)
    td_ep%ra_ctl          = 0.0_cp
    td_ep%ra_b_ctl        = 0.0_cp
    CALL write_ep_data(tg_ep, BCTL_DATA, OUTPUT_FILE, 'make_zero_bg' )
  END SUBROUTINE make_zero_bg

  !> @brief run direct model and saves the trajectory
  !! @details
  !<
  SUBROUTINE make_model_trajectory(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local variables
    INTEGER :: il_nctl
    CHARACTER(LEN=ip_fnl) :: ala_mtFName
    REAL(cp), DIMENSION(:), ALLOCATABLE :: rla_state

    il_nctl = get_ctlSize(td_ep)
    ALLOCATE( rla_state(il_nctl) )

    SELECT CASE (td_ep%aa_solver_action)
      CASE (MAKE_ADMT)
        ala_mtFName  = make_fileName(td_ep, ADMT_DATA, OUTPUT_FILE)
      CASE (MAKE_DMT, RUN_DIRECT)
        ala_mtFName  = make_fileName(td_ep, DMT_DATA, OUTPUT_FILE)
      CASE (MAKE_TDMT)
        ala_mtFName  = make_fileName(td_ep, TDMT_DATA, OUTPUT_FILE)
    END SELECT
    CALL direct( td_ep%ra_ctl, rla_state )
    CALL write_vector_for_plot(rla_state, ala_mtFName)!saving the model trajectory for plot

  END SUBROUTINE make_model_trajectory

  !> @brief make twin obs, and associated model trajectory
  !! @details file name are build according to the action field
  !<
  SUBROUTINE make_twin_obs(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local variables
    INTEGER :: il_nobs, il_nctl, il_ndim
    CHARACTER(LEN=ip_fnl) :: ala_obsFName, ala_mtFName
    TYPE(obs_structure) :: tl_os
    REAL(cp), DIMENSION(:), ALLOCATABLE :: rla_state
    INTEGER , DIMENSION(:,:), ALLOCATABLE :: ila_idx

    il_nobs = CEILING( td_ep%r_mes_fact*td_ep%r_nz_ratio*get_ctlSize(td_ep) )
    il_nctl = get_ctlSize(td_ep)
    il_ndim = 1!1D problem
    SELECT CASE(td_ep%aa_solver_action)
      CASE (MAKE_OBS)
        ala_obsFName = make_fileName(td_ep, OBS_DATA, OUTPUT_FILE)
        ala_mtFName  = make_fileName(td_ep, DMT_DATA, OUTPUT_FILE)
      CASE (MAKE_TOBS)
        ala_obsFName = make_fileName(td_ep, TOBS_DATA, OUTPUT_FILE)
        ala_mtFName  = make_fileName(td_ep, TDMT_DATA, OUTPUT_FILE)
    END SELECT
    ALLOCATE( rla_state(il_nctl), ila_idx(il_ndim, il_nobs) )

    CALL set_obsSize(tl_os, il_nobs, il_ndim, .TRUE., .FALSE.)
    CALL set_default_obs(tl_os)
    CALL set_obs_fName( tl_os, TRIM(ala_obsFName) )

    CALL direct( td_ep%ra_ctl, rla_state )
    CALL write_vector_for_plot(rla_state, ala_mtFName)!saving the model trajectory for plot

    CALL random_indexes(ila_idx(1,:), il_nctl)
    CALL set_obs_with_icoord(tl_os, rla_state( ila_idx(1,:) ), ila_idx)
    !CALL print_os(tl_os)
    CALL write_obs(tl_os)
    CALL set_obsSize(tl_os, 0, il_ndim, .FALSE., .FALSE.)!free dynamic space (zero size)
  END SUBROUTINE make_twin_obs

  SUBROUTINE make_twin_obs_from_coord(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local variables
    INTEGER :: il_nobs, il_nctl, il_ndim
    TYPE(obs_structure) :: tl_os
    REAL(cp), DIMENSION(:), ALLOCATABLE :: rla_state
    INTEGER , DIMENSION(:,:), ALLOCATABLE :: ila_idx

    CALL set_obs_fName( tl_os, make_fileName(td_ep, OBS_DATA, INPUT_FILE) )
    CALL read_obs(tl_os)

    il_nobs = get_obsSize(tl_os)
    il_nctl = get_ctlSize(td_ep)
    il_ndim = 1!1D problem
    ALLOCATE( rla_state(il_nctl), ila_idx(il_ndim, il_nobs) )
    CALL direct( td_ep%ra_ctl, rla_state )
    !CALL set_default_obs(tl_os)
    CALL set_obs_data(tl_os, rla_state( tl_os%ia_icoord(1,:) ) )
    !CALL print_os(tl_os)
    CALL set_obs_fName( tl_os, make_fileName(td_ep, AOBS_DATA, OUTPUT_FILE) )
    CALL write_obs(tl_os)
    CALL set_obsSize(tl_os, 0, il_ndim, .FALSE., .FALSE.)
  END SUBROUTINE make_twin_obs_from_coord

END PROGRAM cs2_test
