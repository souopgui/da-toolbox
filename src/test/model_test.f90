PROGRAM test
  USE com_tools
IMPLICIT NONE
  REAL(cp), DIMENSION(:), ALLOCATABLE :: rga_traj , rga_obsgap, rga_times
  REAL(cp), DIMENSION(:), ALLOCATABLE :: rga_trajb, rga_obsgapb !adjoint variables
  TYPE(exchange_param) :: tg_ep
  TYPE(obs_structure)  :: tg_obs
  REAL(cp) :: rg_costb

  !main program
  !WRITE(*, *) 'Entering model program*******************'
  CALL read_ep(tg_ep)
  tg_obs%i_nb_dim  = 1
  tg_obs%obs_fileName    = get_obs_fileName (tg_ep)
  tg_obs%obsgap_fileName = get_obsgap_fileName(tg_ep)
  !CALL print_ep(tg_ep)

  SELECT CASE(tg_ep%aa_action)
    CASE (RUN_COST) !cost function recquired
      CALL init_cost()
      CALL cost(tg_ep%ra_ctl, tg_obs%ra_real_coord(:,1), tg_ep%r_cost)
    CASE (RUN_GRADIENT) !gradient of the cost function recquired
      CALL init_grad()
      CALL adjoint_zeroing()
      CALL costb(tg_ep%ra_ctl, tg_ep%ra_grad, tg_obs%ra_real_coord(:,1), tg_ep%r_cost, rg_costb)
    CASE (MAKE_OBS) !make observations for twin experiments
      !CALL init_make_obs
      !CALL make_obs
    CASE (RUN_DIRECT) ! direct model run recquired
      !CALL init_direct()
      !CALL direct(rga_ctl, rga_time, rga_traj)
      WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_action)//'>'
    CASE (RUN_ADJOINT)! adjoint model run recquired
      !CALL init_adjoint()
      !CALL adjoint(...)
      WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_action)//'>'
    CASE DEFAULT
      CALL stop_progam('Unknow action <'//TRIM(tg_ep%aa_action)//'>')
  END SELECT
  CALL write_ep(tg_ep)
  !WRITE(*, *) 'Ending model program  *******************'
  !end of main program
CONTAINS

  !> initializes direct model
  SUBROUTINE init_direct(id_size, rda_times)
    INTEGER, INTENT(IN) :: id_size
    REAL(KIND=cp), DIMENSION(:) :: rda_times
    ALLOCATE( rga_traj(id_size), rga_times(id_size) )
    rga_times = rda_times
    !must also initialize times ...
  END SUBROUTINE init_direct

  !> initializes the environment for the computation of the cost function
  SUBROUTINE init_cost()
    CALL read_obs(tg_obs)
    ALLOCATE( rga_traj(tg_obs%i_nobs), rga_obsgap(tg_obs%i_nobs) )
  END SUBROUTINE init_cost

  !> initializes the environment for the computation of the gradient of the cost function
  !<
  SUBROUTINE init_grad()

    tg_obs%obsgap_fileName = get_obsgap_fileName(tg_ep)
    CALL read_obsgap(tg_obs)
    ALLOCATE( rga_trajb(tg_obs%i_nobs), rga_obsgap(tg_obs%i_nobs), rga_obsgapb(tg_obs%i_nobs) )
    rga_obsgap = tg_obs%ra_obsgap
  END SUBROUTINE init_grad

  !> zeroing adjoint variables
  !<
  SUBROUTINE adjoint_zeroing()
    rga_trajb     = 0.0_cp
    rga_obsgapb   = 0.0_cp
    tg_ep%ra_grad = 0.0_cp
    rg_costb = 1.0_cp
  END SUBROUTINE adjoint_zeroing

  !> direct model for the straight line problem
  SUBROUTINE direct(rda_ctl, rda_time, rda_traj)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl, rda_time
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_traj
    !local variable
    REAL(cp) :: rl_x0, rl_alpha
    INTEGER :: ibi
    CALL debug('','In direct ')
    rl_x0 = rda_ctl(1) + tg_ep%ra_b_ctl(1)
    rl_alpha = rda_ctl(2) + tg_ep%ra_b_ctl(2)
    CALL debug(rl_x0   , '  rl_x0,   = ')
    CALL debug(rl_alpha, '  rl_alpha = ')
    DO ibi = 1, SIZE(rda_traj, 1)
      rda_traj(ibi) = rl_x0 + rl_alpha*rda_time(ibi)
    END DO
  END SUBROUTINE direct

  !> adjoint model for the straight line problem
  SUBROUTINE adjoint(rda_ctl, rda_ctlb, rda_time, rda_traj, rda_trajb)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl, rda_time
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_ctlb
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_traj
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_trajb
    !local variable
    !REAL(cp) :: rl_x0, rl_alpha
    REAL(cp) :: rl_x0b, rl_alphab
    INTEGER :: ibi
    !initialization of local adjoint variables
    rl_x0b = 0.0_cp
    rl_alphab = 0.0_cp
    !recomputation
    !rl_x0 = rda_ctl(1)
    !rl_alpha = rda_ctl(2)
    !adjoint code
    DO ibi = SIZE(rda_trajb, 1), 1, -1
      rl_x0b = rl_x0b + rda_trajb(ibi)
      rl_alphab = rl_alphab + rda_time(ibi)*rda_trajb(ibi)
      rda_trajb(ibi) = 0.0_cp
    END DO
    rda_ctlb(1) = rda_ctlb(1) + rl_x0b
    rda_ctlb(2) = rda_ctlb(1) + rl_alphab
  END SUBROUTINE adjoint

  !> computes the cost function
  SUBROUTINE cost(rda_ctl, rda_time, rd_cost)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl, rda_time
    REAL(cp), INTENT(OUT) :: rd_cost
    CALL debug('', 'In cost: starting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    CALL direct(rda_ctl, rda_time, rga_traj)
    rga_obsgap = rga_traj - tg_obs%ra_obs
    rd_cost = 0.5*SUM(tg_obs%ra_Rm1*rga_obsgap**2)
    !saving obsgap
    tg_obs%ra_obsgap = rga_obsgap
    CALL debug(rd_cost, 'In cost: rd_cost = ')
    CALL debug(tg_obs%ra_obsgap, 'In cost: writing obsgap, tg_obs%ra_obsgap = ')
    CALL write_obsgap(tg_obs)
    !saving exchange parameters
    CALL debug('', 'In cost: ending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

  END SUBROUTINE cost

  !> computes the gradient of the cost function
  SUBROUTINE costb(rda_ctl, rda_ctlb, rda_time, rd_cost, rd_costb)
    REAL(cp), DIMENSION(:), INTENT(IN) :: rda_ctl, rda_time
    REAL(cp), DIMENSION(:), INTENT(INOUT) :: rda_ctlb
    REAL(cp), INTENT(IN) :: rd_cost
    REAL(cp), INTENT(INOUT) :: rd_costb
    !local variables

    CALL debug('', 'In costb: starting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    rga_obsgapb = rga_obsgapb + tg_obs%ra_Rm1*rga_obsgap*rd_costb
    rd_costb = 0.0_cp
    rga_trajb = rga_trajb + rga_obsgapb
    rga_obsgapb = 0.0_cp
    CALL adjoint(rda_ctl, rda_ctlb, rda_time, rga_traj, rga_trajb)
    CALL debug(rda_ctlb, 'In costb: rda_ctlb = ')
    CALL debug('', 'In costb: ending <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
  END SUBROUTINE costb

END PROGRAM test
