!> \file simulator.f90
!! Main program, driver for DA assimilation libraries and numerical models
!<

!> Main program
!<
PROGRAM simul
  USE com_tools
  USE simul_tools
  USE solver_tools
  USE randmod
  USE random
  USE general_tools
  USE cov_vars
IMPLICIT NONE

!variables for the minimization algorithm
!!$  !> brief Name of the subroutine that provide the cost function and its gradient in direct communication
!!$  !! \details In direct communication, simul is the name of the simulator inside m1qn3.
!!$  !! The simulator is the subroutine that computes the value of the (cost) function f and
!!$  !! the value of its gradient g at the current iterate. When m1qn3 needs these values, it executes the instruction
!!$  !! call simulator (...)
!!$  !! In reverse communication, the empty subroutine simul_rc, provided in the standard distribution can be used
!!$  !! EXTERNAL simul_rc
!!$  !<
!!$  EXTERNAL simulator
!!$  !> \brief Calling name of the subroutine that computes the inner product <u, v> of two vectors u and v of Rn.
!!$  !! \details This subroutine is supposed to have the following declaration statement:
!!$  !! subroutine prosca (n, u, v, ps, izs, rzs, dzs).
!!$  !<
!!$  EXTERNAL prosca
!!$  !> \brief Calling name of the subroutine that makes the change of variable
!!$  !! \details This subroutine is used only in DIS mode. It is supposed to have the following declaration statement:
!!$  !!   subroutine ctonb (n, x, y, izs, rzs, dzs).
!!$  !<
!!$  EXTERNAL ctonb
!!$  !> \brief Calling name of the subroutine that does the operation reverse to the one done by ctonb
!!$  !! \details This subroutine is used only in DIS mode. It is supposed to have the following declaration statement:
!!$  !!  subroutine ctcab (n, y, x, izs, rzs, dzs).
!!$  !<
!!$  EXTERNAL ctcab
  !> \brief Positive integer variable that gives the numbers of vector pairs used to approximate the Hessian.
  !! \detail
  !<
  INTEGER :: ig_npairs = 2
  !> \brief Positive integer variable that gives the size of the control vector or the dimension of the problem.
  !! \detail
  !<
  INTEGER :: ig_nctl
  !> \brief Positive integer variable that gives the size of the CTL of the preconditionned problem
  !<
  INTEGER :: ig_npctl

  !> \brief Control vector
  !! \details It must be an array of dimension ig_npctl. On entry, it is supposed to be the initial value.
  !! On return, it is the value of the final point calculated by m1qn3.
  !<
  REAL(KIND=cp), DIMENSION(:), ALLOCATABLE  :: rga_pctl

  !> \brief
  !! \details On entry, it is supposed to be the value of the cost function evaluated at the initial value of the control variable
  !! On return, it is the value of f at the final point.
  !<
  REAL(KIND=cp) :: rg_cost

  !> \brief Array variable for the gradient of the cost function
  !! \details On entry, it is supposed to be the value of the gradient of the cost function at the initial control variable.
  !! On return with omode = 1 (see below), it is the value of the gradient of the cost function at the final point.
  !! For other output modes, its value is undetermined.
  !<
  REAL(KIND=cp), DIMENSION(:), ALLOCATABLE  :: rga_grad

  !> \brief Resolution in the control vector for the l∞ norm (positive)
  !! \details two points whose distance in the sup-norm is less than rg_xmin will be considered as indistinguishable by the optimizer.
  !<
  REAL(KIND=cp) :: rg_xmin

  !> \brief Estimation of the expected decrease in the cost function during the first iteration (positive)
  !! \detail
  !<
  REAL(KIND=cp) :: rg_df1

  !> \brief Stopping criterion that is based on the norm of the gradient of the cost function
  !! \details value range is ]0, 1[
  !!
  !<
  REAL(KIND=cp) :: rg_epsg

  !> \brief norm that is used to test optimality (see the argument rg_epsg).
  !! \details It can be one of the following strings:
  !! - 'two' denotes the Euclidean or l_2 norm
  !! - 'sup' denotes the sup or l_infty norm
  !! - 'dfn' denotes the norm k associated with the scalar product defined in the user-supplied subroutine prosca
  !<
  CHARACTER(LEN=3) :: aga_norm

  !> \brief variable that controls the outputs on channel io
  !! \detail
  !! 0 : No print.
  !! >= 1 Initial and final printouts, error messages.
  !! >= 3 One line of printout per iteration that gives: the index k of the iteration going from the point xk to the point xk+1
  !!      the number of time the simulator has been called, the value f(xk) of the objective function and the directional derivative
  !! >= 4 Print outs from mlis3 during the line-search: see the write-up on mlis3 in modulopt library.
  !! >= 5 Some more printouts at the end of iteration k (see m1qn3 doc for details)
  !<
  INTEGER :: ig_impres

  !> \brief variable that will be taken as the channel number for the outputs
  !! \detail
  !<
  INTEGER :: ig_m1qn3io

  !> \brief  Input mode of m1qn3 that tunes its behavior.
  !! \details Integer array of dimension 3 with the following values are meaningful.
  !! imode(1) determines the scaling mode of m1qn3.
  !!   - M1qn3 will run in DIS (recommended) mode if imode(1) = 0
  !!   - M1qn3 will run in SIS mode if imode(1) = 1
  !! imode(2) specifies the starting mode of m1qn3.
  !!   - A cold start is performed if imode(2) = 0: the first descent direction is then -g1.
  !!   - A warm start is performed if imode(2) = 1: the first descent direction is -g1 (W1 Hessian approx)
  !! imode(3) specifies in direct communication when the simulator has to be called with indic = 1 or
  !!   similarly in reverse communication, when m1qn3 returns to the calling subroutine with indic = 1.
  !!   When imode(3) = 0, the simulator is never called with indic = 1
  !!   When imode(3) > 0, the simulator is called with indic = 1 every imode(3) iteration(s), starting at iteration 1.
  !<
  INTEGER, DIMENSION(3) :: iga_imode

  !> \brief output mode of m1qn3
  !! \details Meaningful values
  !! = 0: The simulator asks to stop by returning the value indic = 0.
    !! = 1: This is the normal way of stopping for m1qn3: the test on the gradient is satisfied (see the meaning of rg_epsg).
  !! = 2: One of the input arguments is not well initialized. This can be:
  !!  - n <= 0, niter <= 0, nsim <= 0, dxmin <= 0.0 or epsg mot in ]0, 1[,
  !!  - ndz < 5n + 1 (in SIS mode) or ndz < 6n + 1 (in DIS mode): not enough storage in memory
  !!  - the contents of iz is not correct for a warm restart,
  !!  - the starting point is almost optimal (the norm of the initial gradient is less than 10−20).
  !! = 3: The line-search is blocked on tmax = 1020 (see section 4.4 and then documentation on mlis3 in modulopt library).
  !! = 4: The maximal number of iterations is reached.
  !! = 5: The maximal number of simulations is reached.
  !! = 6: Stop on dxmin during the line-search (see section 4.4 of m1qn3 doc).
  !! = 7: Either hg, di is nonnegative or hy, si is nonpositive (see section 4.4 of m1qn3 doc).
  !! For additional information and comments, see section 4  of m1qn3 doc.
  !<
  INTEGER :: ig_omode

  !> \brief Maximal number of iterations accepted from m1qn3.
  !! \details m1qn3 uses this variable to return the number of iterations really done
  !<
  INTEGER :: ig_niter

  !> \brief Maximal number of simulations accepted from m1qn3.
  !! \details m1qn3 uses this variable to return the number of simulations really done
  !<
  INTEGER :: ig_nsimul

  !> \brief working array for m1qn3
  !! \detail
  !<
  INTEGER, DIMENSION(5) :: iga_iz

  !> \brief working array for m1qn3
  !! \details this array is of size ndz
  !<
  REAL(KIND = dp), DIMENSION(:), ALLOCATABLE :: rga_dz

  !> \brief Size of the working array rga_dz for m1qn3
  !! \details In SIS mode, m1qn3 needs a working area formed of at least 3 vectors of
  !! dimension n (dk, gk and an auxiliary vector) and it needs for each update (Hessian approx)
  !! one scalar and two vectors. Therefore, if m is the desired number of updates for forming the matrix Wk(Hessian approx),
  !! it is necessary to have: ig_ndz >= 3n +m(2n + 1). m1qn3 based its calculation on this value to determine m.
  !!  if ndz is less than 5n + 1, m1qn3 stops with omode = 2.
  !! In DIS mode, m1qn3 needs an additional vector of dimension n for storing Dk. So, take ig_ndz >= 4n + m(2n + 1) and m >= 1.
  !<
  INTEGER :: ig_ndz

  !> \brief Specifies whether direct or reverse communication is desired for m1qn3
  !! \details In reverse communication, it is used to communicate with the call loop
  !! values :
  !! < 0: implies that m1qn3 will stop immediately using the instruction stop, to prevent entering an infinite call loop in reverse communication,
  !!      due to the fact that the calling program has not left the call loop when m1qn3 returns a negative value;
  !! = 0: indicates that m1qn3 has to work in direct communication;
  !! = 1: indicates that m1qn3 has to work in reverse communication.
  !! it is used by m1qn3 on return to send informations to the calling loop, values
  !! < 0: when m1qn3 has terminated, in which case the call loop must be interrupted;
  !! = 1: the call loop must be pursued.
  !<
  INTEGER :: ig_reverse

  !> \brief Indicates the state of the computation required by m1qn3
  !! \details values
  !! < 0: the computation of f and g required by m1qn3 on its last return was not possible at the given control variable
    !!      this indicates to m1qn3 to adjust the step-size;
  !! = 0: m1qn3 has to stop, for example because some events that m1qn3 cannot understand (not in the field of optimization) has occurred;
  !! > 0: the required computation has been done.
  !! m1qn3 also uses this variable to send information to the calling loop (reverse mode) or simulator (direct mode)
  !! values :
  !! = 1: means that the calling program can do anything except changing the values of indic, n, x, f, and g;
  !!      this value of indic is used by m1qn3 every imode(3) iteration(s), when imode(3) > 0, and never, when imode(3) = 0;
  !! = 4: means that the calling program has to compute f and g, to put them in f and g, and to call back m1qn3.
  !<
  INTEGER :: ig_indic

  !> \brief Integer working array for simulator, prosca, ctonb, and ctcab
  !! \details not used in practice
    !<
  INTEGER, DIMENSION(2)  :: iga_wa !(izs in m1qn3

  !> \brief Real working array for simulator, prosca, ctonb, and ctcab
  !! \details not used in practice
  !<
  REAL(KIND=sp), DIMENSION(2)  :: rga_wa !(rzs in m1qn3)

  !> \brief Real working array for simulator, prosca, ctonb, and ctcab
  !! \details not used in practice
  !<
  REAL(KIND=dp), DIMENSION(2)  :: dga_wa !(dzs in m1qn3)

  CHARACTER(LEN=*), PARAMETER :: ala_namelist = "simulator_namelist"
  !CHARACTER(LEN=ip_snl) :: ala_inputfName
  !setup variables
  LOGICAL :: lm_B_matrix
  CHARACTER(LEN=ip_snl) :: ama_B_fName, ama_Binv_Fname, ama_L_fName, ama_U_fName, ama_Bctl_Fname, ama_Tctl_Fname
  !variables for the test of the gradient
  REAL(cp) :: rm_gradtest_sFactor !> initial factor for the test of the gradient
  INTEGER  :: im_gradtest_nFactor !> number of factor (simulations) for the test of the gradient

  CHARACTER(len=ip_snl):: ala_start, ala_end, date, time
  REAL(cp) :: rl_cpuInitial, rl_cpuFinal
  integer(dp) :: il_countInitial, il_countFinal, il_rate


  Call system_clock(il_countInitial, il_rate)
  CALL CPU_TIME(rl_cpuInitial)
  CALL DATE_AND_TIME(DATE=date, TIME=time)!, ZONE, VALUES
  !DATE has form ccyymmdd. TIME has form hhmmss.sss
  ala_start = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//'--'//time(1:2)//':'//time(3:4)//':'//time(5:6)//time(7:10)

  CALL debug( "In simulator; loading "//TRIM(ala_namelist), tag=dALLWAYS )
  CALL load_namelist(ala_namelist)
  ig_npctl = ig_nctl
  !CALL print_namelist_variables()
  !CALL dpause()
  CALL debug( "In simulator; after loading "//TRIM(ala_namelist), tag=dALLWAYS )
  !!!!!
  CALL debug(tg_ep%aa_simul_action, "In simulator; calling  init_solver for ", tag=dALLWAYS )
  CALL init_solver(tg_ep, ig_nctl, npctl=ig_npctl) !Initialize the solver
  CALL debug("In simulator ---------------------------------- ", tag=dALLWAYS)
  CALL debug(tg_ep%aa_simul_action, "In simulator; after  init_solver for ", tag=dALLWAYS )
  tg_ep%prefix = toLower( tg_ep%aa_simul_action )
  !CALL stop_program("Prescribed stop, this stop is located in the simulator main program")
  !
  SELECT CASE(tg_ep%aa_simul_action)
    CASE (RUN_ASSIM) !data assimilation experiment
      CALL assim_run()
    CASE (RUN_IPF) !Implicit particle filter experiment
      CALL ipf_run(tg_ep%i_ipf_nparticle)!
    CASE (GRAD_TEST) !test of the gradient
      CALL gradient_test()
    CASE (RUN_COST) !cost function recquired
      !CALL read_ep_data(tg_ep, BCTL_DATA, INPUT_FILE )
      WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_solver_action)//'>'
    CASE (RUN_GRADIENT) !gradient of the cost function recquired
      !CALL read_ep_data(tg_ep, BCTL_DATA, INPUT_FILE )
      !WRITE(*, *) 'nothing to do for <'//TRIM(tg_ep%aa_solver_action)//'>'
    CASE (MAKE_OBS) !make observations for twin experiments
      tg_ep%aa_solver_action=MAKE_OBS
      tg_ep%l_save_pdata = .TRUE.
      CALL debug( "In simulator; calling run_solver -> MAKE_OBS", tag=dALLWAYS )
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver   -> MAKE_OBS", tag=dALLWAYS )
    CASE (RUN_DIRECT) ! direct model run recquired
      CALL debug( "In simulator; calling run_solver -> RUN_DIRECT", tag=dALLWAYS )
      tg_ep%aa_solver_action=RUN_DIRECT
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver -> RUN_DIRECT", tag=dALLWAYS )
    CASE (RUN_STARTUP) ! direct model run recquired
      CALL debug( "In simulator; calling run_solver -> STARTUP", tag=dALLWAYS )
      tg_ep%aa_solver_action=RUN_STARTUP
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver -> STARTUP", tag=dALLWAYS )
    CASE (RUN_SPINUP) ! direct model run recquired
      CALL debug( "In simulator; calling run_solver -> SPINUP", tag=dALLWAYS )
      tg_ep%aa_solver_action=RUN_SPINUP
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver -> SPINUP", tag=dALLWAYS )
    CASE (RUN_ENKF) ! direct model run recquired
      CALL debug( "In simulator; calling run_solver -> ENKF", tag=dALLWAYS )
      tg_ep%aa_solver_action=RUN_ENKF
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver -> ENKF", tag=dALLWAYS )
    CASE (RUN_ADJOINT)! adjoint model run recquired
      tg_ep%aa_solver_action=RUN_ADJOINT
      CALL run_solver(tg_ep)
    CASE (MAKE_CTL)! making a control vector
      CALL debug( "In simulator; calling run_solver -> MAKE_CTL", tag=dALLWAYS )
      tg_ep%aa_solver_action=MAKE_CTL
      !CALL print_ep(tg_ep)
      CALL run_solver(tg_ep)
      CALL debug('back to the main program..........')
    CASE (MAKE_BG)! make default background
      tg_ep%aa_solver_action=MAKE_BG
      CALL run_solver(tg_ep)
      CALL debug( "In simulator; after run_solver -> MAKE_BG", tag=dALLWAYS )
!     CASE (RUN_SETUP)! make default background
!       CALL setup()
    CASE (USER_DEFINED)! make default background
      tg_ep%aa_solver_action=USER_DEFINED
      CALL run_solver(tg_ep)
    CASE DEFAULT
      CALL stop_program('Unknown action <'//TRIM(tg_ep%aa_simul_action)//'>')
  END SELECT
  CALL debug(get_ctlSize(tg_ep), "In simulator, get_ctlSize(tg_ep) = ", tag=dALLWAYS )
  CALL debug(get_pctlSize(tg_ep), "In simulator, get_pctlSize(tg_ep) = ", tag=dALLWAYS )
  CALL set_ctlsize(tg_ep, 0, ndim=0, npctl=0)!put this in finalize_solver and replace by call to finilize_solver
  CALL debug( "In simulator; set_ctlsize(tg_ep, 0)", tag=dALLWAYS )
  CALL finalize_solver(tg_ep)

  CALL DATE_AND_TIME(DATE=date, TIME=time)!, ZONE, VALUES
  CALL CPU_TIME(rl_cpuFinal)
  Call system_clock(il_countFinal, il_rate)
  !DATE has form ccyymmdd. TIME has form hhmmss.sss
  ala_end = date(1:4)//'/'//date(5:6)//'/'//date(7:8)//'--'//time(1:2)//':'//time(3:4)//':'//time(5:6)//time(7:10)
  CALL debug(ala_start, 'Starting date--time = ', tag=dALLWAYS)
  CALL debug(ala_end  , 'Ending date--time   = ', tag=dALLWAYS)
  CALL printElapseTime("CPU time:", rl_cpuInitial, rl_cpuFinal)
  CALL printElapseTime("System time:", il_countInitial, il_countFinal, il_rate)
  CALL debug('End of the main program @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@', tag=dALLWAYS)

CONTAINS

  !> run data assimilation experiment
  !! \param[in] clean (optional parameter to) says if the environment should be cleaned or not
  !<
  SUBROUTINE assim_run(clean)
    LOGICAL, INTENT(IN), OPTIONAL :: clean
    !local var
    LOGICAL :: ll_clean

    IF ( PRESENT(clean) ) THEN
      ll_clean = clean
    ELSE
      ll_clean = .FALSE.
    END IF

    OPEN(UNIT=ip_ctl_out  , FILE='rt_ctl_evolution.dat', STATUS='REPLACE',FORM='FORMATTED')
    OPEN(UNIT=ip_grad_out , FILE='rt_grad_evolution.dat',STATUS='REPLACE',FORM='FORMATTED')
    OPEN(UNIT=ip_cost_evol, FILE='rt_cost_evolution.dat', STATUS='REPLACE',FORM='FORMATTED')
    tg_ep%l_save_pdata = .FALSE.
    CALL simulator_init_assim()!allocate global variables and run cost_grad at the initial point

    CALL debug(ig_niter, 'Before call to m1qn3, ig_niter = ', tag=dALLWAYS)
!!$    CALL dpause()
    rg_df1 = rg_cost!*0.75_cp
    SELECT CASE(ig_reverse)
      CASE(M1QN3_REVERSE)
        CALL m1qn3_reverse_driver()
      CASE(M1QN3_DIRECT)
        CALL m1qn3_direct_driver()
    END SELECT
    tg_ep%ra_pctl = rga_pctl !final point
    tg_ep%ra_grad = rga_grad!grad at final point if omode=1, undetermined otherwise
    tg_ep%r_cost  = rg_cost !cost at the final point
    CALL debug(ig_niter, 'After call to m1qn3,ig_niter = ', tag=dALLWAYS)
    CLOSE(ip_ctl_out)
    CLOSE(ip_grad_out)
    CLOSE(ip_cost_evol)
    IF (ll_clean) THEN
      CALL clean_assim()
      CALL debug('', 'After clean_assim ', tag=dALLWAYS)
    END IF

    !Saving the analysed trajectory
    tg_ep%l_restart_from_previous = .FALSE.
    tg_ep%l_save_pdata = .TRUE.
    tg_ep%aa_solver_action = MAKE_ADMT
    CALL run_solver(tg_ep)
    !saving the analysed ctl
    CALL write_ep_data(tg_ep, ACTL_DATA, OUTPUT_FILE, 'Simulator for ANALYSIS after DA', prefix=tg_ep%prefix)
    CALL save_ctl_plot_data(tg_ep, ACTL_DATA)
    !CALL write_ep_data(tg_ep, BCTL_DATA, OUTPUT_FILE, 'Simulator for BACKGROUND after DA', prefix=tg_ep%prefix)
    !CALL save_ctl_plot_data(tg_ep, BCTL_DATA)
    CALL write_ep_data(tg_ep, EP_DATA  , OUTPUT_FILE, 'Simulator for restart after DA', prefix=tg_ep%prefix)
  END SUBROUTINE assim_run

  !> \brief run the implicit particle filter
  !! \param[in] id_nparticle number of particles to be run
  !! It consists in running the data assimilation process, then sample some particles to compute accurate estimate
  !<
  SUBROUTINE ipf_run(id_nparticle)
    INTEGER, INTENT(IN) :: id_nparticle
    !local variables
    REAL(KIND=cp), DIMENSION(ig_npctl, id_nparticle) :: rla_particles!>all the particles
    REAL(KIND=cp), DIMENSION(ig_npctl) :: rla_xs
    REAL(KIND=cp), DIMENSION(ig_npctl) :: rla_xi
    REAL(KIND=cp), DIMENSION(id_nparticle) ::&
      rla_w & !weight
      , rla_nw& !normalized weight
      , rla_lambda!
    INTEGER, DIMENSION(id_nparticle) :: ila_niter
    REAL(KIND=cp) :: rl_w, rl_4DVAR_cost
    INTEGER :: ibi, il_ntry
    LOGICAL :: ll_converged
    CHARACTER(LEN=ip_fnl) :: ipf_w_fileName, ipf_nw_fileName, ipf_l_fileName, ipf_i_fileName
    REAL(KIND=cp), DIMENSION(ig_nctl)      :: rla_4DVAR_ctl
    REAL(KIND=cp), DIMENSION(ig_npctl)     :: rla_4DVAR_pctl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!those variables are used for error analysis
!     CHARACTER(LEN=ip_fnl) :: ala_Tctl_Fname, ipf_err_fileName, assim_err_fileName
!     REAL(KIND=cp), DIMENSION(ig_nctl)      :: rla_tmp, rla_truth, rla_4DVAR_ctl
!     REAL(KIND=cp), DIMENSION(ig_npctl)     :: rla_4DVAR_pctl
!     REAL(cp), DIMENSION(:,:), ALLOCATABLE  :: rla_Tctl_vec
!     REAL(KIND=cp), DIMENSION(id_nparticle) :: rla_error
!     REAL(KIND=cp) :: rl_4DVAR_error
!     INTEGER :: il_Tctl_nrow, il_Tctl_ncol
     INTEGER :: il_np
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rla_w  = 0.0_cp
    rla_nw = 0.0_cp

    !saving the background trajectory
    tg_ep%prefix = 'ipf_b'
    tg_ep%l_save_pdata = .TRUE.
    tg_ep%aa_solver_action = MAKE_ADMT
    CALL run_solver(tg_ep)
    !CALL dpause("after runing background trajectory")

    tg_ep%prefix = 'ipf_assim'
    CALL assim_run(clean=.FALSE.) !do not clean the environment as it is used for implicit particle filter
    CALL debug(100, 'In IPF_run = ', tag=dALLWAYS)

    rla_4DVAR_ctl  = tg_ep%ra_dctl
    tg_ep%ra_dctl = 0.0_cp
    rla_4DVAR_pctl = tg_ep%ra_pctl
    tg_ep%ra_pctl = 0.0_cp

    CALL debug( 200, 'In IPF_run = ', tag=dALLWAYS )

    ig_indic = SIMUL_COSTGRAD
    CALL simulator( ig_indic, ig_npctl, rla_4DVAR_pctl, rl_4DVAR_cost, rga_grad, iga_wa, rga_wa, dga_wa )

    tg_ep%prefix = toLower( tg_ep%aa_simul_action )
    tg_ep%l_save_pdata = .FALSE.
    !CALL dpause()
    CALL init_normal_rand( mu=0.0_cp, sigma=1.0_cp )
    ibi = 1
    il_ntry = 0
    DO WHILE(ibi<=id_nparticle)
      CALL debug(ibi, 'In ipf_run, computing particle ', tag=dALLWAYS)
      !CALL normal_rand(rla_xi)
      CALL rand_normal(rla_xi)
      CALL debug_minmax(rla_xi, 'rla_xi = ', tag=dTRACE)
      !/!\ rla_4DVAR_pctl is the preconditionned control vector resulting from 4DVAR
      CALL sample_and_weight(rla_4DVAR_pctl, rl_4DVAR_cost, rla_xi, rla_xs, rl_w&
            ,ll_converged, rla_lambda(ibi), ila_niter(ibi) &
      )
      IF(ll_converged)THEN
        !storing sample and weight
        rla_particles(:, ibi) = rla_xs
        rla_w(ibi) = rl_w
        ibi = ibi + 1
        il_ntry = 0
      ELSE
        IF( il_ntry<=5 )THEN
          CALL debug(ibi, 'In ipf_run, not converged, resampling particle: ', tag=dALLWAYS)
          il_ntry = il_ntry+1
        ELSE
          CALL debug( il_ntry, 'In ipf_run, to much ressampling: ', tag=dALLWAYS )
          CALL debug( ibi-1, 'In ipf_run, number of succesful particles: ', tag=dALLWAYS )
          CALL stop_program
        END IF
      END IF
    END DO
    !normalization of the weights
    CALL debug(id_nparticle, 'In ipf_run, normalising weights, nb particle ', tag=dALLWAYS)
    CALL normalizeWeiths(rla_w, rla_nw)
    CALL debug(rla_w , '  non normalised weights: ', tag=dALLWAYS)
    CALL debug(rla_nw, '  normalised weights    : ', tag=dALLWAYS)
    CALL debug(SUM(rla_nw), '  sum of normalised weights    : ', tag=dALLWAYS)

    !Saving particles
    DO ibi = 1, id_nparticle
      tg_ep%ra_pctl = rla_particles(:, ibi)
      CALL p2ctl(tg_ep)
      WRITE(tg_ep%prefix, FMT='(A,I4.4)') 'particle', ibi
      CALL save_ctl_plot_data(tg_ep, ACTL_DATA)
    END DO

    !computation of the weighted average
    tg_ep%ra_pctl = 0.0_cp
    CALL debug('In ipf_run, computing weighted average ', tag=dALLWAYS)
    DO ibi = 1, id_nparticle
      tg_ep%ra_pctl = tg_ep%ra_pctl + rla_nw(ibi) * rla_particles(:, ibi)
    END DO

    CALL debug('', 'In ipf_run, computation of particles completed ', tag=dALLWAYS)
    tg_ep%l_save_pdata = .TRUE.
    tg_ep%aa_solver_action = MAKE_ADMT
    tg_ep%prefix = toLower( tg_ep%aa_simul_action )
    CALL run_solver(tg_ep)
    CALL write_ep_data(tg_ep, ACTL_DATA, OUTPUT_FILE, 'Simulator for ANALYSIS after IPF', prefix=tg_ep%prefix)
    CALL save_ctl_plot_data(tg_ep, ACTL_DATA)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !saving  weights
    ipf_w_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='ipf_weights.dat')
    CALL write_vector_for_plot(rla_w, ipf_w_fileName)
    ipf_nw_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='ipf_normalized_weights.dat')
    CALL write_vector_for_plot(rla_nw, ipf_nw_fileName)
    !saving  lambda
    ipf_l_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='ipf_lambda.dat')
    CALL write_vector_for_plot(rla_lambda, ipf_l_fileName)
    !saving  niter
    ipf_i_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='ipf_niters.dat')
    CALL write_vector_for_plot(REAL(ila_niter,cp), ipf_i_fileName)
    !**************
    !This part of the code is used for error analysis only
!     !**************
!     ala_Tctl_Fname = make_fileName(OTHER_DATA, INPUT_FILE, basename=ama_Tctl_Fname)
!     !reading the true external control vector
!     CALL readInfo(ala_Tctl_Fname, il_Tctl_nrow, il_Tctl_ncol)
!     IF( ((il_Tctl_nrow/=1).AND.(il_Tctl_ncol/=1)).OR.(il_Tctl_nrow*il_Tctl_ncol/=ig_nctl) )THEN
!       CALL debug(ala_Tctl_Fname, 'Bad size for the external control vector in')
!       CALL debug( 'It should be a 1D vector of size > 0 ',tag=dALLWAYS )
!       CALL debug( (/il_Tctl_nrow, il_Tctl_ncol/), 'Now it is an array of SHAPE: ', tag=dALLWAYS )
!       CALL stop_program( )
!     ELSE
!       CALL debug(ig_nctl, 'In ipf_run: the size of the control vector is: ', tag=dALLWAYS)
!       ALLOCATE( rla_Tctl_vec(il_Tctl_nrow, il_Tctl_ncol) )
!       CALL readMatrix(ala_Tctl_Fname, rla_Tctl_vec)
!       rla_truth = RESHAPE(rla_Tctl_vec, (/ig_nctl/))
!       DEALLOCATE( rla_Tctl_vec )
!     END IF
!     !End of reading of the true external control vector
!
!     rla_tmp = rla_4DVAR_ctl - rla_truth
!     rl_4DVAR_error = SQRT( DOT_PRODUCT(rla_tmp,rla_tmp) )
!     !incremental weighted average and associated error
!     DO ibk =1, id_nparticle
!       !CALL debug(ibk, 'In ipf_run, incremental weighted average ', tag=dALLWAYS)
!       CALL normalizeWeiths( rla_w(1:ibk), rla_nw(1:ibk) )
!       !CALL debug(SUM(rla_nw(1:ibk) ), '  sum of normalised weights    : ', tag=dALLWAYS)
!       tg_ep%ra_pctl = 0.0_cp
!       DO ibi = 1, ibk
!         tg_ep%ra_pctl = tg_ep%ra_pctl + rla_nw(ibi) * rla_particles(:, ibi)
!       END DO
!       CALL p2ctl(tg_ep)
!       rla_tmp = tg_ep%ra_b_ctl + tg_ep%ra_dctl - rla_truth
!       rla_error(ibk) = SQRT( DOT_PRODUCT(rla_tmp,rla_tmp) )
!     END DO
!     ipf_err_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='ipf_error.dat')
!     CALL write_vector_for_plot(rla_error, ipf_err_fileName)
!     assim_err_fileName = make_fileName(OTHER_DATA, OUTPUT_FILE, basename='4D-VAR_error.dat')
!     CALL write_vector_for_plot((/rl_4DVAR_error/), assim_err_fileName)
!     saving the error
!
!     CALL debug(rla_error,      'rla_error      = ', tag=dALLWAYS)
!     CALL debug(rl_4DVAR_error, 'rl_4DVAR_error = ', tag=dALLWAYS)
    !!! Analysis for selected numbers of particles, multiple of 5
    DO il_np=5,id_nparticle,5
      CALL debug(il_np, 'In ipf_run, Analysis for selected number of particles ', tag=dALLWAYS)
      CALL normalizeWeiths( rla_w(1:il_np), rla_nw(1:il_np) )
      CALL debug(SUM(rla_nw(1:il_np) ), '  sum of normalised weights    : ', tag=dALLWAYS)
      tg_ep%ra_pctl = 0.0_cp
      DO ibi = 1, il_np
        tg_ep%ra_pctl = tg_ep%ra_pctl + rla_nw(ibi) * rla_particles(:, ibi)
      END DO
      WRITE(tg_ep%prefix, FMT='(A,I4.4)') 'ipf', il_np
      CALL p2ctl(tg_ep)
!       CALL write_ep_data(tg_ep, ACTL_DATA, OUTPUT_FILE, 'Simulator for ANALYSIS after IPF', prefix=tg_ep%prefix)

      CALL save_ctl(tg_ep, ACTL_DATA, OUTPUT_FILE, 'Simulator for ANALYSIS after IPF', prefix=tg_ep%prefix)
      CALL save_ctl_plot_data(tg_ep, ACTL_DATA)
    END DO
    !!!End of the part of the code that is used for error analysis only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL clean_assim()
  END SUBROUTINE ipf_run

  !> \brief Computes a sample and the associated weight for implicit particle filte
  !! \param[in] rda_mu analazed displacement from the background control vector, result od the 4D-VAR process
  !! \param[in] rd_phi value of the cost function at rda_mu
  !! \param[in] rda_xi reference variable for IPF
  !! \param[out] rda_xs computed sample
  !! \param[out] rd_w weight associated with the sample rda_xs
  !! \param[out] ld_converged Says if the minimization has converged or not
  !! \param[out] rd_lambda parameter of the newton iteration
  !! \param [out] id_niter number of iterations run before convergence
  !<
  SUBROUTINE sample_and_weight( rda_mu, rd_phi, rda_xi, rda_xs, rd_w, ld_converged, rd_lambda, id_niter )
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_mu, rda_xi
    REAL(KIND=cp), INTENT(IN) :: rd_phi
    REAL(KIND=cp), DIMENSION( SIZE(rda_mu) ), INTENT(OUT) :: rda_xs
    REAL(KIND=cp), INTENT(OUT) :: rd_w
    LOGICAL, INTENT(OUT) :: ld_converged
    REAL(KIND=cp), INTENT(OUT) :: rd_lambda
    INTEGER, INTENT(OUT) :: id_niter
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_mu) ) :: rla_eta, rla_grad
    REAL(KIND=cp) :: rl_cost, rl_funcG, rl_lambda, rl_rand, rl_rho, rl_dldrho!,rl_g
    INTEGER :: il_iter, il_nxi

    CALL debug( 0, "In sample_and_weight ", tag=dALLWAYS )
    il_nxi = SIZE(rda_xi)
    rd_w = 1.0_cp
    rl_rho = DOT_PRODUCT(rda_xi, rda_xi)
    rla_eta = rda_xi/SQRT(rl_rho)
    !initialization of the Newton iteration
    CALL RANDOM_NUMBER(rl_rand)
    rl_lambda = SQRT(rl_rand)*1d-1
    rda_xs = rda_mu + rl_lambda*rla_eta
    ig_indic = SIMUL_COSTGRAD
    CALL simulator(ig_indic, ig_npctl, rda_xs, rl_cost, rla_grad, iga_wa, rga_wa, dga_wa )
    il_iter = 0
    ld_converged = .FALSE.
    rl_funcG = rl_cost - rd_phi - 0.5_cp*rl_rho
    !========
      CALL debug_minmax(rda_xi, "rda_xi = ", tag=dALLWAYS)
      CALL debug_minmax(rla_eta, "rla_eta = ", tag=dALLWAYS)
      CALL debug_minmax(rda_mu, "rda_mu = ", tag=dALLWAYS)
      CALL debug_minmax(rla_grad, "rla_grad = ", tag=dALLWAYS)
      CALL debug_minmax(rda_xs, "rda_xs = ", tag=dALLWAYS)
      CALL debug(rl_lambda, "rl_lambda = ", tag=dALLWAYS)
      CALL debug(rl_funcG, "rl_funcG = ", tag=dALLWAYS)
      !CALL debug( (/rl_cost, rd_phi, rl_rho/), '(/rl_cost, rd_phi, rl_rho/)', tag=dALLWAYS )
      PRINT*,"rl_cost = ",rl_cost
      PRINT*,"rd_phi  = ",rd_phi
      PRINT*,"rl_rho  = ",rl_rho
      !CALL dpause()

    !=========
    DO WHILE( (.NOT.ld_converged).AND.(il_iter<tg_ep%i_ipf_maxiter) )
      il_iter = il_iter + 1
      rl_lambda = rl_lambda - rl_funcG/( DOT_PRODUCT(rla_grad, rla_eta) )
      !CALL debug(il_iter, '  In sample_and_weight, loop ', tag=dALLWAYS)
      rda_xs = rda_mu + rl_lambda*rla_eta
      CALL debug_minmax( rla_grad, "rla_grad = ", tag=dALLWAYS )
      CALL debug_minmax( rda_xs, "rda_xs = ", tag=dALLWAYS )
      CALL debug( rl_lambda, "rl_lambda = ", tag=dALLWAYS )
      CALL debug( rl_funcG, "rl_funcG = ", tag=dALLWAYS )
      !CALL dpause(  )
      ig_indic = SIMUL_COSTGRAD
      CALL simulator( ig_indic, ig_npctl, rda_xs, rl_cost, rla_grad, iga_wa, rga_wa, dga_wa )
      SELECT CASE ( ig_indic )
        CASE ( SIMUL_IMPOSSIBLE )!impossible to compute the cost function and/or its gradient at the given ctl
          CALL debug('', ' In sample_and_weight: J and/or its gradJ impossible at the given point', tag=dALLWAYS)
          !CALL dpause( 'Forcing il_iter to tg_ep%i_ipf_maxiter' )
          il_iter = tg_ep%i_ipf_maxiter!force the stopping criteria with failure(ld_converged is FALSE by default)
        CASE DEFAULT
          rl_funcG = rl_cost - rd_phi - 0.5_cp*rl_rho
          IF(rl_funcG < tg_ep%r_ipf_tol)THEN !converges
            ld_converged = .TRUE.
            CALL debug(il_iter, '  In sample_and_weight, converged after ', tag=dALLWAYS)
            !compute the weight of the particle
            rl_dldrho = 1.0_cp/DOT_PRODUCT(rla_grad, rla_eta)
            rd_w = (1- (REAL(il_nxi,cp)/2.0_cp) )*log(rl_rho) + REAL(il_nxi-1, cp)*log(ABS(rl_lambda))&
            +log(ABS(rl_dldrho))
          ELSE
            !rl_lambda = rl_lambda - rl_funcG/( DOT_PRODUCT(rla_grad, rla_eta) )
          END IF
      END SELECT
      !CALL dpause()
    END DO
    rd_lambda = rl_lambda
    id_niter = il_iter
  END SUBROUTINE sample_and_weight

  !> \brief normalizes weight for implicit particle filter
  !! \param[in] rda_w non normalized weight
  !! \param[out] rda_nw normalized weight
  !<
  SUBROUTINE normalizeWeiths(rda_w, rda_nw)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_w
    REAL(KIND=cp), DIMENSION( SIZE(rda_w) ), INTENT(OUT) :: rda_nw
    REAL(KIND=cp) :: rl_wmin, rl_nwsum
    REAL(KIND=cp), DIMENSION( SIZE(rda_w) ) :: rla_tmp
    !INTEGER :: ibi

    rl_wmin  = MINVAL(rda_w)
    rla_tmp = rda_w - rl_wmin
    rda_nw   = EXP( rla_tmp )
    rl_nwsum = SUM(rda_nw)
    rda_nw   = rda_nw/rl_nwsum
  END SUBROUTINE normalizeWeiths

  !> Computes the optimal initial condition in direct mode, this is the system state at initial time
  !! \details  here, this routine is used as the driver for the minimizer algorithm
  !! It runs m1qn3 in direct  mode
  !<
  SUBROUTINE m1qn3_direct_driver()

    CALL debug('', 'In m1qn3_direct_driver', tag=dALLWAYS)
    ig_reverse = M1QN3_DIRECT
    CALL m1qn3(simulator, prosca, ctonb, ctcab, ig_npctl, rga_pctl, rg_cost, rga_grad,&
            rg_xmin, rg_df1, rg_epsg, aga_norm, ig_impres, ig_m1qn3io,&
            iga_imode, ig_omode, ig_niter, ig_nsimul, iga_iz, rga_dz, ig_ndz,&
            ig_reverse, ig_indic, iga_wa, rga_wa, dga_wa &
      )

    CALL m1qn3_print_stat(ig_omode, ig_niter, ig_nsimul)
  END SUBROUTINE m1qn3_direct_driver


  !> Computes the optimal initial condition in reverse mode, this is the system state at initial time
  !!
  !! \details here this routine is used as the driver for the minimizer algorithm
  !! It runs m1qn3 in reverse mode
  !<
  SUBROUTINE m1qn3_reverse_driver ()
    IMPLICIT NONE
    !local variable

    !initialization
    !CALL make_twin_obs(u, nlocal, ne_local, t_local, scl, scl_fltwt, iter)

    !1 - computing or loading the initial control variable

    !2 - computing the cost function and its gradient at the initial control variable
    !CALL simulator_init_assim

    CALL debug('', 'In m1qn3_reverse_driver', tag=dALLWAYS)

    !Setting the communication mode
    ig_reverse = M1QN3_REVERSE !set reverse communication on
    DO WHILE (ig_reverse == M1QN3_REVERSE)
        CALL m1qn3(&
            simulator, prosca, ctonb, ctcab, ig_npctl, rga_pctl, rg_cost, rga_grad,&
            rg_xmin, rg_df1, rg_epsg, aga_norm, ig_impres, ig_m1qn3io,&
            iga_imode, ig_omode, ig_niter, ig_nsimul, iga_iz, rga_dz, ig_ndz,&
            ig_reverse, ig_indic, iga_wa, rga_wa, dga_wa &
        )
        IF (ig_reverse == M1QN3_REVERSE) THEN
              CALL simulator(ig_indic, ig_npctl, rga_pctl, rg_cost, rga_grad, iga_wa, rga_wa, dga_wa )
        END IF
    END DO
    CALL m1qn3_print_stat(ig_omode, ig_niter, ig_nsimul)
    !post processing
  END SUBROUTINE m1qn3_reverse_driver

  !> \brief gradient test
  !!
  !<
  SUBROUTINE gradient_test()
    REAL(KIND=cp), PARAMETER :: rp_fact  = 2.0_cp
    REAL(KIND=cp), DIMENSION(im_gradtest_nFactor) :: rla_alpha, rla_cost, rla_test
    REAL(KIND=cp) :: rl_alpha, rl_cost0, rl_norm0
    INTEGER :: ibi

    CALL simulator_init_assim()
    rl_norm0 = DOT_PRODUCT(rga_grad, rga_grad)
    rl_cost0 = rg_cost
    OPEN (UNIT=ip_grad_test, STATUS='REPLACE', FILE=TRIM('grad_test.dat'), ACCESS='SEQUENTIAL', ACTION='WRITE')
    OPEN (UNIT=28, STATUS='REPLACE', FILE=TRIM('rt_diagnostic.dat'), ACCESS='SEQUENTIAL', ACTION='WRITE')

    ig_indic = SIMUL_COSTONLY
    rl_alpha = rm_gradtest_sFactor
    DO ibi = 1, im_gradtest_nFactor
      !PRINT*, 'gradient test : loop ', ibi
      rga_pctl = -rl_alpha*rga_grad
      rg_cost = -999.0_cp
      CALL simulator(ig_indic, ig_npctl, rga_pctl, rg_cost, rga_grad, iga_wa, rga_wa, dga_wa )
      rla_cost(ibi)  = rg_cost
      rla_alpha(ibi) = rl_alpha
      rla_test(ibi) = ABS( rla_cost(ibi) - rl_cost0 ) / ( rla_alpha(ibi)*rl_norm0 )
      WRITE(*, *) "rla_alpha(ibi), rla_test(ibi) = ",rla_alpha(ibi), rla_test(ibi)
      WRITE(*, *) "rl_cost0, rl_norm0            = ",rl_cost0, rl_norm0
      WRITE(ip_grad_test, *) rla_alpha(ibi), rla_test(ibi), 1.0_cp
      WRITE(28, *) rla_cost(ibi), rl_cost0, rl_norm0
      rl_alpha = rp_fact*rl_alpha
    END DO
    CLOSE(ip_grad_test)
  END SUBROUTINE gradient_test

  !> \brief Initializes variables for data assimilation
  !<
  SUBROUTINE simulator_init_assim()
    IMPLICIT NONE
    !INTEGER, INTENT(IN) :: id_nctl
    !LOGICAL, DIMENSION(3) :: lla_param
    !REAL(cp), DIMENSION(:, :), ALLOCATABLE :: rla_obs
    !INTEGER :: il_nbCol

    CALL debug('', 'Entering simulator_init_assim', tag=dALLWAYS)
    IF(ig_npctl<=0)THEN
        CALL stop_program(' In simulator_init_assim : zero or negative size control vector')
    END IF
    ig_ndz = 4*ig_npctl + ig_npairs*( 2*ig_npctl + 1 )
    CALL debug('', ' In simulator_init_assim : allocating space for optimization process', tag=dALLWAYS)
    ALLOCATE(&
          rga_pctl (ig_npctl),&
          rga_grad(ig_npctl),&
          rga_dz  (ig_ndz )&
    )
    rga_dz = 0.0_cp

    rga_pctl = tg_ep%ra_pctl

    CALL init_simul_counter()

    !initialization of the minimizer
    ig_indic = SIMUL_COSTGRAD
    CALL simulator(ig_indic, ig_npctl, rga_pctl, rg_cost, rga_grad, iga_wa, rga_wa, dga_wa )
    CALL debug('', 'After the initialization of the minimizer', tag=dALLWAYS)
    !CALL dpause()
    SELECT CASE (ig_indic)
    CASE (SIMUL_IMPOSSIBLE)!impossible to compute the cost function and/or its gradient at the given ctl
      CALL stop_program(' Impossible to compute the cost function and/or its gradient at the initial point')
    CASE (SIMUL_ASK_STOP)!something that can not be handle by the minimizer happened
      CALL stop_program(' Something that can not be handled by the minimizer happened')
    CASE DEFAULT
    !
    END SELECT
    CALL debug('', 'Exiting simulator_init_assim ', tag=dALLWAYS)
  END SUBROUTINE simulator_init_assim

  SUBROUTINE clean_assim()
    CALL debug('', ' In clean_assim : deallocating assim variables', tag=dALLWAYS)
    DEALLOCATE(rga_pctl, rga_grad, rga_dz)
  END SUBROUTINE clean_assim

!   !> \brief setup the true CTL and the background covanriance matrix and its inverse from external data
!   !!
!   !<
!   SUBROUTINE setup()
!     !setup variables
! !     LOGICAL :: ll_B_matrix
!     CHARACTER(LEN=ip_fnl) :: ala_B_fName, ala_Binv_Fname, ala_Bctl_Fname, ala_Tctl_Fname
!     INTEGER :: il_B_nrow, il_B_ncol, il_Binv_nrow, il_Binv_ncol, il_Bctl_nrow, il_Bctl_ncol, il_Tctl_nrow, il_Tctl_ncol, il_nctl
!     REAL(cp), DIMENSION(:,:), ALLOCATABLE :: rla_B_vec, rla_Bctl_vec, rla_Tctl_vec!, rla_Binv_vec
!
!     !Reset the dynamically allocated fields of tg_ep
!     CALL set_ctlsize(tg_ep, 0, ndim=0, nchanged=0)
!     !builds the file names
!     ala_B_fName    = make_fileName(OTHER_DATA, INPUT_FILE, basename=ama_B_fName   )
!     ala_Binv_Fname = make_fileName(OTHER_DATA, INPUT_FILE, basename=ama_Binv_Fname)
!     ala_Bctl_Fname = make_fileName(OTHER_DATA, INPUT_FILE, basename=ama_Bctl_Fname)
!     ala_Tctl_Fname = make_fileName(OTHER_DATA, INPUT_FILE, basename=ama_Tctl_Fname)
!
!     !read the information from external data files
!
!     !processing the true external control vector
!     CALL readInfo(ala_Tctl_Fname, il_Tctl_nrow, il_Tctl_ncol)
!     il_nctl = MAX(il_Tctl_nrow, il_Tctl_ncol)
!     tg_ep%l_B_matrix = lm_B_matrix
!     CALL set_ctlsize(tg_ep, il_nctl)
!
!     IF( ((il_Tctl_nrow/=1).AND.(il_Tctl_ncol/=1)).OR.(il_nctl<0) )THEN
!       CALL debug(ala_Tctl_Fname, 'Bad size for the external control vector in')
!       CALL debug( 'It should be a 1D vector of size > 0 ',tag=dALLWAYS )
!       CALL debug( (/il_Tctl_nrow, il_Tctl_ncol/), 'Now it is an array of size: ',tag=dALLWAYS )
!       CALL stop_program( )
!     ELSE
!       CALL debug(il_nctl, 'In setup: the size of the control vector is: ',tag=dALLWAYS)
!       ALLOCATE( rla_Tctl_vec(il_Tctl_nrow, il_Tctl_ncol) )
!       CALL readMatrix(ala_Tctl_Fname, rla_Tctl_vec)
!     END IF
!
!     !processing the background external control vector
!     CALL readInfo(ala_Bctl_Fname, il_Bctl_nrow, il_Bctl_ncol)
!
!     IF( ((il_Bctl_nrow/=1).AND.(il_Bctl_ncol/=1)).OR.(il_Bctl_nrow*il_Bctl_ncol/=il_nctl) )THEN
!       CALL debug(ala_Bctl_Fname, 'Bad size for the external control vector in')
!       CALL debug( 'It should be a 1D vector of size > 0 ',tag=dALLWAYS )
!       CALL debug( (/il_Bctl_nrow, il_Bctl_ncol/), 'Now it is an array of size: ', tag=dALLWAYS )
!       CALL stop_program( )
!     ELSE
!       CALL debug(il_nctl, 'In setup: the size of the control vector is: ', tag=dALLWAYS)
!       ALLOCATE( rla_Bctl_vec(il_Bctl_nrow, il_Bctl_ncol) )
!       CALL readMatrix(ala_Bctl_Fname, rla_Bctl_vec)
!     END IF
!
!     !processing external background covariance matrix
!     CALL readInfo(ala_B_fName, il_B_nrow, il_B_ncol)
!     IF( lm_B_matrix )THEN
!       IF( (il_B_nrow/=il_B_ncol).OR.(il_B_nrow/=il_nctl) )THEN
!         CALL debug('Bad size for the external cov matrix in '//TRIM(ala_B_fName),tag=dALLWAYS )
!         CALL debug( (/il_nctl,il_nctl/),'It should be a square matrix of size: ',tag=dALLWAYS )
!         CALL debug( (/il_B_nrow, il_B_ncol/),'Now it is a matrix of size: ',tag=dALLWAYS )
!         CALL stop_program( )
!       ELSE
!         CALL debug('In setup: full background covariance matrix: ',tag=dALLWAYS)
!         CALL readMatrix(ala_B_fName, tg_ep%ra_B_mat)
!       END IF
!     ELSE
!       IF( ((il_B_ncol/=1).AND.(il_B_nrow/=1)).OR.(il_B_ncol*il_B_nrow/=il_nctl) )THEN
!         CALL debug('Bad size for the external cov matrix in '//TRIM(ala_B_fName),tag=dALLWAYS )
!         CALL debug( (/il_nctl/),'It should be a square diag matrix store as a vector of size: ',tag=dALLWAYS )
!         CALL debug( (/il_B_nrow, il_B_ncol/),'Now it is a matrix of size: ',tag=dALLWAYS )
!         CALL stop_program( )
!       ELSE
!         CALL debug('In setup: diagonal background covariance matrix stored as a vector: ',tag=dALLWAYS)
!         ALLOCATE( rla_B_vec(il_B_nrow, il_B_ncol) )
!         CALL readMatrix(ala_B_fName, rla_B_vec)
!         tg_ep%ra_sigma_ctl = RESHAPE(rla_B_vec, (/il_nctl/))
!         DEALLOCATE( rla_B_vec )
!       END IF
!     END IF
!     !
!     !processing external inverse background covariance matrix
!     CALL readInfo(ala_Binv_Fname, il_Binv_nrow, il_Binv_ncol)
!     IF( lm_B_matrix )THEN
!       IF( (il_Binv_nrow/=il_Binv_ncol).OR.(il_Binv_nrow/=il_nctl) )THEN
!         CALL debug('Bad size for the external cov matrix in '//TRIM(ala_Binv_Fname),tag=dALLWAYS )
!         CALL debug( (/il_nctl,il_nctl/), 'It should be a square matrix of size : ',tag=dALLWAYS )
!         CALL debug( (/il_Binv_nrow, il_Binv_ncol/), 'Now it is a matrix of size: ',tag=dALLWAYS )
!         CALL stop_program( )
!       ELSE
!         CALL debug('In setup: full inverse background covariance matrix: ',tag=dALLWAYS)
!         CALL readMatrix(ala_Binv_fName, tg_ep%ra_Binv_mat)
!         CALL debug(' reading inverse background covariance matrix done',tag=dALLWAYS)
!       END IF
!
!     ELSE
!       CALL debug('In setup: diagonal inverse background covariance matrix stored as a vector: ', tag=dALLWAYS)
!       CALL debug('  No need to read it.', tag=dALLWAYS)
! !       IF( ((il_Binv_ncol/1).AND.(il_Binv_nrow/=1)).OR.(il_Binv_ncol*il_Binv_nrow/=il_nctl) )THEN
! !         CALL debug('Bad size for the external cov matrix in '//TRIM(ala_Binv_Fname), tag=dALLWAYS )
! !         CALL debug( (/il_nctl/), 'It should be a square diag matrix store as a vector of size: ', tag=dALLWAYS )
! !         CALL debug( (/il_Binv_nrow, il_Binv_ncol/), 'Now it is a matrix of size: ', tag=dALLWAYS )
! !         CALL stop_program( )
! !       ELSE
! !         CALL debug('In setup: diagonal inverse background covariance matrix stored as a vector: ', tag=dALLWAYS)
! !         ALLOCATE( rla_Binv_vec(il_Binv_nrow, il_Binv_ncol) )
! !         CALL readMatrix(ala_Binv_fName, rla_Binv_vec)
! !         tg_ep%ra_Bm1 = RESHAPE(rla_Binv_vec, (/il_nctl/))
! !       END IF
!     END IF
!
!     !saving data
!     tg_ep%ra_b_ctl = 0.0_cp
!     tg_ep%ra_dctl   = RESHAPE(rla_Tctl_vec, (/il_nctl/))
!     CALL write_ep_data(tg_ep, CTL_DATA, OUTPUT_FILE,  'Simulator for DIRECt and OBS from setup')
!     CALL save_ctl_plot_data(tg_ep, CTL_DATA)
!
!     tg_ep%ra_dctl = 0.0_cp
!     tg_ep%ra_b_ctl = RESHAPE(rla_Bctl_vec, (/il_nctl/))
!     CALL write_ep_data(tg_ep, BCTL_DATA, OUTPUT_FILE, 'Simulator for DA from setup')
!     CALL save_ctl_plot_data(tg_ep, BCTL_DATA)
!
!     DEALLOCATE( rla_Tctl_vec )
!     DEALLOCATE( rla_Bctl_vec )
!     CALL debug( ' Exiting setup',tag=dALLWAYS )
!   END SUBROUTINE setup

  SUBROUTINE load_namelist(fName)
    INTEGER, PARAMETER           :: ip_numnam = 68
    CHARACTER(LEN=*), INTENT(IN) :: fName
    !INTEGER                      :: il_obs_level, il_nobs_x, il_nobs_t

    !REAL(cp) :: rl_sigmaR
    !REAL(cp) :: rl_max_coef, rl_nz_ratio, rl_mes_fact, rl_cs_reg!CS parameters
    !REAL(cp) :: rl_wGD, rl_wb, rl_wGrad, rl_sigmaB !regularization weighting parameters
    !implicit particle filte
    INTEGER  :: il_ipf_nparticle, il_ipf_maxiter
    REAL(cp) :: rl_ipf_tol
    REAL(cp) :: starting_factor
    INTEGER  :: nb_factor
    !
    LOGICAL :: ll_restart_from_previous, ll_run_from_ctl!, ll_amplitude, ll_location, ll_sigma, ll_useGD
    !LOGICAL :: debug_default, debug_io, debug_memory, debug_netcdf, debug_result, debug_trace
    CHARACTER(LEN=ip_snl) :: ctl_bName, obs_bName, dmt_bName, ims_bName, bctl_bName, ogap_bName, input_dir, output_dir, ala_reverse
    CHARACTER(LEN=ip_snl) :: ala_action, ala_solver_path
    !setup variables
    LOGICAL :: ll_B_matrix
    CHARACTER(LEN=ip_snl) :: B_fName, Binv_Fname, L_fName, U_fName, B_ctl_Fname, T_ctl_Fname

    NAMELIST/NAM_general/&
      ala_action,&
      ig_nctl   ,&
      input_dir ,&
      output_dir,&
      ctl_bName ,&
      obs_bName ,&
      dmt_bName ,&
      ims_bName ,&
      bctl_bName,&
      ogap_bName,&
      ala_solver_path,&
      ll_run_from_ctl,&
      lg_solve_all_together,&
      ll_restart_from_previous

    NAMELIST/NAM_m1qn3/&
      rg_xmin   ,&
      rg_epsg   ,&
      aga_norm  ,&
      ig_impres ,&
      ig_m1qn3io,&
      iga_imode ,&
      ig_niter  ,&
      ig_nsimul ,&
      ig_npairs ,&
      ala_reverse


    NAMELIST/NAM_setup/&
      ll_B_matrix ,&
      B_fName     ,&
      Binv_Fname  ,&
      L_fName     ,&
      U_fName     ,&
      B_ctl_Fname ,&
      T_ctl_Fname

    NAMELIST/NAM_GRADTEST/&
      starting_factor,&
      nb_factor

    NAMELIST/NAM_IPF/&
      rl_ipf_tol      ,&
      il_ipf_nparticle,&
      il_ipf_maxiter

!     NAMELIST/NAM_CS/&
!       rl_max_coef,&
!       rl_nz_ratio,&
!       rl_mes_fact,&
!       rl_cs_reg

    OPEN(ip_numnam, FILE=fName, FORM='FORMATTED', STATUS='OLD')
    !If one is going to read many blocs in the same file, it is recommended to rewind
    REWIND(ip_numnam)
    READ(ip_numnam, NAM_general)!reading the block NAM_general
    REWIND(ip_numnam)
    READ(ip_numnam, NAM_m1qn3  )!reading the block NAM_m1qn3
    REWIND(ip_numnam)
    READ(ip_numnam, NAM_setup  )!reading the block NAM_setup

    REWIND(ip_numnam)
    READ(ip_numnam, NAM_GRADTEST)!reading the block NAM_GRADTEST

    !REWIND(ip_numnam)
    !READ(ip_numnam, NAM_CS  )!reading the block NAM_CS

    IF(ala_action.EQ.RUN_IPF)THEN
      REWIND(ip_numnam)
      READ(ip_numnam, NAM_IPF  )!reading the block NAM_IPF
    END IF
    CLOSE(ip_numnam)
    SELECT CASE (ala_reverse)
      CASE ('REVERSE', 'reverse')
        ig_reverse = M1QN3_REVERSE
      CASE ('DIRECT', 'direct')
        ig_reverse = M1QN3_DIRECT
      CASE DEFAULT
        CALL stop_program( 'In load_namelist; bad value of m1qn3 communication: '//ala_reverse)
    END SELECT

    tg_ep%aa_solver_path  = TRIM(ala_solver_path)
    tg_ep%aa_simul_action = TRIM(ala_action)
    tg_ep%prefix = toLower( TRIM(ala_action) )

    IF(ala_action.EQ.RUN_IPF)THEN
      tg_ep%i_ipf_nparticle = il_ipf_nparticle
      tg_ep%i_ipf_maxiter   = il_ipf_maxiter
      tg_ep%r_ipf_tol       = rl_ipf_tol
    END IF
    tg_ep%l_run_from_ctl          = ll_run_from_ctl
    tg_ep%l_first_simul           = .TRUE.
    tg_ep%l_restart_from_previous = ll_restart_from_previous

    !setup variables
    lm_B_matrix = ll_B_matrix
    ama_B_fName    = B_fName
    ama_Binv_Fname = Binv_Fname
    ama_L_fName    = L_fName
    ama_U_fName    = U_fName
    ama_Bctl_Fname = B_ctl_Fname
    ama_Tctl_Fname = T_ctl_Fname

    tg_ep%r_cost      = 0.0_cp
    !setting the gradient test parameters
    rm_gradtest_sFactor = starting_factor
    im_gradtest_nFactor = nb_factor

    CALL set_data_fileNames( ctl_bName, bctl_bName, obs_bName, dmt_bName, ims_bName, ogap_bName )
    CALL set_data_dir(input_dir, output_dir)
    !CALL print_eep()
    CALL load_debugFlag(fName)
  END SUBROUTINE load_namelist

  SUBROUTINE m1qn3_print_stat(id_omode, id_niter, id_nsim)
    INTEGER, INTENT(IN) :: id_omode, id_niter, id_nsim
    CALL debug('', 'm1qn3_print_stat:: output from m1qn3 $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',tag=dALLWAYS)
    CALL debug(id_niter, '  **Number of iterations = ',tag=dALLWAYS)
    CALL debug(id_nsim , '  **Number of simulations = ',tag=dALLWAYS)
    CALL debug(id_omode, '  **M1QN3 ends with output mode = ',tag=dALLWAYS)
    SELECT CASE(id_omode)
    CASE (0)
        CALL debug('', '  The simulator ask to stop by returning the value indic = 0',tag=dALLWAYS)
    CASE (1)
        CALL debug('', '  Normal way of stopping for m1qn3: the test on the gradient is satisfied',tag=dALLWAYS)
    CASE (2)
        CALL debug('', '  One of the input arguments is not well initialized',tag=dALLWAYS)
        !CALL debug('', '    - n <= 0, niter <= 0, nsim <= 0, dxmin <= 0.0 or epsg not in ]0, 1[;',tag=dALLWAYS)
        IF(ig_npctl<=0)&
            CALL debug(ig_npctl,   '     * n <= 0, current value = ',tag=dALLWAYS)
        IF(ig_niter<=0)&
            CALL debug(ig_niter,  '     * niter <= 0, current value = ',tag=dALLWAYS)
        IF(ig_nsimul<=0)&
            CALL debug(ig_nsimul, '     * nsim <= 0, current value = ',tag=dALLWAYS)
        IF(rg_xmin<=0)&
            CALL debug(rg_xmin,   '     * dxmin <= 0, current value = ',tag=dALLWAYS)
        IF(rg_epsg<=0)&
            CALL debug(rg_epsg,   '     * epsg not in ]0, 1[, current value = ',tag=dALLWAYS)
        IF( (iga_imode(1)==0).AND.(ig_ndz< 6*ig_npctl+1) )THEN!DIS mode
          CALL debug('', '    - Running DIS mode: this mode assumes ndz < 6n + 1',tag=dALLWAYS)
          CALL debug((/ig_ndz, 6*ig_npctl+1/), '     * actual values, ndz, 6n+1 = ',tag=dALLWAYS)
        ELSE IF( (iga_imode(1)==1).AND.(ig_ndz< 6*ig_npctl+1) )THEN!SIS mode
          CALL debug('', '    - Running SIS mode: this mode assumes ndz < 5n + 1',tag=dALLWAYS)
          CALL debug((/ig_ndz, 5*ig_npctl+1/), '     * actual values, ndz, 5n+1 = ',tag=dALLWAYS)
        END IF
        CALL debug('', '    - the contents of iz is not good for a warm restart,',tag=dALLWAYS)
        CALL debug('', '    - the starting point is almost optimal',tag=dALLWAYS)
    CASE (3)
        CALL debug('', '  The line search is blocked on tmax',tag=dALLWAYS)
    CASE (4)
        CALL debug('', '  The maximal number of iterations is reached',tag=dALLWAYS)
    CASE (5)
        CALL debug('', '  The maximal number of simulation is reached',tag=dALLWAYS)
    CASE (6)
        CALL debug('', '  Stop on dxmin during the line-search',tag=dALLWAYS)
    CASE (7)
        CALL debug('', '  Non negative dot product <g,d> or <y,s>',tag=dALLWAYS)
    CASE DEFAULT
        CALL debug('', '  Unknown output mode',tag=dALLWAYS)
    END SELECT
    CALL debug('', 'm1qn3_print_stat:: $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',tag=dALLWAYS)
  END SUBROUTINE m1qn3_print_stat

  SUBROUTINE print_namelist_variables()
    CHARACTER(LEN=ip_snl)   :: array_rformat, array_iformat
    array_rformat = "(A,"//TRIM( NUM2STR(ig_npctl) )//RFORMAT//")"
    array_iformat = "(A,"//TRIM( NUM2STR(ig_npctl) )//IFORMAT//")"

    CALL debug( ''             , 'Print_namelist_variables------------------', tag=dALLWAYS )
    CALL debug( ''             , ' NAM_general *****************************', tag=dALLWAYS )
    CALL debug( ''             , ' NAM_m1qn3 *******************************', tag=dALLWAYS )
    CALL debug( rg_xmin        , '   rg_xmin   = ', tag=dALLWAYS )
    CALL debug( rg_epsg        , '   rg_epsg   = ', tag=dALLWAYS )
    CALL debug( aga_norm       , '   aga_norm  = ', tag=dALLWAYS )
    CALL debug( ig_impres      , '   ig_impres = ', tag=dALLWAYS )
    CALL debug( ig_m1qn3io     , '   ig_m1qn3io= ', tag=dALLWAYS )
    CALL debug( iga_imode      , '   iga_imode = ', tag=dALLWAYS )
    CALL debug( ig_niter       , '   ig_niter  = ', tag=dALLWAYS )
    CALL debug( ig_nsimul      , '   ig_nsimul = ', tag=dALLWAYS )
    CALL debug( ig_npairs      , '   ig_npairs = ', tag=dALLWAYS )
    CALL debug( ''             , ' NAM_obs *********************************', tag=dALLWAYS )
    CALL debug( ''             , ' NAM_control *****************************', tag=dALLWAYS )
    CALL debug( ig_nctl        , '   ig_nctl       = ', tag=dALLWAYS )
    CALL debug( ''             ,' Other *******************************', tag=dALLWAYS )
    CALL debug( lg_solve_all_together       , '   lg_solve_all_together = ', tag=dALLWAYS )
    CALL debug( ''             , 'End of print_namelist_variables ##########', tag=dALLWAYS )
  END SUBROUTINE print_namelist_variables

END PROGRAM simul
