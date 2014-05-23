!>file simul_tools.f90
!!\brief global communication tools
!! The subroutine simulator had been put in its own file for conformance with gfortran. It doesn compile when th call to m1qn3 and the simulator routine are in the same file
!<
MODULE simul_tools
  USE com_tools
  USE solver_tools
IMPLICIT NONE

  INTEGER, PARAMETER ::  &
       M1QN3_DIRECT   = 0,& !> direct communication in M1QN3, used for parameter reverse
       M1QN3_REVERSE  = 1,& !> reverse communication in M1QN3, used for parameter reverse
       SIMUL_NOTHING  = 1,& !> nothing required by M1QN3, used for parameter indic
       SIMUL_COSTGRAD = 4,& !> M1QN3 requires cost function and its gradient, used for parameter indic
       SIMUL_COSTONLY = 11,& !> only the cost function is required, used for parameter indic
       SIMUL_ASK_STOP = 0 ,& !> something that M1QN3 can't not understand happened
       SIMUL_IMPOSSIBLE=-1  !> impossible to compute the cost function and/or its gradient


  TYPE(exchange_param), SAVE :: tg_ep
  INTEGER, PRIVATE, SAVE :: im_iterCount, im_simulCount

  !> \brief Shall direct and adjoint problems be solved simultaneously?
  !! this is used for solvers that can solve direct and adjoint problems at the same time, especially space time solvers
  !<
  LOGICAL  :: lg_solve_all_together

CONTAINS

  !> \brief Computes the cost function and its gradient for conmin library
  !! \param [in] id_ctl size of the control vector
  !! \param [in,out] rda_ctl control vector
  !! \param [in,out] rd_cost value of the cost function,
  !! \param [in,out] rda_grad Gradient of the cost function
  SUBROUTINE cal_fg(id_ctl, rda_ctl, rd_cost, rda_grad)
    IMPLICIT NONE
    !parameters for m1qn3 simulator
    INTEGER, INTENT(IN)    :: id_ctl
    REAL(KIND=cp),DIMENSION(id_ctl),INTENT(INOUT) :: rda_ctl
    REAL(KIND=cp),DIMENSION(id_ctl),INTENT(INOUT) :: rda_grad
    REAL(KIND=cp),INTENT(INOUT) :: rd_cost
    !*V [idr]da_wa : working array inustite (contrainte modulopt)
    INTEGER :: il_indic
    INTEGER,DIMENSION(2) :: ila_wa
    REAL(KIND=sp), DIMENSION(2) :: rla_wa
    REAL(KIND=dp), DIMENSION(2) :: dla_wa

    il_indic = SIMUL_COSTGRAD
    CALL simulator(il_indic, id_ctl, rda_ctl, rd_cost, rda_grad, ila_wa, rla_wa, dla_wa)
  END SUBROUTINE cal_fg

  !> \brief Computes the cost function and its gradient
  !! \param [in,out] id_indic Indicates what is required from the simulator
  !! \param [in] id_ctl size of the control vector
  !! \param [in,out] rda_ctl control vector
  !! \param [in,out] rda_grad Gradient of the cost function
  !! \param [in,out] rd_cost value of the cost function,
  !! \param [in] ida_wa Integer working array for simulator
  !! \param [in] rda_wa Simple precision working array for simulator
  !! \param [in] dda_wa Double precision working array for simulator
  !! \detail parameter id_indic has the following meaningful values
  !! - 4 : compute the cost function and its gradient
  !! - 1 : informes the simulator that the minimizer is moving to the next iteration
  !<
  SUBROUTINE simulator(id_indic, id_ctl, rda_ctl, rd_cost, rda_grad, ida_wa, rda_wa, dda_wa)
    IMPLICIT NONE
    !parameters for m1qn3 simulator
    INTEGER, INTENT(INOUT) :: id_indic
    INTEGER, INTENT(IN)    :: id_ctl
    REAL(KIND=cp),DIMENSION(id_ctl),INTENT(INOUT) :: rda_ctl
    REAL(KIND=cp),DIMENSION(id_ctl),INTENT(INOUT) :: rda_grad
    REAL(KIND=cp),INTENT(INOUT) :: rd_cost
    !*V [idr]da_wa : working array inustite (contrainte modulopt)
    INTEGER,DIMENSION(2),INTENT(IN) :: ida_wa
    REAL(KIND=sp), DIMENSION(2),INTENT(IN) :: rda_wa
    REAL(KIND=dp), DIMENSION(2),INTENT(IN) :: dda_wa
    !local variables
    !REAL(KIND=cp) :: rl_costb, rl_tmp
    CHARACTER(LEN=80)   :: array_format
    LOGICAL :: ll_diverged
    INTEGER :: il_NaNCount, il_InfCount

    array_format = "("//TRIM( NUM2STR(id_ctl+1) )//RFORMAT//")"
    !
    CALL nothing(ida_wa, rda_wa, dda_wa)!to avoid warning of unused variables
    !
    !PRINT*, 'Entering Simulator ---------------------------------'
    SELECT CASE (id_indic)
    CASE(SIMUL_COSTGRAD) !required : cost function and its gradient
      !CALL debug('', ' In simul_tools::simulator; calling routine asks for cost and grad ', tag=dALLWAYS)
       !** initializing exchange parameters
       !tg_ep%i_nctl = id_ctl
       tg_ep%ra_pctl = rda_ctl

       IF(get_nb_outOfBound(tg_ep)/=0)THEN
          id_indic = SIMUL_IMPOSSIBLE
          CALL print_ctl_bound(tg_ep)
       ELSE
          !/!\ lg_solve_all_together should be set to .TRUE. only if the solver can solve the direct and the adjoint problem simultaneously
          IF(lg_solve_all_together)THEN
             tg_ep%aa_solver_action=RUN_COSTGRAD
             CALL run_solver(tg_ep)
             il_NaNCount = NaNCount(tg_ep%ra_grad)
             il_InfCount = InfCount(tg_ep%ra_grad)
             ll_diverged = isNaN(tg_ep%r_cost).OR.isInf(tg_ep%r_cost)&
              .OR.((il_NaNCount+il_InfCount).GT.0).OR.tg_ep%l_simul_diverged
          ELSE
            !**   1. Computing the cost function
            !PRINT*, 'integrating direct '

            tg_ep%aa_solver_action=RUN_COST
            CALL run_solver(tg_ep)
            ll_diverged = isNaN(tg_ep%r_cost).OR.isInf(tg_ep%r_cost).OR.tg_ep%l_simul_diverged
            IF( .NOT.ll_diverged ) THEN
              !**   2. computing the gradient
              tg_ep%aa_solver_action=RUN_GRADIENT
              CALL run_solver(tg_ep)
              il_NaNCount = NaNCount(tg_ep%ra_grad)
              il_InfCount = InfCount(tg_ep%ra_grad)
              ll_diverged = ((il_NaNCount+il_InfCount).GT.0).OR.tg_ep%l_simul_diverged
            END IF
          END IF!lg_solve_all_together
          IF( .NOT. ll_diverged ) THEN
            rd_cost  = tg_ep%r_cost
            rda_grad = tg_ep%ra_grad
          ELSE
            id_indic = SIMUL_IMPOSSIBLE
          END IF
       END IF!(il_nb_outOfBound/=0)
       tg_ep%l_first_simul = .FALSE.
       im_simulCount = im_simulCount+1
       WRITE(ip_cost_evol, "(A,I4,A,I4,A,ES18.9E2)")"Iter ", im_iterCount, " Simul ", im_simulCount, ": cost = ", tg_ep%r_cost
       tg_ep%i_nsimul_in_iter = tg_ep%i_nsimul_in_iter + 1
    CASE(SIMUL_COSTONLY)!only the cost function is required, used for the test of the gradient
      !CALL debug('', ' In simul_tools::simulator; calling routine asks for cost only ', tag=dALLWAYS)
      tg_ep%ra_pctl = rda_ctl

      IF(get_nb_outOfBound(tg_ep)/=0)THEN
        id_indic = SIMUL_IMPOSSIBLE
        CALL print_ctl_bound(tg_ep)
        ll_diverged = .TRUE.
      ELSE
        tg_ep%aa_solver_action=RUN_COST
        CALL run_solver(tg_ep)
        il_NaNCount = NaNCount(tg_ep%ra_grad)
        il_InfCount = InfCount(tg_ep%ra_grad)
        ll_diverged = isNaN(tg_ep%r_cost).OR.isInf(tg_ep%r_cost)&
        .OR.((il_NaNCount+il_InfCount).GT.0).OR.tg_ep%l_simul_diverged
      END IF
      IF( .NOT. ll_diverged ) THEN
        rd_cost  = tg_ep%r_cost
      ELSE
        id_indic = SIMUL_IMPOSSIBLE
      END IF
    CASE(SIMUL_NOTHING)! m1qn3 is moving to the next iteration
      CALL debug('', ' In simul_tools::simulator; calling routine asks for nothing ', tag=dALLWAYS)
      tg_ep%aa_solver_action = RUN_NOTHING
      CALL run_solver(tg_ep)
      WRITE(ip_ctl_out, FMT=array_format)tg_ep%ra_dctl+tg_ep%ra_b_ctl, rd_cost
      WRITE(ip_grad_out, FMT=array_format)rda_grad
      CALL debug('', ' In simul_tools::simulator; m1qn3 moving to the next iteration ', tag=dALLWAYS)
      im_iterCount = im_iterCount+1
      tg_ep%i_iter = im_iterCount
      tg_ep%i_nsimul_in_iter = 0
      WRITE(ip_cost_evol, *)"Moving to iteration ", im_iterCount
    CASE DEFAULT !This value of indic is not supported
      CALL stop_program(id_indic, ' In simul_tools::simulator; indic value not supported: ')
    END SELECT!( id_indic == M1QN3_COSTGRAD)
    CALL debug('Exiting Simulator -----------------------------')
  END SUBROUTINE simulator

  !> \brief Initializes the counters variables used by the simulator
  !!
  !<
  SUBROUTINE init_simul_counter()
    im_iterCount = 0
    im_simulCount= -1 !to account for the initialization
  END SUBROUTINE init_simul_counter

END MODULE simul_tools