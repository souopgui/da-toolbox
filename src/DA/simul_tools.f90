!> @file simul_tools.f90
!! @brief global communication tools
!! The subroutine simulator had been put in its own file for conformance with gfortran. It doesn compile when th call to m1qn3 and the simulator routine are in the same file
!<
module simul_tools
    use com_tools
    use solver_tools
implicit none

    integer, parameter ::    &
        M1QN3_DIRECT   = 0 ,& !> direct communication in M1QN3, used for parameter reverse
        M1QN3_REVERSE  = 1 ,& !> reverse communication in M1QN3, used for parameter reverse
        SIMUL_NOTHING  = 1 ,& !> nothing required by M1QN3, used for parameter indic
        SIMUL_COSTGRAD = 4 ,& !> M1QN3 requires cost function and its gradient, used for parameter indic
        SIMUL_COSTONLY = 11,& !> only the cost function is required, used for parameter indic
        SIMUL_ASK_STOP = 0 ,& !> something that M1QN3 can't not understand happened
        SIMUL_IMPOSSIBLE=-1  !> impossible to compute the cost function and/or its gradient


    type(exchange_param), save :: tg_ep
    integer, private, save :: im_iterCount, im_simulCount

    !> \brief Shall direct and adjoint problems be solved simultaneously?
    !! this is used for solvers that can solve direct and adjoint problems at the same time, especially space time solvers
    !<
    logical  :: lg_solve_all_together

contains

    !> \brief Computes the cost function and its gradient for conmin library
    !! \param [in] id_ctl size of the control vector
    !! \param [in,out] rda_ctl control vector
    !! \param [in,out] rd_cost value of the cost function,
    !! \param [in,out] rda_grad Gradient of the cost function
    subroutine cal_fg(id_ctl, rda_ctl, rd_cost, rda_grad)
        implicit none
        !parameters for m1qn3 simulator
        integer, intent(in)    :: id_ctl
        real(kind=cp),dimension(id_ctl),intent(in) :: rda_ctl
        real(kind=cp),dimension(id_ctl),intent(inout) :: rda_grad
        real(kind=cp),intent(inout) :: rd_cost
        !*v [idr]da_wa : working array inustite (contrainte modulopt)
        integer :: il_indic
        integer,dimension(2) :: ila_wa
        real(kind=sp), dimension(2) :: rla_wa
        real(kind=dp), dimension(2) :: dla_wa

        il_indic = SIMUL_COSTGRAD
        call simulator(il_indic, id_ctl, rda_ctl, rd_cost, rda_grad, ila_wa, rla_wa, dla_wa)
    end subroutine cal_fg

    !> @brief Computes the cost function and its gradient
    !! @param [in,out] id_indic Indicates what is required from the simulator
    !! @param [in] id_ctl size of the control vector
    !! @param [in,out] rda_ctl control vector
    !! @param [in,out] rda_grad Gradient of the cost function
    !! @param [in,out] rd_cost value of the cost function,
    !! @param [in] ida_wa Integer working array for simulator
    !! @param [in] rda_wa Simple precision working array for simulator
    !! @param [in] dda_wa Double precision working array for simulator
    !! @details parameter id_indic has the following meaningful values
    !! - 4 : compute the cost function and its gradient
    !! - 1 : informes the simulator that the minimizer is moving to the next iteration
    !<
    subroutine simulator(id_indic, id_ctl, rda_ctl, rd_cost, rda_grad, ida_wa, rda_wa, dda_wa)
        implicit none
        !parameters for m1qn3 simulator
        integer, intent(inout) :: id_indic
        integer, intent(in)    :: id_ctl
        real(kind=cp),dimension(id_ctl),intent(in) :: rda_ctl
        real(kind=cp),dimension(id_ctl),intent(inout) :: rda_grad
        real(kind=cp),intent(inout) :: rd_cost
        !*v [idr]da_wa : working array inustite (contrainte modulopt)
        integer,dimension(2),intent(in) :: ida_wa
        real(kind=sp), dimension(2),intent(in) :: rda_wa
        real(kind=dp), dimension(2),intent(in) :: dda_wa
        !local variables
        !real(kind=cp) :: rl_costb, rl_tmp
        character(len=80)   :: array_format
        logical :: ll_diverged
        integer :: il_NaNCount, il_InfCount

        array_format = "("//trim( NUM2STR(id_ctl+1) )//RFORMAT//")"
        !
        call nothing( ida_wa, rda_wa, dda_wa )!to avoid warning of unused variables
        !
        !PRINT*, 'Entering Simulator ---------------------------------'
        select case ( id_indic )
            case( SIMUL_COSTGRAD ) !required : cost function and its gradient
                !CALL debug('', ' In simul_tools::simulator; calling routine asks for cost and grad ', tag=dALLWAYS)
                !** initializing exchange parameters
                !tg_ep%i_nctl = id_ctl
                tg_ep%ra_pctl = rda_ctl

                if(get_nb_outOfBound(tg_ep)/=0)then
                    id_indic = SIMUL_IMPOSSIBLE
                    call print_ctl_bound(tg_ep)
                else
                    !/!\ lg_solve_all_together should be set to .TRUE. only if
                    !    the solver can solve the direct and the adjoint problem
                    !    simultaneously
                    if(lg_solve_all_together)then
                        tg_ep%aa_solver_action=RUN_COSTGRAD
                        call run_solver(tg_ep)
                        il_NaNCount = NaNCount(tg_ep%ra_grad)
                        il_InfCount = InfCount(tg_ep%ra_grad)
                        ll_diverged = isNaN(tg_ep%r_cost).or.isInf(tg_ep%r_cost)&
                        .or.((il_NaNCount+il_InfCount).gt.0).or.tg_ep%l_simul_diverged
                    else
                        !**   1. Computing the cost function
                        !PRINT*, 'integrating direct '

                        tg_ep%aa_solver_action=RUN_COST
                        call run_solver(tg_ep)
                        ll_diverged = isNaN(tg_ep%r_cost).or.isInf(tg_ep%r_cost).or.tg_ep%l_simul_diverged
                        if( .not.ll_diverged ) then
                        !**   2. computing the gradient
                        tg_ep%aa_solver_action=RUN_GRADIENT
                        call run_solver(tg_ep)
                        il_NaNCount = NaNCount(tg_ep%ra_grad)
                        il_InfCount = InfCount(tg_ep%ra_grad)
                        ll_diverged = ((il_NaNCount+il_InfCount).gt.0).or.tg_ep%l_simul_diverged
                        end if
                    end if!lg_solve_all_together
                    if( .not. ll_diverged ) then
                        rd_cost  = tg_ep%r_cost
                        rda_grad = tg_ep%ra_grad
                    else
                        id_indic = SIMUL_IMPOSSIBLE
                    end if
                end if!(il_nb_outofbound/=0)
                tg_ep%l_first_simul = .FALSE.
                im_simulCount = im_simulCount+1
                write(ip_cost_evol, "(A,I4,A,I4,A,ES18.9E2)")"Iter ", im_iterCount&
                    , " Simul ", im_simulCount, ": cost = ", tg_ep%r_cost
                tg_ep%i_nsimul_in_iter = tg_ep%i_nsimul_in_iter + 1
            case(SIMUL_COSTONLY)!only the cost function is required, used for the test of the gradient
                !CALL debug('', ' In simul_tools::simulator; calling routine asks for cost only ', tag=dALLWAYS)
                tg_ep%ra_pctl = rda_ctl

                if(get_nb_outOfBound(tg_ep)/=0)then
                    id_indic = SIMUL_IMPOSSIBLE
                    call print_ctl_bound(tg_ep)
                    ll_diverged = .true.
                else
                    tg_ep%aa_solver_action=RUN_COST
                    call run_solver(tg_ep)
                    il_NaNCount = NaNCount(tg_ep%ra_grad)
                    il_InfCount = InfCount(tg_ep%ra_grad)
                    ll_diverged = isNaN(tg_ep%r_cost).or.isInf(tg_ep%r_cost)&
                    .or.((il_NaNCount+il_InfCount).gt.0).or.tg_ep%l_simul_diverged
                end if
                if( .not. ll_diverged ) then
                    rd_cost  = tg_ep%r_cost
                else
                    id_indic = SIMUL_IMPOSSIBLE
                end if
            case(SIMUL_NOTHING)! m1qn3 is moving to the next iteration
                call debug('', ' In simul_tools::simulator; calling routine asks for nothing ', tag=dALLWAYS)
                tg_ep%aa_solver_action = RUN_NOTHING
                call run_solver(tg_ep)
                write(ip_ctl_out, fmt=array_format)tg_ep%ra_dctl+tg_ep%ra_b_ctl, rd_cost
                write(ip_grad_out, fmt=array_format)rda_grad
                call debug('', ' In simul_tools::simulator; m1qn3 moving to the next iteration ', tag=dALLWAYS)
                im_iterCount = im_iterCount+1
                tg_ep%i_iter = im_iterCount
                tg_ep%i_nsimul_in_iter = 0
                write(ip_cost_evol, *)"Moving to iteration ", im_iterCount
            case default !This value of indic is not supported
                call stop_program(id_indic, ' In simul_tools::simulator; indic value not supported: ')
        end select!( id_indic == M1QN3_COSTGRAD)
        call debug('Exiting Simulator -----------------------------')
    end subroutine simulator


    !> @brief Computes the gradient of the cost function at a given CTL
    !! @param [in] rda_ctl control vector
    !! @param [in,out] rda_grad Gradient of the cost function
    !! @param [out] diverged says if the simulation diverges
    !<
    subroutine simulator_grad(rda_ctl, rda_grad, diverged)
        implicit none
        real(kind=cp),dimension(:),intent(in) :: rda_ctl
        real(kind=cp),dimension(:),intent(inout) :: rda_grad
        logical, intent(out) :: diverged
        !local variables
        ! other parameters of the simulator
        integer,dimension(2) :: ila_wa
        real(kind=sp), dimension(2) :: rla_wa
        real(kind=dp), dimension(2) :: dla_wa
        real(dp) :: cost
        integer :: indic, nctl

        indic = SIMUL_COSTGRAD
        nctl = size(rda_ctl)

        call simulator(indic, nctl, rda_ctl, cost, rda_grad, ila_wa, rla_wa, dla_wa)
        diverged = (indic==SIMUL_IMPOSSIBLE)
    end subroutine simulator_grad

    !> \brief Initializes the counters variables used by the simulator
    !!
    !<
    subroutine init_simul_counter()
        im_iterCount = 0
        im_simulCount= -1 !to account for the initialization
    end subroutine init_simul_counter

end module simul_tools