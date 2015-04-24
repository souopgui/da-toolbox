!>file solver_tools.f90
!!@brief generic module for solvers. Define the solver for external case.
!!@details the function run_solver should be rewritten for internal case
!! where the solver is simply a routine
!<
MODULE solver_tools
  USE general_constant
  USE debug_tools
  USE com_tools!, ONLY: exchange_param
IMPLICIT NONE
		PRIVATE
		PUBLIC init_solver, finalize_solver, run_solver, save_ctl_plot_data, prosca, ctonb, ctcab

CONTAINS

  !> @brief Initialize the solver
  !! @details important for internal solver where there is a need to
  !! prepare the system, allocate memory, open files, etc.
  !! this subroutine can be used to load background or control variable
  !! is they are store in a specific format
  !<
  SUBROUTINE init_solver(td_ep, id_nctl)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    INTEGER, INTENT(INOUT) :: id_nctl
    CHARACTER(LEN=ip_snl)  :: ala_inputfName

    SELECT CASE(td_ep%aa_simul_action)
      CASE (MAKE_OBS, RUN_DIRECT)!working with ctl
        IF(td_ep%l_run_from_ctl) THEN
            ala_inputfName = make_fileName(CTL_DATA, INPUT_FILE)
            CALL read_ctl_size(td_ep, ala_inputfName)
            CALL debug( get_ctlSize(td_ep), 'In init_solver, get_ctlSize(td_ep) = ')
            CALL debug( MAXVAL(td_ep%ra_ctl), 'In init_solver, MAXVAL(td_ep%ra_ctl) = ')
            CALL debug( MAXVAL(td_ep%ra_b_ctl), 'In init_solver, MAXVAL(td_ep%ra_b_ctl) = ')
            CALL read_ep_data(td_ep, CTL_DATA, INPUT_FILE )
            CALL debug( MAXVAL(td_ep%ra_ctl), 'In init_solver, MAXVAL(td_ep%ra_ctl) = ')
            CALL debug( MAXVAL(td_ep%ra_b_ctl), 'In init_solver, MAXVAL(td_ep%ra_b_ctl) = ')
            id_nctl = get_ctlSize(td_ep)
        ELSE
            CALL set_ctlsize(td_ep, id_nctl)
        END IF
      CASE (RUN_COST, RUN_ADJOINT, RUN_GRADIENT)
        CALL debug('','Running for cost, adjoint or gradient out of assimilation process recquired special precautions')
        CALL debug('','To do so, make sure that recquired temporary data are available')
        CALL dpause('Consider using CTRL-C to abort or ENTER to try')
      CASE (RUN_ASSIM, RUN_IPF, GRAD_TEST)!working with background ctl
        ala_inputfName = make_fileName(BCTL_DATA, INPUT_FILE)
        CALL read_ctl_size(td_ep, ala_inputfName)
        CALL read_ep_data(td_ep, BCTL_DATA, INPUT_FILE )
        id_nctl = get_ctlSize(td_ep)
        td_ep%l_run_from_ctl = .TRUE.
      CASE (RUN_ASSIM2)!working with background ctl
        ala_inputfName = make_fileName(ACTL_DATA, OUTPUT_FILE)
        CALL read_ctl_size(td_ep, ala_inputfName)
        CALL read_ep_data(td_ep, ACTL_DATA, OUTPUT_FILE )
        id_nctl = get_ctlSize(td_ep)
        td_ep%l_run_from_ctl = .TRUE.

        td_ep%ra_sigma_ctl = REAL(10**10, cp)!very large value
        td_ep%ra_b_ctl = td_ep%ra_ctl
        td_ep%ra_ctl   = 0.0_cp
        WHERE( td_ep%ra_sigma_ctl > 0.0_cp ) td_ep%ra_Bm1 = 1.0_cp/td_ep%ra_sigma_ctl**2
        td_ep%aa_simul_action = RUN_ASSIM

      CASE (MAKE_CTL, MAKE_BG, OBS_TO_BG)!initializing from namelist
        CALL set_ctlsize(td_ep, id_nctl)
      CASE DEFAULT!unknown
        CALL stop_program(td_ep%aa_simul_action, 'In init_solver; bad value of simul action: ')
    END SELECT
  END SUBROUTINE init_solver

  !> @brief Finalize the solver
  !! @details important for internal solver where there is a need
  !! to clean the system, free memory, close files, etc.
  !! this subroutine can be used to save the analysed control variable
  !! if a specific format is used, in which case the change of variable
  !! is applied, and equilibration.
  !<
  SUBROUTINE finalize_solver(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep

    CALL nothing(td_ep%r_nothing)
  END SUBROUTINE finalize_solver

  !> @brief Run the solver to respond to the request specified by
  !> the exchange parameter
  !! @details this routine correspond to the external solver to
  !! be run as independent program by system call
  !<
  SUBROUTINE run_solver( td_ep )
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local var
    CHARACTER(len=*), PARAMETER :: apa_rs = 'RUN_SOLVER'


		SELECT CASE(td_ep%aa_solver_action)
			CASE( RUN_COSTGRAD )
				CALL delete_rtFile( COST_DATA )
				CALL delete_rtFile( GRAD_DATA )
			CASE( RUN_COST )
				CALL delete_rtFile(COST_DATA )
			CASE( RUN_GRADIENT )
				CALL delete_rtFile( GRAD_DATA )
			CASE DEFAULT !This value of indic is not supported
			!Nothing to do
		END SELECT

    CALL write_ep_data(td_ep, EP_DATA, RTIME_FILE, apa_rs )
    !CALL debug('', 'run_solver: calling the solver with action = <'//TRIM(td_ep%aa_solver_action)//'>', tag=dALLWAYS)
    CALL system_run( TRIM(td_ep%aa_solver_path) )!running the external solver
		SELECT CASE(td_ep%aa_solver_action)
			CASE(RUN_COSTGRAD)
				CALL read_ep_data(td_ep, COST_DATA, RTIME_FILE )
				CALL read_ep_data(td_ep, GRAD_DATA, RTIME_FILE )
			CASE(RUN_COST)
				CALL read_ep_data(td_ep, COST_DATA, RTIME_FILE )
			CASE(RUN_GRADIENT)
				CALL read_ep_data(td_ep, GRAD_DATA, RTIME_FILE )
				!CALL debug(td_ep%ra_grad, ' td_ep%ra_grad = ')
				!CALL dpause()
			CASE DEFAULT !This value of indic is not supported
			!Nothing to do
		END SELECT
    !CALL debug('', 'In run_solver, resuming with action = <'//TRIM(td_ep%aa_solver_action)//'>', tag=dALLWAYS)
  END SUBROUTINE run_solver

  !> @brief change of variable : transforms the CTL of the preconditioned pb
  !> to the CTL of the original pb
  !! @param [in,out] td_ep exchange parameter
  !<
  SUBROUTINE chavar(td_ep)
    TYPE(exchange_param), INTENT(IN OUT) :: td_ep

    CALL apply_square_cov( td_ep%ra_c_ctl, td_ep%ra_dctl, die=.FALSE. )
  END SUBROUTINE chavar

	!> @brief interface for saving ctl in simple plot format
	!! @param [in] td_ep exchange parameter
	!! @param [in] id_dType data type (CTL_DATA, ...)
	!!
	!! this routine should be adapted to the project
	!<
	SUBROUTINE save_ctl_plot_data(td_ep, id_dType)
		TYPE(exchange_param), INTENT(IN) :: td_ep
		INTEGER, INTENT(IN) :: id_dType
		!local variables
		TYPE(dyn_rVector), DIMENSION(2), SAVE :: tla_save
		CHARACTER(len=ip_snl):: ala_plotFName

		!saving ctl for plotting
		!CALL debug(100, 'In save_ctl_plot_data ',tag=dALLWAYS)
		ala_plotFName = make_pFileName( id_dType, OUTPUT_FILE, prefix=td_ep%prefix )
		SELECT CASE(id_dType)
			CASE(CTL_DATA, ACTL_DATA)
				CALL write_vector_for_plot(td_ep%ra_b_ctl + td_ep%ra_ctl, ala_plotFName)
				!CALL debug(351, 'In save_ctl_plot_data ',tag=dALLWAYS)
			CASE(BCTL_DATA)
				CALL write_vector_for_plot(td_ep%ra_b_ctl, ala_plotFName)
				!CALL debug(352, 'In save_ctl_plot_data ',tag=dALLWAYS)
		END SELECT
		!CALL debug(400, 'In save_ctl_plot_data ',tag=dALLWAYS)
		!end of saving obs for plotting
	END SUBROUTINE save_ctl_plot_data

  !> @brief compute the inner product of \a rda_u et \a rda_v
  !<
  SUBROUTINE prosca(id_ctl, rda_u, rda_v, rd_ps, ida_wa, rda_wa, dda_wa)
    INTEGER , INTENT(IN) :: id_ctl
    REAL(cp), DIMENSION(id_ctl), INTENT(IN) :: rda_u, rda_v
    REAL(cp), INTENT(OUT) :: rd_ps
    INTEGER,DIMENSION(2),INTENT(IN) :: ida_wa
    REAL(KIND=sp), DIMENSION(2),INTENT(IN) :: rda_wa
    REAL(KIND=dp), DIMENSION(2),INTENT(IN) :: dda_wa

    rd_ps = DOT_PRODUCT(rda_u, rda_v)
    CALL nothing(ida_wa, rda_wa, dda_wa)
  END SUBROUTINE prosca

  !> @brief Change of basis
  SUBROUTINE ctonb(id_ctl, rda_x, rda_y, ida_wa, rda_wa, dda_wa)
    INTEGER , INTENT(IN) :: id_ctl
    REAL(cp), DIMENSION(id_ctl), INTENT(IN) :: rda_x
    REAL(cp), DIMENSION(id_ctl), INTENT(OUT) :: rda_y
    INTEGER,DIMENSION(2),INTENT(IN) :: ida_wa
    REAL(KIND=sp), DIMENSION(2),INTENT(IN) :: rda_wa
    REAL(KIND=dp), DIMENSION(2),INTENT(IN) :: dda_wa

    rda_y = rda_x
    CALL nothing(ida_wa, rda_wa, dda_wa)
  END SUBROUTINE ctonb

  !> @brief Change of basis, reverse operation associated with ctonb
  SUBROUTINE ctcab(id_ctl, rda_y, rda_x, ida_wa, rda_wa, dda_wa)
    INTEGER , INTENT(IN) :: id_ctl
    REAL(cp), DIMENSION(id_ctl), INTENT(IN) :: rda_y
    REAL(cp), DIMENSION(id_ctl), INTENT(OUT) :: rda_x
    INTEGER,DIMENSION(2),INTENT(IN) :: ida_wa
    REAL(KIND=sp), DIMENSION(2),INTENT(IN) :: rda_wa
    REAL(KIND=dp), DIMENSION(2),INTENT(IN) :: dda_wa

    rda_x = rda_y
    CALL nothing(ida_wa, rda_wa, dda_wa)
  END SUBROUTINE ctcab
END MODULE solver_tools