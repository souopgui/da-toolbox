!> \brief Module user_case, allows user to describe its case by writing th body of subroutines used by the solver.
!! \detail User should provide body of routines defined here. These routines are named user_xxx and are call by the solver.
!! many variables imported from included modules are availaible as global variables.
!! - ne number of equation defining the problems, this is usually equql to the numbers of unknowns variable to solve
!! - ng number of discrete values defining each single unknown variable
!! - dim dimension of the physical space in which the problem is being solved
!! - t current integration time
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! - 
!! 
!<
MODULE user_case
  !
  ! Case elliptic poisson 
  USE precision
  USE wlt_vars
  USE wavelet_filters_mod
  USE elliptic_mod
  USE elliptic_vars
  USE wlt_trns_mod									
  USE wlt_trns_vars
  USE wlt_trns_util_mod
  USE io_3d_vars
  USE util_mod
  USE util_vars
  USE share_consts
  USE pde
  USE sizes
  USE share_kry
  USE vector_util_mod
  USE field
  USE input_file_reader
  USE com_tools
  USE gd_tools
  USE gd_toolsadj
  USE regul_tools
  USE regul_toolsadj
  

  !constants
  !
  ! case specific variables
  !
  INTEGER n_var_vorticity ! start of vorticity in u array
  INTEGER n_var_pressure  ! start of pressure in u array

!!!!!!!!!!!!!!!!!!!
  CHARACTER(LEN=ip_snl) :: aga_cand_filename, aga_restart_filename, aga_direct_fileName,&
       aga_restart_com_fName, aga_direct_com_fName, aga_cand_com_fName
  REAL(pr) :: rg_costb
  INTEGER, PARAMETER, PRIVATE :: ip_gd_niter = 50
  REAL(pr), PRIVATE :: rg_sigma, rg_v0, rg_omega, rg_alpha
  INTEGER :: j_mx_adjoint,&!>maximul level for adjoint model
             j_ctl       ,&!>level for the control vector
             ig_av_config  !>configuration of the advection velocity
  !INTEGER , PRIVATE, DIMENSION(3) :: iga_Mvector !> Wavelets count at coarse scale
  LOGICAL :: lg_direct ,& !> running the direct model
             lg_adjoint,& !> running the direct model
             lg_restart,& !> says if initial value is read from file 
             lg_adapt_on_adjoint,& !> says if the grid should be adapted based on the adjoint variables when solving adjoint problem
             lg_adapt_on_direct    !> says if the grid should be adapted based on the direct variables when solving adjoint problem

  REAL(pr), PRIVATE, DIMENSION(3) :: rga_x0
  !the variable rga_ctl is introduced for the purpose of the change of variable
  REAL(pr), PRIVATE, ALLOCATABLE, DIMENSION(:) :: rga_ctl, rga_ctlb, rga_obsgap
  INTEGER , PRIVATE, ALLOCATABLE, DIMENSION(:) :: iga_obs_idx !> observation indices
  INTEGER , PRIVATE, ALLOCATABLE, DIMENSION(:) :: iga_ctl_idx !> indices of the control vector (initial condition), these are indices in the vectorized system state (trajectory)
  INTEGER , PRIVATE, ALLOCATABLE, DIMENSION(:, :) :: iga_ctl_coord !> Coordinates of the control vector (initial condition); these are coordinated in the grid (non adaptive)
  
  TYPE(exchange_param) :: tg_ep !> exchange parameter structure, used to exchange informations between the minimizer and the solver
  TYPE(obs_structure)  :: tg_obs

  !
  ! local variables
  !
  REAL(pr)  :: nu
  

  LOGICAL :: divgrad

  ! debug flag to be read from .inp file to test database independent interpolation
  LOGICAL, PRIVATE :: test_interpolate

CONTAINS

  !
  ! In user_setup_pde() we setup how many variables are integrated, which are interpolated
  ! to the next times step and if any exeact solution exists to compare to during the run.
  ! We also set the variable names that are used when the result files are written.
  !
  ! The following variables must be setup in this routine:
  !
  !
  ! n_integrated     ! first n_integrated eqns will be acted on for time integration
  ! n_var_additional ! 
  ! n_var
  !
  !
  !
  !
  SUBROUTINE  user_setup_pde ( VERB ) 
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: VERB         ! print debug info
    LOGICAL :: do_verb
    INTEGER :: i, ibi
    CHARACTER(LEN=50) :: ala_tmp
    
    do_verb = .TRUE.
    IF (PRESENT(VERB)) do_verb = VERB
    
    IF (do_verb) THEN
       PRINT * ,''
       PRINT *, '**********************Setting up PDE*****************'
       PRINT * ,'CASE Elliptic '
       PRINT *, '*****************************************************'
    END IF

    n_integrated = 2 ! dim

    n_var_additional = dim-1 !advection velocity, physical space only

    n_var = n_integrated + n_var_additional !--Total number of variables

    n_var_exact = 1 ! 0 <--> No exact solution

    n_var_pressure  = n_integrated + 1 ! ! no pressure

    !
    !--Allocate logical arrays 
    !  Must be done after setting number of different types of variables
    !  and before setting logical variable mapping arrays
    !
    CALL alloc_variable_mappings


    !
    ! Fill in variable names for default variables for integration meth 2
    ! These may be changed in the call to the user specified routine
    !  user_setup_pde()
    !
    ! In integrated variables (3D)
    WRITE (u_variable_names(1), u_variable_names_fmt)  'u'
    WRITE (u_variable_names(2), u_variable_names_fmt)  'p'
    !Advection velocity
    DO ibi = 1, dim-1
       ala_tmp = ""
       WRITE( ala_tmp, FMT="(A,I1)") 'w', ibi
       WRITE (u_variable_names(n_integrated + ibi), u_variable_names_fmt) TRIM(ala_tmp)
    END DO
    !
    ! setup logical variable mapping arrays
    ! This defined which variables are used for adaptation, Which need to be interpololated to the new grid at 
    ! each time step, which are saved, and for which we have an exect solutiuon we want to check at each time step
    !


    !
    ! Setup which components we will base grid adaptation on.    !
    !** Adaption strategy
    !** Always adapt on direct variables when solving the direct problem
    !** Adjoint variable, adapt on adjoint variable only when solving the adjoint problem
    !**** CASE 1 : solving direct and adjoint problems separately
    !****** CASE 1.1 adjoint run on new mesh, adapt only on adjoint variable
    !******** lg_adapt_on_adjoint = T
    !******** lg_adapt_on_direct  = F
    !****** CASE 1.2 adjoint run on old mesh, adapt only on direct variable
    !******** lg_adapt_on_adjoint = F
    !******** lg_adapt_on_direct  = T, direct variable loaded from file
    !****** CASE 1.3 adjoint run on mixte mesh, adapt on both adjoint and direct variables
    !******** lg_adapt_on_adjoint = T
    !******** lg_adapt_on_direct  = T, direct variable loaded from file
    !**** CASE 2 : solving direct and adjoint problems simultaneously, always adapt on direct variables
    !****** CASE 2.1 Adaption on adjoint variable
    !******** lg_adapt_on_adjoint = T
    !****** CASE 2.2 No adaption on adjoint variable
    !******** lg_adapt_on_adjoint = F

    !We adapt on adjoint var when solving the adjoint problem and when explicitely asked to do so.
    !We adapt on direct  var when solving the direct  problem or  when explicitely asked to do so or if no adaption is explicitely set
    !Initial adaption at first time level
    n_var_adapt(1,0) = (lg_direct .OR. lg_adapt_on_direct )
    n_var_adapt(2,0) = (lg_adjoint.AND.lg_adapt_on_adjoint)
    !--After first time step adapt on
    n_var_adapt(1,1) = (lg_direct .OR. lg_adapt_on_direct ).OR.(.NOT.lg_adapt_on_adjoint)
    n_var_adapt(2,1) = (lg_adjoint.AND.lg_adapt_on_adjoint)
    !Variables that need to be interpoleted to new adapted grid in initial grid adaptation
    n_var_interpolate(1:n_var,0) = .TRUE. 

    !Variables that need to be interpoleted to new adapted grid at each time step
    n_var_interpolate(1:n_var,1) = .TRUE. 

    !
    ! setup which components we have an exact solution for
    !
    n_var_exact_soln(:,0:1) = .TRUE.



    !
    ! setup which variables we will save the solution
    !
    n_var_save(:) = .TRUE. ! save all for restarting code


    ! Setup which variables are required on restart
    !
    n_var_req_restart = .FALSE.
    n_var_req_restart(1:n_integrated) = .TRUE.


    !
    ! Set the maximum number of components for which we have an exact solution
    !
    n_var_exact = 0 ! MAX( COUNT(n_var_exact_soln(:,0)),COUNT(n_var_exact_soln(:,1)) )

    !
    ! Setup a scaleCoeff array of we want to tweak the scaling of a variable
    ! ( if scaleCoeff < 1.0 then more points will be added to grid 
    !
    ALLOCATE ( scaleCoeff(1:n_var), Umn(1:n_var) )
    scaleCoeff = 1.0_pr

    IF (do_verb) THEN
       PRINT *, 'n_integrated = ',n_integrated 
       PRINT *, 'n_var = ',n_var 
       PRINT *, 'n_var_exact = ',n_var_exact 
       PRINT *, '*******************Variable Names*******************'
       DO i = 1,n_var
          WRITE (*, u_variable_names_fmt) u_variable_names(i)
       END DO
       PRINT *, '****************************************************'
    END IF

  END SUBROUTINE  user_setup_pde

  !> \brief read input variables from input file ("case_name"_pde.inp) 
  !!
  !! \detail
  !! Read input from "case_name"_pde.inp file
  !! case_name is in string file_name read from command line
  !! in read_command_line_input()
  !! \note Make sure that the variable file_gen is not overwritten unintentionally after call to this routine. 
  !< 
  SUBROUTINE user_read_input()
    USE io_3d_vars
    USE parallel                      ! par_size
    IMPLICIT NONE
    INTEGER il_dim, il_ctlSize
    CHARACTER(LEN=255) :: ala_command
    CHARACTER(LEN=ip_snl) :: ala_cand_filename, ala_restart_filename
    CHARACTER (LEN=4) :: p_string
    REAL (pr), ALLOCATABLE, DIMENSION(:) :: rla_ctl_tmp

    lg_restart = .FALSE.
    CALL debug(par_rank, 'Entering user_read_input, process ')
    
    j_mx_adjoint = j_mx
    call input_integer ('j_mx_adjoint',j_mx_adjoint,'default', ' j_mx_adjoint: j_IC is set to j_mx_adjoint')
    !call input_real ('nu',nu,'stop', ' nu: viscosity coefficient')
    call input_logical ('lg_adapt_on_adjoint', lg_adapt_on_adjoint, 'stop', ' : adap grid based on adjoint var when solving adjoint pb?')
    call input_logical ('lg_adapt_on_direct',  lg_adapt_on_direct,  'stop', ' : adap grid based on direct var when solving adjoint pb?')
    call input_real ('rg_alpha', rg_alpha, 'stop', ' : amplitude of the bell shape')
    call input_real ('rg_omega', rg_omega, 'stop', ' : 2*PI')
    
    !Parameters relative to the dimentionality of the problem
    SELECT CASE(dim)
       CASE (2)
          CALL input_real_vector ('rga_x0_2d', rga_x0, dim, 'stop', ' : location of the center of the bell shape')
          call input_real ('rg_v0_2d',    rg_v0,    'stop', ' : advection velocity factor')
          call input_real ('rg_sigma_2d', rg_sigma, 'stop', ' : std deviation of the Gaussian function')
          call input_real ('nu_2d', nu, 'stop', ' : diffusion coef')
          call input_integer ('ig_av_config_2d',ig_av_config,'default', ' ig_av_config: configuration of the advection velocity')
       CASE (3)
          CALL input_real_vector ('rga_x0_3d', rga_x0, dim, 'stop', ' : location of the center of the bell shape')
          call input_real ('rg_v0_3d',    rg_v0,    'stop', ' : advection velocity factor')
          call input_real ('rg_sigma_3d', rg_sigma, 'stop', ' : std deviation of the Gaussian function')
          call input_real ('nu_3d', nu, 'stop', ' : diffusion coef')
          call input_integer ('ig_av_config_3d',ig_av_config,'default', ' ig_av_config: configuration of the advection velocity')
    END SELECT
    
    j_ctl = min(j_zn,j_mx)
    !call input_real ('', , 'stop', ' : obs error standard deviation')
    !CALL input_real_vector ('', , , 'stop', ' : control params')!model

    CALL read_ep_data(tg_ep, EP_DATA, RTIME_FILE)
    il_ctlSize = SIZE(tg_ep%ra_ctl)
    ALLOCATE( rga_ctl(il_ctlSize) )
    IF( par_size ==1 )THEN !sequential mode
       aga_restart_filename = TRIM(res_path)//TRIM(tg_ep%restart_filename)//'.res'
       aga_cand_filename    = TRIM(res_path)//TRIM(tg_ep%restart_candidate)//'.res'
       aga_direct_fileName  = TRIM(res_path)//TRIM(tg_ep%direct_fileName)//'.res'
    ELSE!=> parallel mode
       WRITE(p_string, '(I4)') par_rank
       p_string = ADJUSTL(p_string)
       aga_restart_filename = TRIM(res_path)//TRIM(tg_ep%restart_filename)//'.p'//TRIM( p_string )//'.res'
       aga_restart_com_fName= TRIM(res_path)//TRIM(tg_ep%restart_filename)//'.com.res'
       aga_cand_filename    = TRIM(res_path)//TRIM(tg_ep%restart_candidate)//'.p'//TRIM( p_string )//'.res'
       aga_cand_com_fName   = TRIM(res_path)//TRIM(tg_ep%restart_candidate)//'.com.res'
       aga_direct_fileName  = TRIM(res_path)//TRIM(tg_ep%direct_fileName)//'.p'//TRIM( p_string )//'.res'
       aga_direct_com_fName = TRIM(res_path)//TRIM(tg_ep%direct_fileName)//'.com.res'
    END IF
    
    IF(tg_ep%l_useGD.AND.tg_ep%l_run_from_ctl) THEN
       CALL gd_var_change(tg_ep%ra_ctl, rga_ctl)
    ELSE
       rga_ctl = tg_ep%ra_ctl
    END IF
    !!!!!!!!!!!!!!!!!!!!!!
    IF ( tg_ep%aa_solver_action == MAKE_ADMT ) THEN
       !saving ctl
       ALLOCATE( rla_ctl_tmp( SIZE(tg_ep%ra_ctl) ) )
       rla_ctl_tmp  = tg_ep%ra_ctl
       tg_ep%ra_ctl = rga_ctl
       CALL write_ep_data( tg_ep, ACTL_DATA, OUTPUT_FILE, 'Solver, CTL for ADMT' )
       !restauring CTL
       tg_ep%ra_ctl = rla_ctl_tmp
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!
    tg_obs%obs_fName  = make_fileName(tg_ep, OBS_DATA,  INPUT_FILE)
    tg_obs%ogap_fName = make_fileName(tg_ep, OGAP_DATA, RTIME_FILE)

    
    lg_direct  = .FALSE.
    lg_adjoint = .FALSE.
    SELECT CASE(tg_ep%aa_solver_action)
    CASE (RUN_NOTHING)!Nothing to do, the minimizer is moving to the next iteration
       !copy the candidate restart file to the restart file
       !todo
       !*****This is not the good way to do it, should take into account some other informations from the minimizer
       !*****If there is only one simulation in the iteration, this is finest
       !*****Otherwise, should run on more simulation (solver)    
       CALL debug('', 'In user_read_input: copying the restart file here')
       ala_command = 'cp -p '//TRIM(aga_cand_filename)//' '//TRIM(aga_restart_filename)
       IF( (par_size>1).AND.(par_rank==0) ) THEN
          CALL SYSTEM( 'cp -p '//TRIM(aga_cand_com_fName)//' '//TRIM(aga_restart_com_fName) )
       END IF
       !CALL debug(ala_command, 'Running: ')
       CALL SYSTEM( ala_command )
       !CALL debug('', 'Done')
       IF(par_rank==0) CALL debug('', 'In user_read_input: minimizer moving to the next iteration')
       CALL parallel_finalize
       STOP
    CASE (RUN_COSTGRAD) !cost function and its gradient recquired simultaneously
       lg_direct  = .TRUE.
       lg_adjoint = .TRUE.
       CALL init_costgrad()
       file_gen = TRIM(file_gen)//'costgrad.'
    CASE (RUN_COST) !cost function recquired
       lg_direct = .TRUE.
       CALL init_cost()
       file_gen = TRIM(file_gen)//'cost.'
    CASE (RUN_GRADIENT) !gradient of the cost function recquired
       lg_adjoint = .TRUE.
       CALL init_grad()
       CALL adjoint_zeroing()
       file_gen = TRIM(file_gen)//'grad.'

       !setting the initialization from the solution of the direct problem
       !set the restarts parameters
       IC_restart_mode = 3! restart from IC
       IC_restart_station = 0! should set other value to keep trace of the evolution
       IF( par_size == 1 )THEN !sequential mode
          IC_filename = TRIM(res_path)//TRIM(tg_ep%direct_fileName)//'.res'!
       ELSE
          IC_filename = TRIM(res_path)//TRIM(tg_ep%direct_fileName)//'.p0.res'!
       END IF
       lg_restart = .TRUE.
    CASE (MAKE_OBS) !make observations for twin experiments
       lg_direct = .TRUE.
       CALL init_twin_obs
       file_gen = TRIM(file_gen)//'obs.'
    CASE (MAKE_CTL)! generates a control vector, the true trajectory and the true obs
      CALL make_zero_ctl(tg_ep)
      CALL parallel_finalize
      STOP
    CASE (MAKE_BG)! generates and saves a zero background control vector
      CALL make_zero_bg(tg_ep)
      !CALL make_bg_from_ctl(tg_ep, 0.8)
      CALL parallel_finalize
      STOP
    CASE (RUN_DIRECT) ! direct model run recquired 
       lg_direct = .TRUE.
       file_gen = TRIM(file_gen)//'direct.'
       !CALL debug(TRIM(tg_ep%aa_solver_action), 'In user_read_input: nothing to do for: '
    CASE (MAKE_ADMT) ! analysed direct model trajectory, after minimization 
       lg_direct = .TRUE.
       file_gen = TRIM(file_gen)//'analysis.'
       !CALL debug(TRIM(tg_ep%aa_solver_action), 'In user_read_input: nothing to do for: '
    CASE (RUN_ADJOINT)! adjoint model run recquired
       lg_adjoint = .TRUE.
       !CALL init_adjoint()
       !CALL adjoint(...)
       file_gen = TRIM(file_gen)//'adjoint.'
       !CALL debug(TRIM(tg_ep%aa_solver_action), 'In user_read_input: nothing to do for: ')
    CASE DEFAULT
       CALL debug(TRIM(tg_ep%aa_solver_action), 'In user_read_input: nothing to do for: ')
       CALL parallel_finalize
       STOP
    END SELECT
    
    IF(.not.lg_direct .AND. lg_adjoint) j_IC = MAX(j_zn,j_mx_adjoint)

    IF( ( (tg_ep%aa_solver_action==RUN_COST).OR.(tg_ep%aa_solver_action==RUN_COSTGRAD) ).AND.&
        ( tg_ep%l_restart_from_previous.AND.(.NOT.tg_ep%l_first_simul) )&
      )THEN
       !set the restarts parameters
       IC_restart_mode = 3! restart from IC
       IC_restart_station = 0! should set other value to keep trace of the evolution
       IF( par_size == 1 )THEN !sequential mode
          IC_filename = TRIM(res_path)//TRIM(tg_ep%restart_fileName)//'.res'!
       ELSE
          IC_filename = TRIM(res_path)//TRIM(tg_ep%restart_fileName)//'p0.res'!
       END IF
       lg_restart = .TRUE.
    END IF

    IF(lg_restart)THEN
       !--------------------------- RESETTING PARAMETERS ASSOCIATED WITH ic_restart_mode  ------------------- 
       ! for transition to IC_restart_mode
       IC_restart = .FALSE.
       IC_from_file = .FALSE.
       IF (IC_restart_mode /= -1) THEN ! IC_restart_mode is not used in .inp file
          IF (par_rank.EQ.0) THEN
             PRINT *, ' '
             PRINT *, 'Thank you for using IC_restart_mode instead of IC_restart/IC_restart_from_file'
             PRINT *, ' '
          END IF
          SELECT CASE (IC_restart_mode)
          CASE (1:2)                          ! hard/soft restart
             IC_restart = .TRUE.
          CASE (3)                            ! restart from IC
             IC_from_file = .TRUE.
          END SELECT
       END IF
       IF (IC_restart.OR.IC_from_file) THEN
          IC_file_fmt = 0
          IC_adapt_grid = .TRUE.
          IC_SINGLE_LOOP = .FALSE.
          IC_filename = ADJUSTL(IC_filename)
       END IF
    END IF

    CALL debug(par_rank, 'Exiting user_read_input ***********************************')
  END SUBROUTINE user_read_input
  
  !> Set the exact solution for comparison to the simulated solution
  !!
  !! \param[in,out] u state variable for exact solution, array to fill in the exact solution
  !! \param[in] nlocal number of active wavelets   
  !! \param[in] ne_local total number of equations
  !! \param[in] t time of current time step 
  !! \param[in] l_n_var_exact_soln_index index into the elements of u for which we need to find the exact solution
  !<
  SUBROUTINE  user_exact_soln (u, nlocal,  t_local, l_n_var_exact_soln)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) ::  t_local
    REAL (pr), DIMENSION (nlocal,n_var_exact), INTENT (INOUT) :: u
    LOGICAL , INTENT (IN) :: l_n_var_exact_soln(n_var)

  END SUBROUTINE  user_exact_soln

  !< \brief Computes initial condition for the first variable of the direct problem
  !! \param[in] id_ng size of the state vector
  !! \details
  !! uses the global variables x and rg_sigma
  !>
  SUBROUTINE direct_initial_condition(rda_u)
    IMPLICIT NONE
    REAL(KIND=pr), INTENT(INOUT), DIMENSION(:) :: rda_u
    INTEGER :: idim
    REAL(KIND=pr), DIMENSION(nwlt) :: dummy
    CALL extract_ctl_idx(j_lev)
    IF(.NOT.lg_restart)THEN
       IF(tg_ep%l_run_from_ctl)THEN
          !CALL debug(10, 'In direct_initial_condition: here ...... ')
          WHERE(iga_ctl_idx > 0) 
             rda_u(iga_ctl_idx) = rga_ctl + tg_ep%ra_b_ctl
          END WHERE
          !CALL debug(20, 'In direct_initial_condition: here ...... ')
       ELSE
          dummy = 0.0_pr
          DO idim = 1, dim-1
             dummy = dummy + ( x(:, idim) - rga_x0(idim) )**2
          END DO
          WHERE(iga_ctl_idx > 0)
             rda_u(iga_ctl_idx) = rg_alpha*exp(-dummy(iga_ctl_idx)/(2.0_pr * rg_sigma**2) )
          END WHERE
       END IF
    END IF
  END SUBROUTINE direct_initial_condition

  !> Computes initial condition, this is the system state at initial
  ! time
  !! \param[in, out] u sytem state
  !! \param[in] nlocal number of grid point in each single variable
  !!  of the system state
  !! \param[in] ne_local number of equation in the system of PDE
  !! \param[in] t_local time associated to the initial condition
  !! \param[in] scl [to be specified]
  !! \param[in] scl_fltwt [to be specified]
  !! \param[in, out] iter [to be specified]
  !!
  !! \detail here this routine is used as the driver for the
  !!  minimizer algorithm
  !<
  SUBROUTINE user_initial_conditions (u, nlocal, ne_local, t_local,&
       & scl, scl_fltwt, iter)
    USE parallel
    IMPLICIT NONE
    INTEGER  , INTENT (IN) :: nlocal, ne_local
    INTEGER  , INTENT (INOUT) :: iter ! iteration of call while
    !  adapting initial grid
    REAL (pr), DIMENSION (nlocal,ne_local), INTENT (INOUT) :: u
    REAL (pr), DIMENSION (nlocal,ne_local) :: f
    REAL (pr), INTENT (IN)   :: scl(1:n_var),scl_fltwt
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION(ne_local) :: scl_u
    INTEGER, DIMENSION(ne_local) :: clip
    INTEGER :: eq1, eq2, il_neq!number of equations to be solve, just
    ! for test purposes.
    LOGICAL, SAVE :: ll_fisrtCall = .TRUE.
    clip = 0 !no clipping for Dirichlet BC
    IF (ll_fisrtCall) THEN
       CALL direct_initial_condition( u(:, 1) )
       u(:, 2) = 0.0_pr
    END IF
    ll_fisrtCall = .FALSE.
    scl_u = 1.0_pr

    f = 0.0_pr
    f = Laplace_rhs (u, nlocal, ne_local, HIGH_ORDER)

    IF(lg_direct .AND. lg_adjoint)THEN  !both direct and adjoint
       ! simulation required
       eq1 = 1
       eq2 = 2
    ELSE IF(lg_direct)THEN !only direct simulation required, when
       ! called by the routine <cost>
       eq1 = 1
       eq2 = 1
    ELSE IF(lg_adjoint) THEN !only adjoint required, when called by
       ! the routine <costb>
       eq1 = 2
       eq2 = 2
    ELSE
       CALL debug('', 'Nothing to do********************************')
       CALL parallel_finalize
       STOP
    END IF

    il_neq = eq2 - eq1 + 1
    !CALL debug('', 'In user_initial_conditions : calling Linsolve *
    !**2')
    CALL Linsolve ( u(:, eq1:eq2), f(:, eq1:eq2), tol2, nlocal,&
         & il_neq, clip, Laplace, Laplace_diag, SCL=scl_u )  !
    !CALL debug('', 'In user_initial_conditions : after Linsolve ***')
  END SUBROUTINE user_initial_conditions

  !--********************************  
  ! Arguments
  ! u         - field on adaptive grid
  ! nlocal      - number of active points
  ! ne_local       - number of equations
  ! t         - current time
  !
  ! Global variables used
  ! ifn       - number of equations
  ! x()
  ! xO
  ! y0
  !--******************************** 
  SUBROUTINE user_algebraic_BC (Lu, u, nlocal, ne_local, jlev, meth)
    !--Defines boundary condition type
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (IN)    :: u

    INTEGER :: ie, ii, shift
    REAL (pr), DIMENSION (nlocal,ne_local) :: du, d2u


    !
    ! There are periodic BC conditions
    !


  END SUBROUTINE user_algebraic_BC

  !> \brief
  !!
  !<
  SUBROUTINE user_algebraic_BC_diag (Lu_diag, nlocal, ne_local, jlev,&
       & meth)
    !--Defines boundary condition type
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu_diag
    REAL (pr), DIMENSION (nlocal,dim) :: du, d2u

    !
    ! There are periodic BC conditions
    !

  END SUBROUTINE user_algebraic_BC_diag

  !> \brief
  !!
  !<
  SUBROUTINE user_algebraic_BC_rhs (rhs, ne_local, nlocal, jlev)
    !--Sets rhs for boundary conditions
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: ne_local, nlocal, jlev
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: rhs

    INTEGER :: ie, ii, shift
    !
    ! There are periodic BC conditions
    !

  END SUBROUTINE user_algebraic_BC_rhs

  !> \brief
  !!
  !<
  SUBROUTINE user_project (u, p, nlocal, meth)
    !--Makes u divergence free
    IMPLICIT NONE

    INTEGER, PARAMETER :: ne_local = 1
    INTEGER, INTENT (IN) :: meth, nlocal
    REAL (pr), DIMENSION (nlocal),     INTENT(INOUT) :: p
    REAL (pr), DIMENSION (nlocal,n_integrated), INTENT(INOUT) :: u

    !------------ not used in this case

  END SUBROUTINE user_project

  !> \brief Computes the advection velocity
  !! \param[in] id_jlev level of computation, level refers to wavelet
  !!  decomposition
  !! \param[in] id_nlocal 
  !! \detail computes the advection velocity field on the grid
  !!  associate to jlev wavelet decomposition
  !<
  FUNCTION advection_velocity(id_jlev, id_nlocal)RESULT(rla_v)
    USE constant
    USE parallel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id_jlev
    INTEGER, INTENT(IN) :: id_nlocal
    REAL (pr), DIMENSION ( id_nlocal, dim-1 ) :: rla_v !spatial
    ! advection velocity field. No component in the time direction,
    !  to be adapted according to the number of spatial dimensions
    INTEGER :: ibi, face_type, nloc, j, wlt_type, j_df, k, i, ii, id,&
         & vshift
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    REAL(KIND=pr) :: rl_omega, rl_rvm
    REAL(KIND=pr), DIMENSION(dim) :: rla_vc!cordinate of the vortex center

    rl_omega = 2*rp_pi
  
    i_p_face(0) = 1
    DO ibi=1,dim
       i_p_face(ibi) = i_p_face(ibi-1)*3
    END DO
    SELECT CASE (ig_av_config)
       CASE (2,3)
          rla_vc = ( xyzlimits(1,:) + xyzlimits(2,:) )/2.0_pr
          rl_rvm = SQRT( ( xyzlimits(2,1) - xyzlimits(1,1) )**2 + ( xyzlimits(2,2) - xyzlimits(1,2) )**2 )/10.0_pr
          !CALL debug(rla_vc, 'In advection_velocity, rla_vc = ')
          !CALL debug(rl_rvm, 'In advection_velocity, rl_vm  = ')
       CASE DEFAULT
          
    END SELECT
    !computing advection velocity
    DO j = 1, id_jlev
       DO wlt_type = MIN(j-1,1),2**dim-1
          DO face_type = 0, 3**dim - 1
             !face = INT(MOD(face_type,i_p_face(1:dim))
             !/i_p_face(0:dim-1))-1
             DO j_df = j, j_lev
                DO k = 1, indx_DB(j_df,wlt_type,face_type,j)%length  
                   i = indx_DB(j_df,wlt_type,face_type,j)%p(k)%i 
                   ii = i+ indx_DB(id_jlev,wlt_type,face_type,j)%shift
                   id = 1
                   vshift=(id-1)*id_nlocal
                   SELECT CASE(ig_av_config)
                      CASE(1)
                         rla_v(vshift+ii, :) = rg_v0*sin( rl_omega&
                              &*x(i,dim) )
                      CASE (2)!vortex for 2D+time
                         IF(dim==3) THEN
                            rla_v(vshift+ii, 1) = u_init( x(i,1), x(i,2), rla_vc(1), rla_vc(2), rl_rvm)
                            rla_v(vshift+ii, 2) = v_init( x(i,1), x(i,2), rla_vc(1), rla_vc(2), rl_rvm)
                         ELSE
                            CALL debug('In advection_velocity: vortex configuration of advection velocity defined only for 2D+time')
                            CALL parallel_finalize
                            STOP
                         END IF
                      CASE (3)!rotation for 2D+time
                         IF(dim==3) THEN
                            rla_v(vshift+ii, 1) = -( x(i,2) -  rla_vc(2) )*0.5
                            rla_v(vshift+ii, 2) = ( x(i,1) -  rla_vc(1) )*0.5
                         ELSE
                            CALL debug('In advection_velocity: rotation configuration of advection velocity defined only for 2D+time')
                            CALL parallel_finalize
                            STOP
                         END IF
                      CASE DEFAULT
                         CALL debug('In advection_velocity : Undefined configuration for advection velocity')
                         CALL parallel_finalize
                         STOP
                   END SELECT
                END DO
             END DO
          END DO
       END DO
    END DO

  END FUNCTION advection_velocity

  !> \brief direct model equation
  !! \param[in] id_jlev level of computation, level refers to wavelet decomposition
  !! \param[in] rda_u state variable
  !! \param[in] id_ng size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] id_neq number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] id_meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION direct(id_jlev, rda_u, id_ng, id_neq, id_meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id_jlev, id_ng, id_neq, id_meth
    REAL (pr), DIMENSION ( id_ng*id_neq ), INTENT (INOUT) :: rda_u
    REAL (pr), DIMENSION ( id_ng*id_neq ) :: direct
    REAL (pr), DIMENSION ( id_ng, dim-1 ) :: rla_v!spatial advection velocity field. No component in the time direction, to be adapted according to the number of spatial dimensions

    INTEGER :: i, ii, ie, shift, idim
    INTEGER :: meth_central, meth_backward, meth_forward
    REAL (pr), DIMENSION (id_neq,id_ng, dim) :: dub, du, d2u
    REAL (pr), DIMENSION (id_neq,id_ng, dim) :: du_dummy! passed when only calculating 1st derivative.
    INTEGER :: face_type, nloc, j, wlt_type, j_df, k
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    INTEGER :: ibi

    !CALL debug('', ' direct starting <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
    meth_central  = id_meth + BIASING_NONE
    meth_backward = id_meth + BIASING_BACKWARD
    meth_forward  = id_meth + BIASING_FORWARD

    CALL c_diff_fast (rda_u, dub, du_dummy, id_jlev, id_ng, meth_backward, 10, id_neq, 1, id_neq)!first order derivatives for direct model, for time dimension
    CALL c_diff_fast (rda_u, du , d2u     , id_jlev, id_ng, meth_central , 11, id_neq, 1, id_neq)!first order der and laplacian for spatial dimensions
    rla_v = advection_velocity(id_jlev, id_ng)
    !****first equation of the  direct model
    ie = 1
    shift=(ie-1)*id_ng
    direct(shift+1:shift+Nwlt_lev(id_jlev,1)) = dub(ie, 1:Nwlt_lev(id_jlev,1), dim) !time dimension
    DO idim = 1, dim-1
       direct(shift+1:shift+Nwlt_lev(id_jlev,1)) = direct(shift+1:shift+Nwlt_lev(id_jlev,1)) &
            - nu*d2u(ie, 1:Nwlt_lev(id_jlev,1), idim)&
            + rla_v( 1:Nwlt_lev(id_jlev,1), idim )*du( ie ,1:Nwlt_lev(id_jlev,1),idim)
    END DO

    !--Boundary points
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1

       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
          CALL get_all_indices_by_face (face_type, id_jlev, nloc, iloc)
          
          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                direct(shift+iloc(1:nloc)) = rda_u(shift+iloc(1:nloc))    !Boundaries conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == -1 )THEN ! IC only 
                direct(shift+iloc(1:nloc)) = rda_u(shift+iloc(1:nloc))    !Initial Conditions, no boundary points
             END IF
          END IF

       END IF
    END DO
    !CALL debug('', ' direct ending <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
  END FUNCTION direct

  !> \brief adjoint model equation
  !! \param[in] id_jlev level of computation, level refers to wavelet decomposition
  !! \param[in] rda_u state variable
  !! \param[in] id_ng size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] id_neq number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] id_meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION adjoint(id_jlev, rda_u, id_ng, id_neq, id_meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id_jlev, id_ng, id_neq, id_meth
    REAL (pr), DIMENSION ( id_ng*id_neq ), INTENT (INOUT) :: rda_u
    REAL (pr), DIMENSION ( id_ng*id_neq ) :: adjoint
    REAL (pr), DIMENSION ( id_ng, dim-1 ) :: rla_v!spatial advection velocity field. No component in the time direction, to be adapted according to the number of spatial dimensions

    INTEGER :: i, ii, ie, shift, idim
    INTEGER :: meth_central, meth_backward, meth_forward
    REAL (pr), DIMENSION (id_neq,id_ng, dim) :: duf, du, d2u
    REAL (pr), DIMENSION (dim-1 ,id_ng, dim) :: dv!
    REAL (pr), DIMENSION (dim-1,id_ng, dim) :: du_dummy! passed when only calculating 1st derivative.
    INTEGER :: face_type, nloc, j, wlt_type, j_df, k
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    INTEGER :: ibi

    !CALL debug('', ' adjoint starting vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ')
    meth_central  = id_meth + BIASING_NONE
    meth_backward = id_meth + BIASING_BACKWARD
    meth_forward  = id_meth + BIASING_FORWARD

    CALL c_diff_fast (rda_u, duf, du_dummy, id_jlev, id_ng, meth_forward , 10, id_neq, 1, id_neq)!first order derivatives for adjoint model, for time dimension
    CALL c_diff_fast (rda_u, du , d2u     , id_jlev, id_ng, meth_central , 11, id_neq, 1, id_neq)!first order der and laplacian for spatial dimensions
    
    rla_v = advection_velocity(id_jlev, id_ng)

    CALL c_diff_fast (rla_v, dv,  du_dummy, id_jlev, id_ng, meth_central , 10, dim-1, 1, dim-1)

    !****first equation of the adjoint model
    ie = 1
    shift=(ie-1)*id_ng
    adjoint(shift+1:shift+Nwlt_lev(id_jlev,1)) = duf(ie, 1:Nwlt_lev(id_jlev,1), dim) !time dimension
    DO idim = 1, dim-1
       adjoint(shift+1:shift+Nwlt_lev(id_jlev,1)) = adjoint(shift+1:shift+Nwlt_lev(id_jlev,1)) &
            + nu*d2u(ie, 1:Nwlt_lev(id_jlev,1), idim)&
            + rla_v( 1:Nwlt_lev(id_jlev,1), idim )*du( ie ,1:Nwlt_lev(id_jlev,1),idim)&
            + rda_u(shift+1:shift+Nwlt_lev(id_jlev,1))*dv( idim, 1:Nwlt_lev(id_jlev,1), idim)
    END DO
    !--Boundary points
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1

       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
          CALL get_all_indices_by_face (face_type, id_jlev, nloc, iloc)

          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                adjoint(shift+iloc(1:nloc)) = rda_u(shift+iloc(1:nloc))    !Boundaries conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == 1 )THEN ! IC only 
                adjoint(shift+iloc(1:nloc)) = rda_u(shift+iloc(1:nloc))   !Initial Conditions, no boundary points
             END IF
          END IF

       END IF
    END DO
    !CALL debug('', ' adjoint ending vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ')
  END FUNCTION adjoint
  
  !> \brief Computes the left hand side of time indepemdent problem written as F(u) = 0
  !! \param[in] jlev level of computation, level refers to wavelet decomposition
  !! \param[in] u state variable
  !! \param[in] nlocal size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] ne_local number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION Laplace (jlev, u, nlocal, ne_local, meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: jlev, nlocal, ne_local, meth
    REAL (pr), DIMENSION ( nlocal*ne_local ), INTENT (INOUT) :: u
    REAL (pr), DIMENSION ( nlocal*ne_local ) :: Laplace
    INTEGER :: il_dshift, il_ashift
    
    il_dshift = -nlocal
    !CALL debug('', 'Entering Laplace ************************************************')
    !direct model
    IF(lg_direct)THEN
       il_dshift = 0
       Laplace(il_dshift+1:il_dshift+nlocal) = direct(jlev, u(il_dshift+1:il_dshift+nlocal), nlocal, 1, meth)
    ELSE
       il_dshift = -nlocal
    END IF

    !adjoint model
    IF(lg_adjoint)THEN
       il_ashift = il_dshift + nlocal
       Laplace(il_ashift+1:il_ashift+nlocal) = adjoint(jlev, u(il_ashift+1:il_ashift+nlocal), nlocal, 1, meth)
    END IF
    
    !CALL debug('', 'Ending Laplace ************************************************')
  END FUNCTION Laplace

  !> \brief Computes the diagonal part of direct model equations
  !! \param[in] id_jlev level of computation, level refers to wavelet decomposition
  !! \param[in] id_ng size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] id_neq number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] id_meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION direct_diag(id_jlev, id_ng, id_neq, id_meth)
    USE precision
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: id_jlev, id_ng, id_neq, id_meth
    REAL (pr), DIMENSION (1:id_ng*id_neq) :: direct_diag
    REAL (pr), DIMENSION ( id_ng, dim-1 ) :: rla_v

    INTEGER :: i, ii, ie, shift, idim
    INTEGER :: meth_central, meth_backward, meth_forward
    REAL (pr), DIMENSION (id_ng,dim) :: dub, du, d2u, du_dummy
!    REAL (pr), DIMENSION (:, :), POINTER :: du_dummy ! passed when only calculating 1st derivative.

    INTEGER :: face_type, nloc, j, wlt_type, j_df, k
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    INTEGER :: ibi

    !CALL debug('', ' direct_diag starting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    meth_central  = id_meth + BIASING_NONE
    meth_backward = id_meth + BIASING_BACKWARD
    meth_forward  = id_meth + BIASING_FORWARD

    CALL c_diff_diag (dub, du_dummy, id_jlev, id_ng, meth_backward, meth_backward, 10)
    CALL c_diff_diag (du , d2u     , id_jlev, id_ng, meth_central , meth_central , 11)

    rla_v = advection_velocity(id_jlev, id_ng)

    !****first equation : direct model
    ie = 1
    shift=(ie-1)*id_ng
    direct_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) = dub(1:Nwlt_lev(id_jlev,1),dim)
    DO idim = 1,dim-1
       direct_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) = direct_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) &
            - nu*d2u(1:Nwlt_lev(id_jlev,1),idim)&
            + rla_v( shift+1:shift+Nwlt_lev(id_jlev,1), idim )*du( 1:Nwlt_lev(id_jlev,1), idim)
    END DO
    !--Boundary points
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO

    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points

          CALL get_all_indices_by_face (face_type, id_jlev, nloc, iloc)

          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                direct_diag(shift+iloc(1:nloc)) = 1.0_pr    !Boundaries conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == -1 )THEN ! IC only 
                direct_diag(shift+iloc(1:nloc)) = 1.0_pr    !Initial Conditions, no boundary points
             END IF
          END IF
       END IF
    END DO

    !CALL debug('', ' Ending direct_diag >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
  END FUNCTION direct_diag

  !> \brief Computes the diagonal part of adjoint model equations
  !! \param[in] id_jlev level of computation, level refers to wavelet decomposition
  !! \param[in] id_ng size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] id_neq number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] id_meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION adjoint_diag(id_jlev, id_ng, id_neq, id_meth)
    USE precision
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: id_jlev, id_ng, id_neq, id_meth
    REAL (pr), DIMENSION (id_ng*id_neq) :: adjoint_diag
    REAL (pr), DIMENSION ( id_ng, dim ) :: rla_v

    INTEGER :: i, ii, ie, shift, idim
    INTEGER :: meth_central, meth_backward, meth_forward
    REAL (pr), DIMENSION (id_ng,dim) :: duf_diag, du_diag, d2u_diag, du_dummy

    INTEGER :: face_type, nloc, j, wlt_type, j_df, k
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    INTEGER :: ibi

    !CALL debug('', ' adjoint_diag starting ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    meth_central  = id_meth + BIASING_NONE
    meth_backward = id_meth + BIASING_BACKWARD
    meth_forward  = id_meth + BIASING_FORWARD

    CALL c_diff_diag (duf_diag, du_dummy, id_jlev, id_ng, meth_forward, meth_forward, 10)!first order derivatives for adjoint model
    CALL c_diff_diag (du_diag, d2u_diag     , id_jlev, id_ng, meth_central, meth_central, 11)
    
    rla_v = advection_velocity(id_jlev, id_ng)
    !****second equation : adjoint model
    ie = 1
    shift=(ie-1)*id_ng
    adjoint_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) = duf_diag(1:Nwlt_lev(id_jlev,1),dim) !time dimension
    DO idim = 1,dim-1
       adjoint_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) = adjoint_diag(shift+1:shift+Nwlt_lev(id_jlev,1)) &
            + nu*d2u_diag(1:Nwlt_lev(id_jlev,1),idim)&
            + rla_v( 1:Nwlt_lev(id_jlev,1), idim )*du_diag( 1:Nwlt_lev(id_jlev,1), idim)&
            + du_diag(1:Nwlt_lev(id_jlev,1), idim)
    END DO
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
          
          CALL get_all_indices_by_face (face_type, id_jlev, nloc, iloc)
          
          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                adjoint_diag(shift+iloc(1:nloc)) = 1.0_pr   !Boundaries conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == 1 )THEN ! IC only 
                adjoint_diag(shift+iloc(1:nloc)) = 1.0_pr    !Initial Conditions, t=t_end
             END IF
          END IF
       END IF
    END DO
       
    !CALL debug('', ' Ending adjoint_diag ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
  END FUNCTION adjoint_diag
  
  !> \brief Computes the diagonal of the left hand side of time indepemdent problem written as F(u) = 0
  !! \param[in] jlev level of computation, level refers to wavelet decomposition
  !! \param[in] nlocal size of a algebraic component of the state vector, this is the number of grid point (discretization)
  !! \param[in] ne_local number of algebraic variables in the system state, this is also the number of equations in the direct model
  !! \param[in] meth method to use for computing derivatives
  !! \detail the goal is to separate direct and adjoint model so that each can be checked separately
  !<
  FUNCTION Laplace_diag (jlev, nlocal, ne_local, meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: jlev, nlocal, ne_local, meth
    REAL (pr), DIMENSION (1:nlocal*ne_local) :: Laplace_diag
    INTEGER :: il_dshift, il_ashift

    !CALL debug(ne_local, 'Entering Laplace_diag ***************************** ne_local = ')
    il_dshift = -nlocal

    !direct model
    IF(lg_direct)THEN
       il_dshift = 0
       Laplace_diag(il_dshift+1:il_dshift+nlocal) = direct_diag (jlev, nlocal, 1, meth)
    END IF

    !adjoint model
    IF(lg_adjoint)THEN
       il_ashift = il_dshift + nlocal
       Laplace_diag(il_ashift+1:il_ashift+nlocal) = adjoint_diag(jlev, nlocal, 1, meth)
    END IF

   !CALL debug('', 'Ending Laplace_diag ************************************************')
  END FUNCTION Laplace_diag

!/!\attention, u is the full size u of size nlocal*ne_local
  FUNCTION direct_rhs(rda_u, id_ng, id_neq, id_meth)
    USE penalization
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: id_ng, id_neq, id_meth
    REAL (pr), DIMENSION (id_ng* id_neq), INTENT (INOUT) :: rda_u
    REAL (pr), DIMENSION (id_ng) :: direct_rhs

    INTEGER :: i, ii, ie, meth, idim
    INTEGER :: meth_central, meth_backward, meth_forward 
    REAL (pr), DIMENSION (id_ng, dim) :: Xdel
    REAL (pr), DIMENSION (id_neq, id_ng, dim) :: du, d2u
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    
    !CALL debug('', ' Entering direct_rhs  *******************************')

    direct_rhs  = 0.0_pr
    !--Boundary points
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points         
          CALL get_all_indices_by_face (face_type, j_lev, nloc, iloc)

          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                direct_rhs(iloc(1:nloc)) = 0.0_pr    !Boundaries conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == -1 )THEN ! IC only 
                direct_rhs(iloc(1:nloc)) = 0.0_pr
                DO idim=1,dim-1
                    direct_rhs(iloc(1:nloc)) = direct_rhs(iloc(1:nloc)) + (x(iloc(1:nloc), idim) - rga_x0(idim))**2
                 END DO
                direct_rhs(iloc(1:nloc)) = rg_alpha*exp( - direct_rhs(iloc(1:nloc))/(2.0_pr * rg_sigma**2) )!Initial Conditions, no boundary points
             END IF
          END IF
       END IF
    END DO
    !when runing from ctl
    IF(tg_ep%l_run_from_ctl)THEN
       CALL extract_ctl_idx(j_lev)
       WHERE(iga_ctl_idx > 0) 
          direct_rhs(iga_ctl_idx) = rga_ctl + tg_ep%ra_b_ctl
       END WHERE
    END IF

    !CALL debug('', ' Ending direct_rhs  *******************************')
    !IF(.NOT.lg_m1qn3_initialization) CALL dpause()
  END FUNCTION direct_rhs


  FUNCTION adjoint_rhs(rda_u, id_ng, id_neq, id_meth)
    USE penalization
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: id_ng, id_neq, id_meth
    REAL (pr), DIMENSION (id_ng* id_neq), INTENT (INOUT) :: rda_u
    REAL (pr), DIMENSION (id_ng) :: adjoint_rhs

    INTEGER :: i, ii, ie, meth, idim
    INTEGER :: meth_central, meth_backward, meth_forward 
    REAL (pr), DIMENSION (id_ng, dim) :: Xdel
    REAL (pr), DIMENSION (id_neq, id_ng, dim) :: du, d2u
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    INTEGER :: il_nobs
    
    !CALL debug('', ' Entering adjoint_rhs  *******************************')
    !initialization
    adjoint_rhs = 0.0_pr
    il_nobs = get_obsSize(tg_obs)
    !--Boundary points
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
       IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
          CALL get_all_indices_by_face (face_type, j_lev, nloc, iloc)

          IF(nloc > 0 ) THEN ! D + t = dim
             IF( .NOT.(SUM(ABS(face(1:dim-1))) == 0))THEN ! BC only 
                adjoint_rhs(iloc(1:nloc)) = 0.0_pr      !Boundary conditions
             ELSE IF( SUM(ABS(face(1:dim-1))) == 0 .AND. face(dim) == 1 )THEN ! IC at at t=t_end
                adjoint_rhs(iloc(1:nloc)) = 0.0_pr      !Initial conditions at t=t_end
             END IF
          END IF
       END IF
    END DO
    CALL debug('', ' Forcing term based on obs gap')
    IF(tg_ep%aa_solver_action==RUN_COSTGRAD)THEN
       CALL compute_obsgap()! compute_obsgap makes call to compute_obs_idx
    ELSE!!the obsgap is read from file
       CALL compute_obs_idx()
    END IF
    
    WHERE(iga_obs_idx>0)
       adjoint_rhs(iga_obs_idx) = tg_obs%ra_obsgap*tg_obs%ra_Rm1
    END WHERE
    !CALL debug('', ' Ending adjoint_rhs  *******************************')

  END FUNCTION adjoint_rhs

  FUNCTION Laplace_rhs(u, nlocal, ne_local, meth_in)
    USE penalization
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal, ne_local, meth_in
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: u
    REAL (pr), DIMENSION (nlocal, ne_local) :: Laplace_rhs

    INTEGER :: i, ii, ie, meth, idim
    INTEGER :: meth_central, meth_backward, meth_forward 
    REAL (pr), DIMENSION (nlocal,dim) :: Xdel
    REAL (pr), DIMENSION (ne_local ,nlocal,dim) :: du, d2u
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    REAL(pr) :: r(nlocal), Q(nlocal)
    REAL (pr), DIMENSION(dim) :: x0
    
    !CALL debug('j_lev', 'Entering Laplace_rhs  ******************************* j_lev = ')
    !CALL debug( SHAPE(Laplace_rhs), 'SHAPE(Laplace_rhs) = ')
    !CALL debug('', 'In Laplace_rhs, calling direct_rhs')
    IF(lg_direct )  Laplace_rhs(:, 1) = direct_rhs (u, nlocal, 1, meth_in)
    !CALL debug('', 'In Laplace_rhs, calling adjoint_rhs')
    IF(lg_adjoint) Laplace_rhs(:, 2) = adjoint_rhs(u, nlocal, 1, meth_in)
    !CALL debug('', 'In Laplace_rhs, end of  adjoint_rhs')

    !CALL debug('', 'Laplace_rhs ending ********************************************')
  END FUNCTION Laplace_rhs

  !> \brief Compute the right hand side of the PDE
  !!
  !<
  FUNCTION user_rhs (u_integrated,p)
    IMPLICIT NONE
    REAL (pr), DIMENSION (ng,ne), INTENT(IN) :: u_integrated !1D flat version of integrated variables without BC 
    REAL (pr), DIMENSION (ng), INTENT(IN) :: p
    REAL (pr), DIMENSION (n) :: user_rhs

    INTEGER :: ie, shift
    INTEGER, PARAMETER :: meth=1
    REAL (pr), DIMENSION (ne,ng,dim) :: du, d2u
    REAL (pr), DIMENSION (ng,dim)    :: dp
  
    user_rhs = 0.0_pr
    !--Set operator on boundaries
!!$    IF( Nwlt_lev(j_lev,1) >  Nwlt_lev(j_lev,0) ) &
!!$         CALL user_algebraic_BC_rhs (user_rhs, ne, ng, j_lev)
    
  END FUNCTION user_rhs

  ! find Jacobian of Right Hand Side of the problem
  FUNCTION user_Drhs (u, u_prev_timestep_loc, meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: meth
    !u_prev_timestep passed local to cast it into 2dim array to work with vector derivatives
    REAL (pr), DIMENSION (ng,ne) :: u, u_prev_timestep_loc
    REAL (pr), DIMENSION (n) :: user_Drhs

    INTEGER :: ie, shift
    REAL (pr), DIMENSION (2*ne,ng,dim) :: du  ! 1st derivatives for u (in du(1:ne,:,:)) and u_prev_timestep (in du(ne+1:2*ne,:,:))
    REAL (pr), DIMENSION (ne  ,ng,dim) :: d2u ! 2nd derivatives for u 
    REAL (pr), DIMENSION (ne  ,ng,dim) :: du_dummy ! passed when only calculating 1st derivative.
    !Find batter way to do this!! du_dummy with no storage..
    
    user_Drhs = 0.0_pr
  END FUNCTION user_Drhs


  !
  ! Uses u_prev_timestep, which is global flat (1D ) array set in time_step_cn()
  !
  FUNCTION user_Drhs_diag (meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: meth
    REAL (pr), DIMENSION (n) :: user_Drhs_diag

    INTEGER :: ie, shift
    REAL (pr), DIMENSION (ng,dim) :: du, d2u
    REAL (pr), DIMENSION (ne,ng,dim) :: du_prev_timestep
    REAL (pr), DIMENSION (ne,ng,dim) :: du_dummy ! passed when only calculating 1st derivative.

    user_Drhs_diag = 0.0+pr
  END FUNCTION user_Drhs_diag

  !>
  !!
  !<
  FUNCTION user_chi (nlocal, t_local )
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (nlocal) :: user_chi

    user_chi = 0.0_pr
  END FUNCTION user_chi



  !
  ! Calculate any statitics
  !
  ! startup_flag - 0 when adapting to IC, then 1 inmain integration loop
  !
  !
  SUBROUTINE user_stats ( u ,j_mn, startup_flag)
    USE parallel
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: startup_flag
    INTEGER , INTENT (IN) :: j_mn 
    REAL (pr), DIMENSION (nwlt,1:n_var), INTENT (IN) :: u

    !local vars
    TYPE(exchange_param)  :: tl_ep
    CHARACTER(LEN=ip_snl) :: ala_filename, ala_com_fName
    CHARACTER (LEN=4) :: it_string, p_string
    INTEGER :: il_direct_NaNCount, il_direct_InfCount, il_adjoint_NaNCount, il_adjoint_InfCount,&
               il_gradient_NaNCount, il_gradient_InfCount, ibi, ibj
    LOGICAL :: ll_cost_diverged

    ! debug: test database independent interpolation
    IF (test_interpolate) CALL user_interpolate (u)
    !initialization of local variables
    il_direct_NaNCount   = 0
    il_direct_InfCount   = 0
    il_adjoint_NaNCount  = 0
    il_adjoint_InfCount  = 0
    il_gradient_NaNCount = 0
    il_gradient_InfCount = 0
    ll_cost_diverged     = .FALSE.
    IF(startup_flag==-1)THEN
       
       WRITE(it_string,'(I4.4)') (iwrite-1)
       !get the filename
       IF( par_size == 1 )THEN !sequential mode
          ala_filename = TRIM(res_path)//TRIM(file_gen)//it_string//'.res'
       ELSE
          WRITE(p_string, '(I4)') par_rank
          p_string = ADJUSTL(p_string)
          ala_filename = TRIM(res_path)//TRIM(file_gen)//it_string//'.p'//TRIM( p_string )//'.res'
          ala_com_fName= TRIM(res_path)//TRIM(file_gen)//it_string//'.com.res'
       END IF

       !Computing the cost function or the gradient
       IF (j_lev >= j_zn) THEN
          CALL extract_ctl_idx(j_lev)
       END IF
       SELECT CASE(tg_ep%aa_solver_action)
       CASE (RUN_COSTGRAD)! cost and grad simultaneously
          CALL cost(tg_ep%ra_ctl, tg_ep%r_cost)
          CALL adjoint_zeroing()
          CALL costb(tg_ep%ra_ctl, tg_ep%ra_grad, rg_costb)
          CALL compute_grad()
       CASE (RUN_COST) !cost function required
          CALL cost(tg_ep%ra_ctl, tg_ep%r_cost)
       CASE (RUN_GRADIENT) !gradient of the cost function required
          CALL adjoint_zeroing()
          CALL costb(tg_ep%ra_ctl, tg_ep%ra_grad, rg_costb)
          CALL compute_grad()
       CASE (MAKE_OBS) !make observations for twin experiments
          CALL make_twin_obs()
       CASE (RUN_DIRECT) ! direct model run required 
          CALL debug('<'//TRIM(tg_ep%aa_solver_action)//'>', 'In user_post_process: nothing to do for: ')
       CASE (MAKE_ADMT) ! direct model run required 
          CALL debug('<'//TRIM(tg_ep%aa_solver_action)//'>', 'In user_post_process: nothing to do for: ')
       CASE (RUN_ADJOINT)! adjoint model run required
          CALL debug('<'//TRIM(tg_ep%aa_solver_action)//'>', 'In user_post_process: nothing to do for: ')
       CASE DEFAULT
          CALL debug('<'//TRIM(tg_ep%aa_solver_action)//'>', 'In user_post_process: Unknown action: ')
          CALL debug('', '  aborting')
          CALL parallel_finalize
          STOP
       END SELECT
       !End of computing
       
       !Checking the result for NaN and Inf values
       IF(lg_direct)THEN
          il_direct_NaNCount = NaNCount( u(:,1) )
          il_direct_InfCount = InfCount( u(:,1) )
       END IF

       IF(lg_adjoint)THEN
          il_adjoint_NaNCount = NaNCount( u(:,2) )
          il_adjoint_InfCount = InfCount( u(:,2) )
       END IF
       SELECT CASE(tg_ep%aa_solver_action)
          CASE(RUN_COST)
             ll_cost_diverged = isNaN(tg_ep%r_cost).OR.isInf(tg_ep%r_cost)
             
             !save the direct solution for the adjoint
             CALL debug('', 'In user_stat: copying the solution of the direct solution')
             !CALL debug(ala_filename, '  - ala_filename = ')
             !CALL debug(aga_direct_fileName, '  - aga_direct_fileName = ')
             !CALL dpause()
             CALL SYSTEM( 'cp -p '//TRIM(ala_filename)//' '//TRIM(aga_direct_fileName) )
             IF( (par_size>1).AND.(par_rank==0) ) THEN
                CALL SYSTEM( 'cp -p '//TRIM(ala_com_fName)//' '//TRIM(aga_direct_com_fName) )
             END IF
          CASE(RUN_GRADIENT)
             il_gradient_NaNCount = NaNCount(tg_ep%ra_grad)
             il_gradient_InfCount = InfCount(tg_ep%ra_grad)
          CASE(RUN_COSTGRAD)
             il_gradient_NaNCount = NaNCount(tg_ep%ra_grad)
             il_gradient_InfCount = InfCount(tg_ep%ra_grad)
             ll_cost_diverged = isNaN(tg_ep%r_cost).OR.isInf(tg_ep%r_cost)
          CASE DEFAULT
             !nothing to do with other cases
       END SELECT
       !Note : the variable diverged is provided by the solver
       tg_ep%l_simul_diverged = ( diverged.OR.& !solver diverged
                                  (il_direct_NaNCount > 0).OR.(il_direct_InfCount  >0).OR.& !NaN or Inf values in the solution of direct pb
                                  (il_adjoint_NaNCount> 0).OR.(il_adjoint_InfCount >0).OR.& !NaN or Inf values in the solution of adjoint pb
                                  (il_gradient_NaNCount>0).OR.(il_gradient_InfCount>0).OR.& !NaN or Inf values in the gradient
                                  ll_cost_diverged &!The cost function is NaN or Inf
                                )
       !End of checking the result for NaN and Inf values
       !Processing alternate variables and temporary files depending on the solution 
       IF(tg_ep%l_simul_diverged)THEN !
          !debuging messages
          CALL debug('', '+-+-+-+-+-+-+-+-+-+*****************+-+-+-+-+-##############+-+-+-+-+-+-+-+-+')
          CALL debug('', 'In user_stats: Simulation was impossible at this point of the control vector')
          IF(diverged)                 CALL debug('', '  - The solver diverged'                                      )
          IF(ll_cost_diverged)         CALL debug(tg_ep%r_cost, '  - The value of the cost function = '              )
          IF(il_direct_NaNCount   > 0) CALL debug('', '  - There are <NaN> values in the direct trajectory variable' )
          IF(il_direct_InfCount   > 0)THEN
             CALL debug(il_direct_InfCount, '  - There are <Inf> values in the direct trajectory variable:' )
             CALL dpause()
          END IF
          IF(il_adjoint_NaNCount  > 0) CALL debug('', '  - There are <NaN> values in the adjoint trajectory variable')
          IF(il_adjoint_InfCount  > 0) CALL debug('', '  - There are <Inf> values in the adjoint trajectory variable')
          IF(il_gradient_NaNCount > 0) CALL debug('', '  - There are <NaN> values in the gradient variable'          )
          IF(il_gradient_InfCount > 0) CALL debug('', '  - There are <Inf> values in the gradient variable'          )
          CALL debug('', '+-+-+-+-+-+-+-+-+-+*****************+-+-+-+-+-##############+-+-+-+-+-+-+-+-+')
          CALL dpause()
       ELSE!(tg_ep%l_simul_diverged)
          IF(lg_direct)THEN
             !update the candidate restart file and save the initial condition(boundary at time t=0)
             IF(tg_ep%l_first_simul)THEN !copy directly to the restart
                CALL debug('', 'In user_stat: first simul, copying the solution to the restart file')
                CALL debug(ala_filename, '  - ala_filename = ')
                CALL debug(aga_restart_filename, '  - aga_restart_filename = ')
                !CALL dpause()
                CALL SYSTEM( 'cp -p '//TRIM(ala_filename)//' '//TRIM(aga_restart_filename) )
                CALL SYSTEM( 'cp -p '//TRIM(ala_filename)//' '//TRIM(aga_cand_filename) )
                IF( (par_size>1).AND.(par_rank==0) ) THEN
                   CALL SYSTEM( 'cp -p '//TRIM(ala_com_fName)//' '//TRIM(aga_restart_com_fName) )
                   CALL SYSTEM( 'cp -p '//TRIM(ala_com_fName)//' '//TRIM(aga_cand_com_fName) )
                END IF
             ELSE!(tg_ep%l_first_simul)
                CALL debug('', 'In user_stat: copying the solution to the candidate restart file')
                CALL debug(ala_filename, '  - ala_filename = ')
                CALL debug(aga_cand_filename, '  - aga_cand_filename = ')
                CALL SYSTEM( 'cp -p '//TRIM(ala_filename)//' '//TRIM(aga_cand_filename) )
                IF( (par_size>1).AND.(par_rank==0) ) THEN
                   CALL SYSTEM( 'cp -p '//TRIM(ala_com_fName)//' '//TRIM(aga_cand_com_fName) )
                END IF
             END IF !(tg_ep%l_first_simul)

             !Extracting and saving the initial condition when running from scalar parameters
             !Think of the parallel distribution
             IF (.NOT.tg_ep%l_run_from_ctl) THEN
                IF ( (j_lev >= j_zn) )THEN !
                   CALL debug('', 'In user_stats: preparing to save the initial condition')
                   CALL debug(j_lev, 'CALL extract_ctl_idx at level ')
                   CALL extract_ctl_idx(j_lev)!extracting ctl indices
                   CALL debug( COUNT(iga_ctl_idx > 0), 'Nonzero elements in iga_ctl_idx = ')
                   CALL debug(SIZE(iga_ctl_idx), 'CALL set_ctlsize with ')
                   CALL set_ctlsize(tl_ep, SIZE(iga_ctl_idx) )
                   CALL debug('', 'assigning ctl')
                   WHERE(iga_ctl_idx > 0) tl_ep%ra_ctl = u(iga_ctl_idx, 1)
                   !Gathering informations from every processes, this approach is time consumming
                   DO ibi=1,SIZE(iga_ctl_idx)
                      CALL parallel_global_sum( REAL=tl_ep%ra_ctl(ibi) )
                   END DO
                   !u(iga_ctl_idx, 1)
                   IF( par_rank==0 ) THEN
                      CALL debug('', 'saving ctl')
                      CALL write_ep_data(tl_ep, CTL_DATA, RTIME_FILE, 'CTL from scalar parameters' )
                      tl_ep%ra_b_ctl = 0.8*tl_ep%ra_ctl
                      CALL write_ep_data(tl_ep, BCTL_DATA, RTIME_FILE, 'B_CTL from scalar parameters' )
                   END IF
                END IF
             ELSE
             END IF
          END IF !(lg_direct)
       END IF!(tg_ep%l_simul_diverged)

       !saving data for the simulator
       IF(par_rank==0) THEN   
          SELECT CASE(tg_ep%aa_solver_action)
             CASE (RUN_COSTGRAD) !cost function and its gradient  required
                CALL write_ep_data(tg_ep, COST_DATA, RTIME_FILE, tg_ep%aa_solver_action)
                CALL write_ep_data(tg_ep, GRAD_DATA, RTIME_FILE, tg_ep%aa_solver_action)
             CASE (RUN_COST) !cost function required
                CALL write_ep_data(tg_ep, COST_DATA, RTIME_FILE, tg_ep%aa_solver_action)
             CASE (RUN_GRADIENT) !gradient of the cost function required
                CALL write_ep_data(tg_ep, GRAD_DATA, RTIME_FILE, tg_ep%aa_solver_action)
             CASE DEFAULT
                !nothing to do
          END SELECT
       END IF
       !End of saving data for the simulator

       !End of processing
       CALL debug('', 'Exiting user_stats')
       !CALL FLUSH!for output to be written to file
    END IF!(startup_flag==-1)
  END SUBROUTINE user_stats

  !
  ! calculate any additional variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additiona variables are calculated and left in real space.
  !
  SUBROUTINE user_cal_force (u, n, t_local, force, drag, lift)
    !--Calculates drag and lift on obstacle using penalization formula
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: n
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (dim), INTENT (INOUT) :: force
    REAL (pr), INTENT (OUT) :: drag, lift
    REAL (pr), DIMENSION (n,dim) :: u
    drag = 0.0_PR
    lift = 0.0_PR

    !
    ! There is no obstacle in flow
    !

  END SUBROUTINE user_cal_force

  !
  ! calculate any additional variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additiona variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars( t_local, flag )
    IMPLICIT NONE
    REAL (pr), INTENT (IN) ::  t_local
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop


    ! Calculate the vorticity if we are in the main integration loop and we
    ! are saving the solution
    ! 
    ! NOTE we do not calculate vorticity in initial adaptation because
    ! derivatives are not setup yet (to save memmory)
    ! This shoild be changed eventually DG

    u(:, n_integrated+1:n_integrated+dim-1 ) = advection_velocity(j_lev, nwlt) !rg_v0*sin( rg_omega*x(:,2) )

  END SUBROUTINE user_additional_vars

  !
  ! calculate any additional scalar variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  SUBROUTINE user_scalar_vars( flag )
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop



  END SUBROUTINE user_scalar_vars

  !
  !************ Calculating Scales ***************************
  !
  ! Note the order of the components in the scl array
  ! correspond to u_tn, v_tn, w_tn, u_tn-1, v_tn-1, w_tn-1, u_tn-2, v_tn-2, w_tn-2
  ! 
  SUBROUTINE user_scales(flag, use_default, u_loc, nlocal, ne_local, l_n_var_adapt , l_n_var_adapt_index, &
       scl, scl_fltwt) !add
    USE precision
    USE pde
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: flag ! 0 during initial adaptation, then 1 during time advancement
    LOGICAL , INTENT(INOUT) :: use_default
    INTEGER, INTENT (IN) :: nlocal, ne_local
    REAL (pr), DIMENSION (1:nlocal,1:ne_local), INTENT (IN) :: u_loc
    LOGICAL , INTENT (IN) :: l_n_var_adapt(ne_local)
    INTEGER , INTENT (IN) :: l_n_var_adapt_index(1:ne_local)
    REAL (pr), DIMENSION (1:ne_local), INTENT (INOUT) :: scl
    REAL (pr) ,  INTENT(IN) ::scl_fltwt !weight for temporal filter on scl

    !
    ! Ignore the output of this routine and use default scales routine
    !
    use_default = .TRUE. 
    !
    ! NOTE: For a parallel run, synchronize scl(:) across the processors in user_scales.
    !       Use the subroutine scales of default_util.f90 as an example,
    !       subroutine parallel_global_sum of parallel.f90 may be helpful.
    !       Already synchronized: sumdA_global
    !
  END SUBROUTINE user_scales

  SUBROUTINE user_cal_cfl (use_default, u, cfl_out)
    USE precision
    USE sizes
    USE pde
    USE parallel
    IMPLICIT NONE
    LOGICAL , INTENT(INOUT) :: use_default
    REAL (pr),                                INTENT (INOUT) :: cfl_out
    REAL (pr), DIMENSION (nwlt,n_integrated), INTENT (IN)    :: u

    INTEGER                    :: i
    REAL (pr)                  :: floor
    REAL (pr), DIMENSION (dim) :: cfl
    REAL (pr), DIMENSION(dim,nwlt) :: h_arr

    use_default = .FALSE.

    floor = 1e-12_pr
    cfl_out = floor

    CALL get_all_local_h (h_arr)

    DO i = 1, nwlt
       cfl(1:dim) = ABS (u(i,1:dim)+Umn(1:dim)) * dt/h_arr(1:dim,i)
       cfl_out = MAX (cfl_out, MAXVAL(cfl))
    END DO
    CALL parallel_global_sum( REALMAXVAL=cfl_out )

  END SUBROUTINE user_cal_cfl

  !******************************************************************************************
  !************************************* SGS MODEL ROUTINES *********************************
  !******************************************************************************************

  !
  ! Intialize sgs model
  ! This routine is called once in the first
  ! iteration of the main time integration loop.
  ! weights and model filters have been setup for first loop when this routine is called.
  !
  SUBROUTINE user_init_sgs_model( )
    IMPLICIT NONE


    ! LDM: Giuliano

    ! THE sgs forcign is stored for the time step in sgs_mdl_force(1:nlocal,1:n_integrated)
    ! where nlocal should be nwlt.


    !          print *,'initializing LDM ...'       
    !          CALL sgs_mdl_force( u(:,1:n_integrated), nwlt, j_lev, .TRUE.)


  END SUBROUTINE user_init_sgs_model

  !
  ! calculate sgs model forcing term
  ! user_sgs_force is called int he beginning of each times step in time_adv_cn().
  ! THE sgs forcign is stored for the time step in sgs_mdl_force(1:nlocal,1:n_integrated)
  ! where nlocal should be nwlt.
  ! 
  ! Accesses u from field module, 
  !          j_lev from wlt_vars module,
  !
  SUBROUTINE  user_sgs_force (u_loc, nlocal)
    IMPLICIT NONE

    INTEGER,                         INTENT (IN) :: nlocal
    REAL (pr), DIMENSION (nlocal,n_integrated), INTENT (INOUT) :: u_loc

  END SUBROUTINE  user_sgs_force

  FUNCTION user_sound_speed (u, neq, nwlt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nwlt, neq
    REAL (pr), DIMENSION (nwlt,neq), INTENT(IN) :: u
    REAL (pr), DIMENSION (nwlt) :: user_sound_speed

    !user_sound_speed(:) = SQRT( gamma*(gamma-1.0_pr)* &
    !       ( u(:,4)-0.5_pr*SUM(u(:,n_var_mom(1:dim))**2,DIM=2)/u(:,1)) ) ! pressure
    user_sound_speed = 0.0_pr

  END FUNCTION user_sound_speed



  !-----------------------------------------------------------------------------------
  SUBROUTINE user_interpolate ( u )
    ! This is an example of database independent interpolation subroutine usage;
    ! predefined number of points and stack allocation has been used for that example.
    ! Read also the info inside subroutine interpolate of module wavelet_filters_mod
    USE parallel
    IMPLICIT NONE
    REAL (pr), DIMENSION (nwlt, n_var), INTENT (IN) :: u
    REAL (pr), DIMENSION (n_var, nwlt, dim) :: du, d2u
    INTEGER, PARAMETER :: x_size = 5000                  ! number of points to interpolate into
    REAL(pr) :: points(dim,x_size)
    INTEGER :: var(1:n_var)                           ! interpolate for all variables
    REAL (pr) :: res(1:n_var,1:x_size)                ! result of the interpolation
    INTEGER :: i, j, interpolation_order, var_size
    LOGICAL, SAVE :: BEEN_HERE = .FALSE.              ! file initialization flag
    INTEGER, PARAMETER :: iunit = 91                  ! ...
    INTEGER, SAVE :: set_counter = 0
    REAL(pr), PARAMETER :: MARGIN = 0.3_pr            ! margin
    REAL (4)  :: t0(0:2), t1(0:2)
    INTEGER :: cproc, ixyz(1:dim)

    BEEN_HERE = .FALSE. ! force overwriting

    ! set output file for debugging
    IF (BEEN_HERE) THEN
       ! open file for appending
       OPEN (UNIT=iunit, FILE='int'//TRIM(par_rank_str)//'.agr', FORM='formatted', STATUS='old', POSITION='append', ERR=1)
    ELSE
       ! create new file
       ! and write XMGR header
       OPEN (UNIT=iunit, FILE='int'//TRIM(par_rank_str)//'.agr', FORM='formatted', STATUS='unknown', ERR=1)
       WRITE (iunit,'("@g0 hidden false")')
       BEEN_HERE = .TRUE.
    END IF

    
    ! compute first derivatives to be passed to interpolate()
    ! for 2.5 order interpolation
    CALL c_diff_fast (u, du, d2u, j_lev, nwlt, HIGH_ORDER, 10, n_var, 1, n_var)
    DO i=1,n_var
       WRITE (*,'("   MAXVAL(ABS(u(1:nwlt,",I2,"))) = ",E17.10)') i, MAXVAL(ABS(u(1:nwlt,i)))
       WRITE (*,'("   MAXVAL(ABS(du(1:nwlt,",I2,"))) = ",E17.10)') i, MAXVAL(ABS(du(i,1:nwlt,1:dim)))
    END DO
    

    ! interpolate for all n_var variables,
    ! set var_size and var(:) respectively
    var_size = n_var
    DO i=1,n_var
       var(i) = i
    END DO


    ! predefine points to interpolate into
    ! diagonal has been used for that example
    DO i=1,x_size
       DO j=1,dim ! xx(0:nxyz(:),1:dim)
          points(j,i) = MARGIN + xx(0,j) + (xx(nxyz(j),j) - xx(0,j) - 2*MARGIN)*(i-1)/(1.0*(x_size-1))
       END DO
    END DO


    ! perform interpolation for different orders
    ! and write XMGR file ordered by the first coordinate
    DO interpolation_order = 1,1 !0,2 !0,2
       CALL CPU_TIME( t0(interpolation_order) )
       CALL interpolate( u, du, nwlt, n_var, &
            x_size, points, interpolation_order, var_size, var, &
            res, &                                    ! output interpolation result
            VERB = .TRUE.)                            ! show interpolation tree statistics (default T)
       CALL CPU_TIME( t1(interpolation_order) )
!!$       DO j=1, n_var
!!$          DO i=1,x_size
!!$             ixyz(1:dim) = INT( (points(1:dim,i)-xx(0,1:dim))/(xx(nxyz(j),j) - xx(0,j))*nxyz(1:dim) )
!!$             CALL DB_get_proc_by_coordinates (ixyz, j_lev, cproc)
!!$             IF (cproc.EQ.par_rank) &
!!$                  WRITE (iunit,'(2(E12.5,1X))') points(1,i), res(j,i)
!!$          END DO
!!$          WRITE (iunit,'("&")')
!!$       END DO
    END DO
    set_counter = set_counter + 1
    
    
    ! close XMGR file
    CLOSE (iunit)
    WRITE (*, '("   DEBUG: 0,1,3 order interpolation sets written ",I5," times")') set_counter
    DO interpolation_order = 0,2
       WRITE (*, '("CALL interpolate (USING CPU_TIME) = ", es12.5)') t1(interpolation_order) - t0(interpolation_order)
    END DO
    !CALL dpause()
    IF (set_counter.EQ.2) THEN
       CALL parallel_finalize
       STOP
    END IF
    
    
    ! error handling
    RETURN
1   PRINT *, 'ERROR while opening XMGR file'
    STOP 'in user_interpolate'
    
  END SUBROUTINE user_interpolate
  
  SUBROUTINE  user_pre_process
    IMPLICIT NONE
    
  END SUBROUTINE user_pre_process

!
  SUBROUTINE  user_post_process
    IMPLICIT NONE
  END SUBROUTINE user_post_process

!!Data assimilation routines
  !SUBROUTINE init

  !> \brief extract the gradient of the cost function with respect to the control vector
  !!
  !<
  SUBROUTINE extract_gradient()
    !CALL extract_ctl_idx(j_lev)
    
  END SUBROUTINE extract_gradient

  !> \brief Computes the indices of the control vector (initial condition)
  !! \details This subroutine comoputes the indices of the control vector component in the vecorized state variable (trajectory)
  !! It uses the global variable iga_ctl_coord to compute and store the coordinates of the control vector.
  !! These are the dim-dimensional index of the control vector (initial condition) on the non-adaptive mesh at j_face level of resolution
  !<
  SUBROUTINE extract_ctl_idx(jlev_in)
    USE wlt_trns_vars
    USE parallel

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: jlev_in
    INTEGER, DIMENSION(0:dim) :: i_p_face, i_p_xyz, i_p_face_xyz
    INTEGER, DIMENSION(dim) :: face, ixyz, ixyz_face
    INTEGER :: i, j, k, ii, j_df, j_face, iface, face_type, wlt_type, il_count, il_idxsize, il_current_size, cproc
    INTEGER :: NUMERATOR, DENOMINATOR

    j_face = min(j_zn,j_mx)
    !CALL debug( (/j_face, jlev_in, j_zn,  j_mx/), 'In extract_ctl_idx, j_face, jlev_in, j_zn,  j_mx = ')
    il_idxsize = PRODUCT(mxyz(1:dim-1)*2**(j_face-1)+1-prd(1:dim-1))
    IF(.NOT.ALLOCATED(iga_ctl_idx))THEN
       il_current_size = -1
    ELSE
       il_current_size = SIZE(iga_ctl_idx)
    END IF

    IF(il_current_size /= il_idxsize)THEN
       IF(il_current_size > 0) DEALLOCATE(iga_ctl_idx)
       ALLOCATE( iga_ctl_idx(il_idxsize) )
    END IF
    iga_ctl_idx = 0
    
    i_p_face(0) = 1
    DO i=1,dim
       i_p_face(i) = i_p_face(i-1)*3
    END DO
    i_p_xyz(0) = 1
    DO i=1,dim
       i_p_xyz(i) = i_p_xyz(i-1)*(1+nxyz(i))
    END DO
    i_p_face_xyz(0) = 1
    DO i=1,dim
       i_p_face_xyz(i) = i_p_face_xyz(i-1)*( mxyz(i)*2**(j_face-1)+1-prd(i) )
    END DO

    il_count = 0
    
    IF(jlev_in >= j_face) THEN!!!(jlev_in >= j_face)
       NUMERATOR = 1
       DENOMINATOR = 2**(jlev_in-j_face)!!!DENOMINATOR = 2**(j_lev-j_face)
    ELSE
       NUMERATOR = 2**(j_face-jlev_in)!!!NUMERATOR = 2**(j_face-j_lev)
       DENOMINATOR = 1
    END IF
    !WRITE(*, FMT="(A6,X,2A11,X,2A8,X,A6,X,A6)") 'ixyz', (/'ixyz_face1', 'ixyz_face2'/), (/'iface1', 'iface2'/), 'il_count', 'i'
    DO face_type = 0, 3**dim - 1
       face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
       IF(face(dim) == -1) THEN
          DO j = 1, MIN(jlev_in,j_face)!!!!DO j = 1, MIN(jlev_in,j_face)
             DO wlt_type = MIN(j-1,1),2**dim-1
                DO j_df = j, j_lev!!!!DO j_df = j, j_lev
                   DO k = 1, indx_DB(j_df,wlt_type,face_type,j)%length
                      i = indx_DB(j_df,wlt_type,face_type,j)%p(k)%i   ! i is 1...nwlt index on the adapted mesh
                      ii = indx_DB(j_df,wlt_type,face_type,j)%p(k)%ixyz ! 1-D index on non-adapted mesh
                      ixyz(1:dim) = INT(MOD(ii-1,i_p_xyz(1:dim))/i_p_xyz(0:dim-1)) ! dim-dimensional index on the non-adapted mesh at jlev_in level of resolution
!Innocent: the only way I found to compile this routine with db_wrk was to check for parallel mode
!there maybe a more appropriate variable to check this
#ifdef MULTIPROC
                      CALL DB_get_proc_by_coordinates (ixyz, jlev_in, cproc)!get the id of the proc on which the node is stored
                      IF (cproc.EQ.par_rank) THEN !check if the node is store on the local processor
#endif
                         ixyz_face(1:dim) = ixyz(1:dim)*NUMERATOR/DENOMINATOR ! dim-dimensional index on the non-adapted mesh at j_face level of resolution
                         iface = 1+SUM(ixyz_face(1:dim-1)*i_p_face_xyz(0:dim-2))
                         !il_count = il_count + 1
                         !WRITE(*, FMT="(I6,X,2I11,X,2I8,X,I6,X,I6)") ixyz, ixyz_face, iface, il_count, i
                         iga_ctl_idx(iface) = i + indx_DB(jlev_in,wlt_type,face_type,j)%shift
#ifdef MULTIPROC
                      END IF
#endif
                   END DO
                END DO
             END DO
          END DO
       END IF
    END DO
  END SUBROUTINE extract_ctl_idx
  
  
  !> initializes the environment for the computation of the cost function and the gradient simultaneously
  SUBROUTINE init_costgrad()
    IMPLICIT NONE
    INTEGER :: il_nobs
    !initialization

    CALL init_cost()
    il_nobs = get_obsSize(tg_obs)
  END SUBROUTINE init_costgrad

  !> initializes the environment for the computation of the cost function
  SUBROUTINE init_cost()
    IMPLICIT NONE

    CALL read_obs(tg_obs)
    CALL allocate_cost()
    IF (tg_obs%i_obs_level > J_MN)THEN
       CALL print_var('', 'In init_cost: warning; obs_level greater than J_MIN')
       CALL print_var(tg_obs%i_obs_level, '  obs_level = ')
       CALL print_var(J_MN              , '  J_MN      = ')
       CALL print_var('', ' forcing J_MIN to obs_level to ensure that observation points are always on the grid')
       J_MN = tg_obs%i_obs_level
    END IF
  END SUBROUTINE init_cost

  SUBROUTINE allocate_cost()
    IMPLICIT NONE
    INTEGER :: il_nobs
    !initialization
    il_nobs = get_obsSize(tg_obs)

    IF( ALLOCATED(rga_obsgap ) ) DEALLOCATE(rga_obsgap )
    IF( ALLOCATED(iga_obs_idx) ) DEALLOCATE(iga_obs_idx)
    ALLOCATE(&
         rga_obsgap (il_nobs),&
         iga_obs_idx(il_nobs) &
    )
  END SUBROUTINE allocate_cost

  !> initializes the environment for the computation of the gradient of the cost function
  SUBROUTINE init_grad()
    IMPLICIT NONE
    INTEGER :: il_nobs
    !initialization
    il_nobs = get_obsSize(tg_obs)

    CALL read_obsgap(tg_obs)
    IF( ALLOCATED( rga_ctlb ) ) DEALLOCATE( rga_ctlb )
    ALLOCATE( rga_ctlb( SIZE(rga_ctl) ) )
    !CALL print_os_gap(tg_obs)
    CALL allocate_cost()
  END SUBROUTINE init_grad
  
  !> zeroing adjoint variables
  SUBROUTINE adjoint_zeroing()
    IMPLICIT NONE
    tg_ep%ra_grad = 0.0_cp
    rga_ctlb      = 0.0_cp
    rg_costb = 1.0_cp
  END SUBROUTINE adjoint_zeroing
  
  !> \brief Computes the cost funtion
  !! \param[in] rda_ctl control vector
  !! \param[out] rd_cost contains the value of the cost function
  !<
  SUBROUTINE cost(rda_ctl, rd_cost)
    USE parallel
    IMPLICIT NONE
    !parameters for the cost function
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_ctl
    REAL(KIND=cp), INTENT(OUT)              :: rd_cost
    !local variables
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: rla_ctl, rla_x
    REAL(KIND=cp), DIMENSION(:,:), ALLOCATABLE :: rla_ctl2d
    REAL(KIND=cp) :: rl_ocost, rl_bcost, rl_reg_grad_cost, rl_reg_gd_cost
    REAL(KIND=cp) :: rl_beq, rl_geq, rl_gdeq
    INTEGER :: il_nctl, il_nobs, il_pos, ibi
    INTEGER, DIMENSION(dim) :: ila_nxyzctl
    !statements
    CALL debug('', 'Entering cost ++++++++++++++++++++++++++++')
    
    CALL extract_ctl_idx(j_lev)!computes iga_ctl_idx
    il_nctl = SIZE( iga_ctl_idx) !COUNT( iga_ctl_idx>0 )
    il_nobs = get_obsSize(tg_obs)
    rl_beq  = REAL( il_nobs, pr )/REAL( il_nctl, pr ) ! balance factor for the background term
    rl_geq  = rl_beq                                  ! balance factor for the gradient regularization term
    rl_gdeq = rl_beq                                 ! balance factor for the GD regularization term
    rl_ocost         = 0.0_pr
    rl_bcost         = 0.0_pr
    rl_reg_grad_cost = 0.0_pr
    rl_reg_gd_cost   = 0.0_pr
    
    CALL compute_obsgap()
    rl_ocost         = 0.5*SUM(tg_obs%ra_Rm1*tg_obs%ra_obsgap**2)
    !Gathering obs term from every processors, this approach is time consumming
    CALL parallel_global_sum( REAL=rl_ocost )

    !Computing the regularization term by the MASTER processor
    IF ( par_rank==0 ) THEN
       ALLOCATE( rla_ctl(il_nctl), rla_x(il_nctl) )
       rla_ctl = rda_ctl + tg_ep%ra_b_ctl
       IF(tg_ep%r_wb>0.0_cp)rl_bcost  = 0.5*SUM(tg_ep%ra_Bm1*rda_ctl**2)
       IF( tg_ep%r_wGD>0.0_pr ) CALL gd_regul(rla_ctl, rl_reg_gd_cost)
       SELECT CASE (dim)
       CASE(2)!1D + time
          IF( tg_ep%r_wGrad>0.0_pr )THEN 
             rl_reg_grad_cost = grad_regul(rla_ctl, CENTERED)!new version
          END IF
       CASE (3)!2D + time
          IF( tg_ep%r_wGrad>0.0_pr )THEN 
             rl_reg_grad_cost = grad_regul(rla_ctl, ila_nxyzctl(1), ila_nxyzctl(2), CENTERED, 1.0_cp, 1.0_cp)!new version
          END IF
       END SELECT
       rd_cost  = rl_ocost&
            + rl_beq*tg_ep%r_wb*rl_bcost&
            + rl_geq*tg_ep%r_wGrad*rl_reg_grad_cost&
            + rl_beq*tg_ep%r_wGD*rl_reg_gd_cost
       CALL debug('', 'In cost : <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
       WRITE(*,FMT="(A,ES13.5E2)")                          'rd_cost                     ==> ', rd_cost
       WRITE(*,FMT="(A,ES13.5E2)")                          '  = rl_ocost                ==> ', rl_ocost
       IF(tg_ep%r_wb>0.0_cp)THEN
          WRITE(*,FMT="(A,ES13.5E2,A,ES13.5E2,A,ES13.5E2)") '  + rl_beq*r_wb*rl_bcost    ==> ', rl_beq, ' *', tg_ep%r_wb   ,' *', rl_bcost
       END IF
       IF( tg_ep%r_wGrad>0.0_pr )THEN
          WRITE(*,FMT="(A,ES13.5E2,A,ES13.5E2,A,ES13.5E2)") '  + rl_geq*r_wGrad*rl_gcost ==> ', rl_geq, ' *', tg_ep%r_wGrad,' *', rl_reg_grad_cost
       END IF
       IF( tg_ep%r_wGD>0.0_pr )THEN
          WRITE(*,FMT="(A,ES13.5E2,A,ES13.5E2,A,ES13.5E2)") '  + rl_gdeq*r_wGD*rl_gd_cost==> ', rl_gdeq, ' *', tg_ep%r_wGD  ,' *', rl_reg_gd_cost
       END IF
    ELSE
       rd_cost = 0.0_pr
    END IF
    CALL parallel_global_sum( REAL=rd_cost )

    CALL debug('', 'Exiting cost ++++++++++++++++++++++++++++')
  END SUBROUTINE cost

  !> \brief Compute the index of the observations
  !> \detail This subroutine assumes that the coordinates of the observation points are already computed and stored in tg_obs%ia_icoord
  SUBROUTINE compute_obs_idx
    USE parallel
    IMPLICIT NONE
    !local variables
    INTEGER :: il_nobs, cproc, ibi
    INTEGER, DIMENSION(dim) :: ixyz
    !initialization
    il_nobs = get_obsSize(tg_obs)
    !statements

    CALL get_indices_by_coordinate(&
         il_nobs, 2**(j_lev-tg_obs%i_obs_level)*tg_obs%ia_icoord, j_lev, iga_obs_idx, 0&
    )! 0 for Extraction, 1 for initialization and -1 for cleaning
    !Zeroing the indexes that are not stored on the local processor
!Innocent: the only way I found to compile this routine with db_wrk was to check for parallel mode
!there maybe a more appropriate variable to check this
#ifdef MULTIPROC
    DO ibi=1,il_nobs
       ixyz = 2**(j_lev-tg_obs%i_obs_level)*tg_obs%ia_icoord(:, ibi)
       CALL DB_get_proc_by_coordinates (ixyz, j_lev, cproc)!get the id of the proc on which the node is stored
       IF (cproc /= par_rank) THEN !check if the node is not store on the local processor
          iga_obs_idx(ibi) = 0
       END IF
    END DO
#endif
    
  END SUBROUTINE compute_obs_idx
  
  !> \brief Computes the cost funtion
  !! \param[in] rda_ctl control vector
  !! \param[out] rd_cost contains the value of the cost function
  !<
  SUBROUTINE compute_obsgap()
    USE parallel
    IMPLICIT NONE
    !local variables
    INTEGER :: il_nobs, cproc, ibi
    INTEGER, DIMENSION(dim) :: ixyz
    !initialization
    il_nobs = get_obsSize(tg_obs)
    !statements
    CALL debug('', 'Entering compute_obsgap ++++++++++++++++++++++++++++')
    CALL compute_obs_idx()
    WHERE(iga_obs_idx>0)
       tg_obs%ra_obsgap = u(iga_obs_idx, 1) - tg_obs%ra_obs
    ELSEWHERE
       tg_obs%ra_obsgap = 0.0_pr
    END WHERE
    
    !Gathering informations from every processes, this approach is time consumming
    DO ibi=1,SIZE(iga_obs_idx)
       CALL parallel_global_sum( REAL=tg_obs%ra_obsgap(ibi) )
    END DO
    !saving by the MASTER processor
    IF( par_rank==0 ) THEN
       CALL write_obsgap(tg_obs)
    END IF
    WHERE(iga_obs_idx==0)
       tg_obs%ra_obsgap = 0.0_pr
    END WHERE
    CALL debug('', 'Exiting  compute_obsgap++++++++++++++++++++++++++++')
  END SUBROUTINE compute_obsgap
  
  !> \brief compute the gradient of the cost funtion
  !<
  SUBROUTINE costb(rda_ctl, rda_ctlb, rd_costb)
    USE parallel
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN)    :: rda_ctl
    REAL(KIND=cp),DIMENSION(:),INTENT(INOUT) :: rda_ctlb
    REAL(KIND=cp),INTENT(INOUT)              :: rd_costb
    !local variables
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: rla_ctl, rla_ctlb, rla_x
    REAL(KIND=cp) :: rl_ocost, rl_bcost, rl_reg_grad_cost, rl_reg_gd_cost
    REAL(KIND=cp) :: rl_ocostb, rl_bcostb, rl_reg_grad_costb, rl_reg_gd_costb
    REAL(KIND=cp) :: rl_beq, rl_geq, rl_gdeq
    INTEGER :: il_nctl, il_nobs, il_pos, ibi
    INTEGER, DIMENSION(dim) :: ila_nxyzctl

    !statements
    CALL debug('', 'Entering costb *************************************')
    !recomputing
    CALL extract_ctl_idx(j_lev)!computes iga_ctl_idx
    il_nctl = SIZE( iga_ctl_idx) !COUNT( iga_ctl_idx>0 )
    il_nobs = get_obsSize(tg_obs)
    rl_beq  = REAL( il_nobs, pr )/REAL( il_nctl, pr ) ! balance factor for the background term
    rl_geq  = rl_beq                                  ! balance factor for the gradient regularization term
    rl_gdeq = rl_beq
    ALLOCATE( rla_ctl(il_nctl), rla_x(il_nctl), rla_ctlb(il_nctl) )
    rla_ctl = rda_ctl + tg_ep%ra_b_ctl
    !end of recomputing
    !zeroing local adjoint variables
    rla_ctlb          = 0.0_pr
    rl_ocostb         = 0.0_pr
    rl_bcostb         = 0.0_pr
    rl_reg_grad_costb = 0.0_pr
    rl_reg_gd_costb   = 0.0_pr
    !Compute the gradient
    rl_ocostb         = rl_ocostb + rd_costb
    rl_bcostb         = rl_bcostb + rl_beq*tg_ep%r_wb*rd_costb
    rl_reg_grad_costb = rl_reg_grad_costb + rl_geq*tg_ep%r_wGrad*rd_costb
    rl_reg_gd_costb   = rl_reg_gd_costb   + rl_gdeq*tg_ep%r_wGD*rl_beq*rd_costb
    rd_costb          = 0.0_pr
    
    ila_nxyzctl(1:dim) = Mxyz(1:dim)*2**(j_zn-1)+1-prd(1:dim)
    IF ( par_rank==0 ) THEN
       SELECT CASE (dim)
       CASE(2)!1D + time
          IF( tg_ep%r_wGrad>0.0_pr )THEN 
             CALL grad_regulAdj(rla_ctl, rla_ctlb, CENTERED, rl_reg_grad_costb)!new version
             !ALL grad_regul_1db(rla_ctl, rla_ctlb, rl_reg_grad_costb)!old version
          END IF
       CASE (3)!2D+time
          IF( tg_ep%r_wGrad>0.0_pr )THEN 
             CALL grad_regulAdj(rla_ctl, rla_ctlb, ila_nxyzctl(1), ila_nxyzctl(2), CENTERED, 1.0_cp, 1.0_cp, rl_reg_grad_costb)!new version
          END IF
          CALL debug('','In costb: no regularization for 2D + time')
       END SELECT
       IF( tg_ep%r_wGD>0.0_pr ) CALL gd_regulb(rla_ctl, rla_ctlb, rl_reg_gd_costb)
       
       rl_reg_gd_costb = 0.0_pr
       IF(tg_ep%r_wb>0.0_cp)THEN 
          rda_ctlb = rda_ctlb + tg_ep%ra_Bm1*rda_ctl*rl_bcostb
          rl_bcostb= 0.0_cp
       END IF
       rda_ctlb = rda_ctlb + rla_ctlb
       rla_ctlb = 0.0_cp
    END IF

    rl_ocostb         = 0.0_pr
    rl_bcostb         = 0.0_pr
    rl_reg_grad_costb = 0.0_pr
    rl_reg_gd_costb   = 0.0_pr

!/!\ alert, the next section is to be moved to compute_grad
    !!!!!!!!!! gradient for the initial condition
    !iga_ctl_idx is already computed
!!$    WHERE(iga_ctl_idx>0)
!!$       rda_ctlb = rda_ctlb - u(iga_ctl_idx, 2)
!!$    END WHERE

!!$    !!!!!!!!!! gradient for scalar parameters
!!$    !gradient with respect to the amplitude
!!$    IF(tg_ep%l_amplitude)THEN
!!$       rda_ctlb(1) = rda_ctlb(1) &
!!$            - SUM( u(:,2)*exp( -(x(:,1)- rga_x0(1) )**2/(2.0_pr*rg_sigma**2) ), mask = x(:,2)==rga_x0(2) )
!!$    END IF
!!$    !gradient with respect to the location of the gaussian
!!$    IF(tg_ep%l_location) THEN
!!$       rda_ctlb(2) = rda_ctlb(2) &
!!$            - SUM( 2.0_pr*u(:,2)*( (x(:,1)- rga_x0(1))/(2.0_pr*rg_sigma**2))&
!!$                   *exp( -(x(:,1)- rga_x0(1) )**2/(2.0_pr*rg_sigma**2) ),&
!!$                   mask = x(:,2)==rga_x0(2)&
!!$                 )
!!$    END IF
!!$    !gradient with respect to the standard deviation of the gaussian
!!$    IF(tg_ep%l_sigma) THEN
!!$       rda_ctlb(3) = rda_ctlb(3) &
!!$            - SUM( 2.0_pr*u(:,2)*( (x(:,1)- rga_x0(1))**2/(2.0_pr*rg_sigma**3))&
!!$                   *exp( -(x(:,1)- rga_x0(1) )**2/(2.0_pr*rg_sigma**2) ),&
!!$                   mask = x(:,2)==rga_x0(2)&
!!$                 )
!!$    END IF
!!$    CALL debug(COUNT(x(:,2)==rga_x0(2)), 'costb: COUNT(x(:,2)==rga_x0(2)) = ')
    CALL debug('', 'Exiting costb *************************************')
  END SUBROUTINE costb

  !> This routine is the adjoint of a part of user_read_input (the part that compute the change of variable) and user_initial_condition
  !! It computes the gradient of the cost function by adding the term issued from costb(regularization, ...) and the adjoint at time zero
  !<
  SUBROUTINE compute_grad()
    USE parallel
    IMPLICIT NONE

    INTEGER :: ibi
    CALL extract_ctl_idx(j_lev)!computes iga_ctl_idx
    WHERE (iga_ctl_idx>0)
       rga_ctlb = - u(iga_ctl_idx, 2)
    END WHERE
    !Gathering informations from every processes, this approach is time consumming
    DO ibi=1,SIZE(iga_ctl_idx)
       CALL parallel_global_sum( REAL=rga_ctlb(ibi) )
    END DO
    
    IF( par_rank==0 ) THEN
       IF(tg_ep%l_useGD.AND.tg_ep%l_run_from_ctl) THEN
          CALL debug('', 'In compute_grad, calling gd_var_changeb')
          CALL gd_var_changeb(tg_ep%ra_ctl, tg_ep%ra_grad, rga_ctlb)
       ELSE
!!$       rga_ctl = tg_ep%ra_ctl
          tg_ep%ra_grad = tg_ep%ra_grad + rga_ctlb
          rga_ctlb = 0.0_pr
       END IF
    END IF
  END SUBROUTINE compute_grad

  !> \brief Compute gradient regularization term for a 1D vector, side differences
  !! \param[in] rda_u vector for wich the gradient regularization term is recquired
  !! \param[in], rda_x coordinates of elements of rda_u
  !<
  SUBROUTINE grad_regul_1d(rda_u, rd_reg, rda_x)
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN) :: rda_u
    REAL(KIND=cp),DIMENSION(:),OPTIONAL, INTENT(IN) :: rda_x
    REAL(KIND=cp),INTENT(OUT) :: rd_reg
    !local variables
    REAL(KIND=cp),DIMENSION( SIZE(rda_u) ) :: rla_ux, rla_x
    INTEGER :: ibi, il_nb
    IF ( PRESENT(rda_x) ) THEN
       rla_x = rda_x
    ELSE
       DO ibi = 1, SIZE(rla_x)
          rla_x(ibi) = REAL(ibi, pr)
       END DO
    END IF
    il_nb = SIZE(rda_u)
    rla_ux(il_nb) = ( rda_u(il_nb) - rda_u(il_nb - 1) )/( rla_x(il_nb) - rla_x(il_nb - 1) ) !left approximation on the right boundaries
    DO ibi = 2, il_nb-1
       rla_ux(ibi) = ( rda_u(ibi+1) - rda_u(ibi) )/( rla_x(ibi+1) - rla_x(ibi) )
    END DO
    rd_reg = 0.5*SUM(rla_ux**2)!/REAL(il_nb, pr)
  END SUBROUTINE grad_regul_1d

  SUBROUTINE grad_regul_1db(rda_u, rda_ub, rd_regb, rda_x)
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN)    :: rda_u
    REAL(KIND=cp),DIMENSION(:),INTENT(INOUT) :: rda_ub
    REAL(KIND=cp),DIMENSION(:),OPTIONAL, INTENT(IN) :: rda_x
    REAL(KIND=cp),INTENT(INOUT) :: rd_regb
    !local variables
    REAL(KIND=cp),DIMENSION( SIZE(rda_u) ) :: rla_x, rla_ux, rla_uxb
    INTEGER :: ibi, il_nb
    !zeroing local adjoint variables
    rla_uxb = 0.0_pr
    !recomputing 
    IF ( PRESENT(rda_x) ) THEN
       rla_x( 1:SIZE(rda_x) ) = rda_x
    ELSE
       DO ibi = 1, SIZE(rla_x)
          rla_x(ibi) = REAL(ibi, pr)
       END DO
    END IF
    il_nb = SIZE(rda_u)
    rla_ux(il_nb) = ( rda_u(il_nb) - rda_u(il_nb - 1) )/( rla_x(il_nb) - rla_x(il_nb - 1) ) !left approximation on the right boundaries
    DO ibi = 2, il_nb-1
       rla_ux(ibi) = ( rda_u(ibi+1) - rda_u(ibi) )/( rla_x(ibi+1) - rla_x(ibi) )
    END DO
    !end of recomputing
    rla_uxb = rla_uxb + rd_regb*rla_ux!/REAL(il_nb, pr)
    rd_regb = 0.0_pr
    
    DO ibi = il_nb-1, 1, -1       
       rda_ub(ibi+1) = rda_ub(ibi+1) + rla_uxb(ibi)/( rla_x(ibi+1) - rla_x(ibi) )
       rda_ub(ibi)   = rda_ub(ibi)   - rla_uxb(ibi)/( rla_x(ibi+1) - rla_x(ibi) )
       rla_uxb(ibi)  = 0.0_pr
    END DO
    rda_ub(il_nb)    = rda_ub(il_nb)     + rla_uxb(il_nb)/( rla_x(il_nb) - rla_x(il_nb - 1) )
    rda_ub(il_nb - 1)= rda_ub(il_nb - 1) - rla_uxb(il_nb)/( rla_x(il_nb) - rla_x(il_nb - 1) )
    rla_uxb(il_nb)   = 0.0_pr
  END SUBROUTINE grad_regul_1db

  !> \brief Compute GD regularization term
  !! \param[in] rda_u vector for wich the GD regularization term is recquired
  !! \param[in], rda_x coordinates of elements of rda_u
  !<
  SUBROUTINE gd_regul(rda_u, rd_reg)
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN) :: rda_u
    REAL(KIND=cp),INTENT(OUT) :: rd_reg
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ) :: rla_gd, rla_res
    
    CALL gd_var_change(rda_u, rla_gd)
    rla_res = 0.5*(rda_u - rla_gd)**2
    rd_reg = SUM(rla_res)
  END SUBROUTINE gd_regul
  
  !> \computes default trust function for generalized diffusion
  !! \param[in,out] rda_phi
  !! \details the default trus function is based on the obs position at time 0
  !<
  SUBROUTINE default_phi_1d(rda_phi)
    IMPLICIT NONE
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)  :: rda_phi
    !local var
    REAL(KIND=cp), DIMENSION( SIZE(rda_phi) ) :: rla_tmp, rla_tphi
    REAL(KIND=cp) :: rl_dx
    INTEGER :: il_obs_level, il_nobs_x, ibi, il_gd_niter
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ila_obs_pos
    
    rda_phi = 1.0_pr
    il_obs_level = tg_ep%i_obs_level
    il_nobs_x = Mxyz(1)*2**(il_obs_level-1) - 1
    ALLOCATE( ila_obs_pos(il_nobs_x) )
!!$    ila_obs_pos = ( 2**(j_zn-il_obs_level) )*(/(ibi, ibi=1, il_nobs_x)/)-2
!!$    rda_phi = 0.20_pr
!!$    rda_phi( ila_obs_pos ) = 0.9_pr
!!$    rla_tmp = 0.20_pr
!!$    rla_tmp( ila_obs_pos ) = 0.8_pr
!!$    rla_tphi = 0.0_cp
!!$    rla_tphi( ila_obs_pos )= 1.0_cp
!!$    rl_dx = 1.0_cp
!!$    il_gd_niter = 2**(j_zn-il_obs_level)
!!$    CALL gd_projection(rla_tmp, rda_phi, rla_tphi, rl_dx, ip_gd_niter)
!!$    CALL debug('','In default_phi ================================== ')
!!$    CALL debug(il_obs_level, 'il_obs_level = ')
!!$    CALL debug(ila_obs_pos,'In default_phi, ila_obs_pos = ')
!!$    CALL dpause('Ending default_phi')
  END SUBROUTINE default_phi_1d
  
  !> \computes default trust function for generalized diffusion
  !! \param[in,out] rda_phi
  !! \details the default trus function is based on the obs position at time 0
  !<
  SUBROUTINE default_phi_2d(rda_phi)
    IMPLICIT NONE
    REAL(KIND=cp), DIMENSION(:,:), INTENT(INOUT)  :: rda_phi
    !local var
    !REAL(KIND=cp), DIMENSION( SIZE(rda_phi,1), SIZE(rda_phi,2) ) :: rla_tmp, rla_tphi
    !REAL(KIND=cp) :: rl_dx, rl_dy
    !INTEGER :: il_obs_level, ibi, il_gd_niter
    !INTEGER, DIMENSION(dim) :: ila_nobs
    !INTEGER, ALLOCATABLE, DIMENSION(:) :: ila_obs_pos
    
    rda_phi = 1.0_pr
  END SUBROUTINE default_phi_2d

  !> \brief change of variable
  !! \param[in] rda_u original variable (old basis)
  !! \param[out] rda_gd variable in the new basis
  !! \detail this subroutine is just an interface to gd_projection
  !<
  SUBROUTINE gd_var_change(rda_u, rda_gd)
    IMPLICIT NONE
    REAL(KIND=cp), DIMENSION(:),            INTENT(IN)  :: rda_u
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ),INTENT(OUT) :: rda_gd
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ) :: rla_phi
    REAL(KIND=cp), DIMENSION( :, : ), ALLOCATABLE :: rla_u2d, rla_gd2d, rla_phi2d
    REAL(KIND=cp) :: rl_dx, rl_dy
    INTEGER, DIMENSION(dim) :: ila_nxyzctl
    
    ila_nxyzctl(1:dim) = Mxyz(1:dim)*2**(j_zn-1)+1-prd(1:dim)

    rl_dx = 1.0_cp!make it simple for now
    rl_dy = 1.0_cp!make it simple for now
    SELECT CASE (dim)
    CASE(2)!1D + time
       CALL default_phi_1d(rla_phi)!make it simple for now
       CALL gd_projection(rda_u, rda_gd, rla_phi, rl_dx, ip_gd_niter)
    CASE (3)!2D+time
       ALLOCATE( rla_u2d  ( ila_nxyzctl(1),ila_nxyzctl(2) ),&
                 rla_gd2d ( ila_nxyzctl(1),ila_nxyzctl(2) ),&
                 rla_phi2d( ila_nxyzctl(1),ila_nxyzctl(2) ) &
            )
       rla_u2d = RESHAPE( rda_u, ila_nxyzctl(1:2) )
       CALL default_phi_2d(rla_phi2d)!make it simple for now
       CALL gd_projection(rla_u2d, rla_gd2d, rla_phi2d, rl_dx, rl_dy, ip_gd_niter)
       rda_gd = RESHAPE( rla_gd2d, SHAPE(rda_gd) )
    END SELECT
  END SUBROUTINE gd_var_change

  !> \brief Compute GD regularization term
  !! \param[in] rda_u vector for wich the GD regularization term is recquired
  !<
  SUBROUTINE gd_regulb(rda_u, rda_ub, rd_regb)
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN) :: rda_u
    REAL(KIND=cp),DIMENSION(:),INTENT(INOUT) :: rda_ub
    REAL(KIND=cp),INTENT(INOUT) :: rd_regb
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ) :: rla_gd
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ) :: rla_gdb, rla_resb!adjoint variables
    !recomputing
    CALL gd_var_change(rda_u, rla_gd)
    !End of recomputing
    !zeroing local adjoint variables
    rla_gdb = 0.0_pr
    rla_resb = 0.0_pr
    !End of zeroing
    rla_resb = rla_resb + rd_regb
    rd_regb = 0.0_pr
    rda_ub = rda_ub + rda_u*rla_resb
    rla_gdb= rla_gdb-rla_gd*rla_resb
    rla_resb = 0.0_pr
    CALL gd_var_changeb(rda_u, rda_ub, rla_gdb)
    
  END SUBROUTINE gd_regulb

  !> 
  !<
  SUBROUTINE gd_var_changeb(rda_u, rda_ub, rda_gdb)
    IMPLICIT NONE
    REAL(KIND=cp),DIMENSION(:),INTENT(IN) :: rda_u
    REAL(KIND=cp),DIMENSION(:),INTENT(INOUT) :: rda_ub, rda_gdb
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_u) ) :: rla_phi
    REAL(KIND=cp), DIMENSION( :, : ), ALLOCATABLE :: rla_u2d, rla_gd2d, rla_phi2d
    REAL(KIND=cp), DIMENSION( :, : ), ALLOCATABLE :: rla_u2db, rla_gd2db
    REAL(KIND=cp) :: rl_dx, rl_dy
    INTEGER, DIMENSION(dim) :: ila_nxyzctl
    !recomputing
    ila_nxyzctl(1:dim) = Mxyz(1:dim)*2**(j_zn-1)+1-prd(1:dim)
    CALL debug('100','In gd_var_changeb --- ')
    rl_dx = 1.0_cp!make it simple for now
    rl_dy = 1.0_cp!make it simple for now
    !End of recomputing
    
    SELECT CASE (dim)
    CASE(2)!1D + time
       !recomputing
       CALL default_phi_1d(rla_phi)!make it simple for now
       !end of recomputing
       CALL gd_projectionadj(rda_u, rda_ub, rda_gdb, rla_phi, rl_dx, ip_gd_niter)
    CASE (3)!2D+time
       ALLOCATE( rla_u2d  ( ila_nxyzctl(1),ila_nxyzctl(2) ), rla_u2db  ( ila_nxyzctl(1),ila_nxyzctl(2) ),&
                 rla_gd2d ( ila_nxyzctl(1),ila_nxyzctl(2) ), rla_gd2db ( ila_nxyzctl(1),ila_nxyzctl(2) ),&
                 rla_phi2d( ila_nxyzctl(1),ila_nxyzctl(2) ) &
            )
       !zeroing local adjoint variables
       rla_u2db  = 0.0_pr
       rla_gd2db = 0.0_pr
       !end of zeroing
       !recomputing
       CALL debug('500','In gd_var_changeb --- ')
       rla_u2d = RESHAPE( rda_u, ila_nxyzctl(1:2) )
       CALL debug('600','In gd_var_changeb --- ')
       CALL default_phi_2d(rla_phi2d)!make it simple for now
       CALL debug('700','In gd_var_changeb --- ')
       !end of recomputing
       rla_gd2db = rla_gd2db + RESHAPE( rda_gdb, SHAPE(rla_gd2db) )
       rda_gdb   = 0.0_pr
       CALL debug('800','In gd_var_changeb --- ')
       CALL gd_projectionadj(rla_u2d, rla_u2db, rla_gd2db, rla_phi2d, rl_dx, rl_dy, ip_gd_niter)
       CALL debug('900','In gd_var_changeb --- ')
       rda_ub   = rda_ub + RESHAPE( rla_u2db, SHAPE(rda_ub) )
       CALL debug('950','In gd_var_changeb --- ')
       rla_u2db = 0.0_pr
    END SELECT
  END SUBROUTINE gd_var_changeb
  
  !> \brief initialize environment for twin obs
  !!
  !<
  SUBROUTINE init_twin_obs()
    USE parallel
    IMPLICIT NONE
    !local variables
    INTEGER  :: ibi, il_val, il_idx, il_nobs, il_obs_level
    CHARACTER(LEN=ip_fnl) :: ala_obsFName, ala_mtFName
    INTEGER , DIMENSION(:,:), ALLOCATABLE :: ila_idx
    INTEGER , DIMENSION(dim) ::ila_nobs,&
         ila_prodinf!local variable used for computation

    CALL debug('', 'Entering init_twin_obs *********************')
    SELECT CASE(tg_ep%aa_solver_action)
      CASE (MAKE_OBS)
        ala_obsFName = make_fileName(tg_ep, OBS_DATA, OUTPUT_FILE)
        ala_mtFName  = make_fileName(tg_ep, DMT_DATA, OUTPUT_FILE)
      CASE (MAKE_TOBS)
        ala_obsFName = make_fileName(tg_ep, TOBS_DATA, OUTPUT_FILE)
        ala_mtFName  = make_fileName(tg_ep, TDMT_DATA, OUTPUT_FILE)
    END SELECT
    !compute obs locations and dates   
    il_obs_level = tg_ep%i_obs_level
    IF( tg_ep%i_obs_level > J_MN )THEN
       CALL print_var('', '+-+-+-+-+-+-+-+-+-+*****************+-+-+-+-+-##############+-+-+-+-+-+-+-+-+')
       CALL print_var('', 'In init_twin_obs: error; obs_level greater than J_MN')
       CALL print_var(tg_ep%i_obs_level, '  obs_level = ')
       CALL print_var(J_MN             , '  J_MIN     = ')
       CALL print_var('', '  obs_level < J_MIN ensures that all observation points are on the grid')
       CALL print_var('', ' Check input and run again, do not forget to set J_MN_init')
       CALL print_var('', '+-+-+-+-+-+-+-+-+-+*****************+-+-+-+-+-##############+-+-+-+-+-+-+-+-+')
       CALL parallel_finalize
       STOP
    END IF
    ila_nobs(1:dim) = Mxyz(1:dim)*2**(il_obs_level-1) - 1! + 1 - prd(1:dim)
    il_nobs   = PRODUCT( ila_nobs )
    
    CALL debug(ila_nobs, '  ila_nobs ... = ')
    CALL debug(il_nobs , '  il_nobs  ... = ')

    CALL set_obsSize(tg_obs, il_nobs, dim, .TRUE., .FALSE.)!integer coord, no real coord
    CALL set_default_obs(tg_obs)
    CALL set_obs_fName( tg_obs, TRIM(ala_obsFName) )
    tg_obs%i_obs_level = il_obs_level
    ALLOCATE( iga_obs_idx(il_nobs) )
    !/!\ the next part (up to the end of the subrountine) must be replaced by call to build_regular_int_coord
    !CALL build_regular_int_coord(tg_obs%ia_icoord, ncoords=ila_nobs)

    !the index i of  variable ila_prodinf stores the product of the (i-1) first elts of ila_nobs, ila_obs actually give the shape of the non vectorized obs array
    ila_prodinf(1) = 1
    DO ibi = 2,dim
       ila_prodinf(ibi) = ila_prodinf(ibi-1)*ila_nobs(ibi-1)
    END DO
    DO il_idx = 1,il_nobs
       il_val = il_idx-1 !work with 0 based indexes
       DO ibi = dim,2,-1
          tg_obs%ia_icoord(ibi, il_idx) = il_val/ila_prodinf(ibi) + 1 !add 1 to move to 1 base indexes
          il_val = mod( il_val, ila_prodinf(ibi) )
       END DO
       tg_obs%ia_icoord(1, il_idx) = il_val+1; !add 1 to move to 1 base indexes
    END DO
    CALL debug('', 'Exiting init_twin_obs *********************')
  END SUBROUTINE init_twin_obs
  
  !> \brief Make twin obs for twin experiments in data assimilation
  !!
  !<
  SUBROUTINE make_twin_obs()
    USE parallel
    IMPLICIT NONE
    INTEGER :: il_nobs, il_obs_level, cproc, ibi
    INTEGER, DIMENSION(dim) :: ixyz
    !initialization
    il_nobs      = get_obsSize(tg_obs)
    il_obs_level = tg_obs%i_obs_level
    CALL debug('', 'Entering make_twin_obs *********************')
    CALL compute_obs_idx()

    WHERE(iga_obs_idx>0)
       tg_obs%ra_obs   = u(iga_obs_idx, 1)
       tg_obs%ra_sigma = 1.0_pr
    ELSEWHERE
       tg_obs%ra_sigma = 0.0_pr
    END WHERE
    !Gathering informations from every processes, this approach is time consumming
    DO ibi=1,SIZE(iga_obs_idx)
       CALL parallel_global_sum( REAL=tg_obs%ra_obs(ibi) )
       CALL parallel_global_sum( REAL=tg_obs%ra_sigma(ibi) )
    END DO
    !saving observation par the MASTER processor
    IF( par_rank==0 ) THEN
       CALL write_obs(tg_obs)
       CALL debug('', '  In make_twin_obs, obs contains :')
!!$    CALL print_os(tg_obs)
       CALL debug('', 'Exiting make_twin_obs *********************')
    END IF
  END SUBROUTINE make_twin_obs

  !> \brief Generates a zero  control vector
  !! \param[inout] td_ep exchange parameter
  !<
  SUBROUTINE make_zero_ctl(td_ep)
    USE parallel
    IMPLICIT NONE
    TYPE(exchange_param), INTENT(INOUT) :: td_ep

    CALL set_ctlsize(td_ep, PRODUCT( mxyz(1:dim-1)*2**(j_zn-1)+1-prd(1:dim-1) ) )
    td_ep%ra_ctl          = 0.0_cp
    !saving by the MASTER processor
    IF( par_rank==0 ) THEN
       CALL write_ep_data(td_ep, TCTL_DATA, OUTPUT_FILE, 'make_default_ctl' )
    END IF
  END SUBROUTINE make_zero_ctl
  
  !> \brief Generates a zero background control vector
  !! \param[inout] td_ep exchange parameter
  !<
  SUBROUTINE make_zero_bg(td_ep)
    USE parallel
    IMPLICIT NONE
    TYPE(exchange_param), INTENT(INOUT) :: td_ep

    CALL debug('', 'Entering make_zero_bg *******************************')
    CALL set_ctlsize(td_ep, PRODUCT( mxyz(1:dim-1)*2**(j_zn-1)+1-prd(1:dim-1) ) )
    td_ep%ra_ctl          = 0.0_cp
    td_ep%ra_b_ctl        = 0.0_cp
    !saving by the MASTER processor
    IF( par_rank==0 ) THEN
       CALL write_ep_data(td_ep, BCTL_DATA, OUTPUT_FILE, 'make_zero_bg' )
    END IF
    CALL debug('', 'Exiting make_zero_bg  *******************************')
  END SUBROUTINE make_zero_bg
  
  !> \brief compute a background by scaling the ctl
  !! \param[inout] td_ep exchange parameter
  !! \param[in] rd_scale scale factor
  !<
  SUBROUTINE make_bg_from_ctl(td_ep, rd_scale)
    USE parallel
    IMPLICIT NONE
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    REAL(pr), INTENT(IN) :: rd_scale

    CALL debug('', 'Entering make_bg_from_ctl *******************************')
    CALL set_ctlsize(td_ep, PRODUCT( mxyz(1:dim-1)*2**(j_zn-1)+1-prd(1:dim-1) ) )
    CALL read_ep_data(td_ep, CTL_DATA, INPUT_FILE)
    td_ep%ra_b_ctl = rd_scale*td_ep%ra_ctl
    !saving by the MASTER processor
    IF( par_rank==0 ) THEN
       CALL write_ep_data(td_ep, BCTL_DATA, OUTPUT_FILE, 'make_bg_from_ctl' )
    END IF
    CALL debug('', 'Exiting make_zero_bg  *******************************')
  END SUBROUTINE make_bg_from_ctl

  SUBROUTINE print_obs(ida_obs_coord, rda_uobs)
    USE parallel
    IMPLICIT NONE
    INTEGER , DIMENSION(:,:), INTENT(IN) :: ida_obs_coord
    REAL(pr), DIMENSION(:), INTENT(IN) :: rda_uobs
    INTEGER ::  ibj
    
    !print by the MASTER processor
    IF( par_rank==0 ) THEN
       WRITE(*,*) ' Num  Coord1  Coord2  Value'
       WRITE(*,*) '--------------------------'
       
       DO ibj = 1, SIZE(rda_uobs)
          WRITE(*, '(I5, I8, I8, E10.2)') ibj, ida_obs_coord(:, ibj), rda_uobs(ibj)
       END DO
    END IF
  END SUBROUTINE print_obs
!!!!!!!!!!!!!!!!!!!!
END MODULE user_case
