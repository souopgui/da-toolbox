!> \file com_tools.f90
!! Defines a set of tools for communication between models and Data Assimilation libraries/drivers
!<


!> Communication module
!!
!<
MODULE com_tools
  USE general_constant
  USE general_tools
  USE debug_tools
  USE conv_obs
  USE ncconv_obs
IMPLICIT NONE

  INTEGER, PARAMETER :: INT_COORD  = 301,& !< integer coordinates
                        REAL_COORD = 310,& !< real coordinates
                        BOTH_COORD = 311   !< both integer and real coordinates

  !enumeration for the identification of the program or routine that saves exchange parameters
  !todo this should be turn into enumeration type
  INTEGER, PARAMETER ::&
       SIMUL_PROGRAM  = 0,& !< value corresponding to the simulator
       SOLVER_PROGRAM = 1111  !< value corresponding to the solver

  INTEGER, PARAMETER ::&
       ip_ctl_out  = 43,& !< id of the output file for the control vector
       ip_cost_evol= 45,& !< id of the output file for the evolution of the cost function
       ip_grad_out = 47,& !< id of the output file for the gradient
       ip_grad_test= 49   !< id of the output file for the test of the gradient

!-----Constants for the action to be performed by the simulator-----
  !> \brief action value for make obs simulation
  !<
  CHARACTER(LEN=*), PARAMETER :: MAKE_OBS= "OBS"
  CHARACTER(LEN=*), PARAMETER :: MAKE_AOBS= "AOBS"!<analysed observation
  CHARACTER(LEN=*), PARAMETER :: MAKE_TOBS= "TOBS"!<true observation
  !> \brief action value for make obs simulation, given coordinates
  !<
  CHARACTER(LEN=*), PARAMETER :: MAKE_OBS_FROM_COORD= "OBS_FROM_COORD"
  !> \brief action value for user defined operation
  !! when the simulator action is set to this value, the operation is completely defined by the program writer
  !<
  CHARACTER(LEN=*), PARAMETER :: USER_DEFINED   = "USER"
  !> \brief action value for data startup
  !! when the simulator action is set to this value, the direct model is run on a larger domain and boundary conditions are saved for future experiments as well as initial condition.
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_STARTUP   = "STARTUP"
  !> \brief action value for data spinup
  !! When the simulator action is set to this value, the direct model is run from rest or predefined initial condition up to a statistical equilibrium. The end point is used as initial condition for subsequent experiments.
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_SPINUP   = "SPINUP"
  !> \brief action value for SETUP
  !! When the simulator action is set to this value, external data are read to set the true CTL, the background, the error covariance B and its inverse Binv.
  !! caution must be taken when using this option, make sure that the routine init_solver know about this option
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_SETUP   = "SETUP"
  !> \brief action value for the test of the gradient
  !<
  CHARACTER(LEN=*), PARAMETER :: GRAD_TEST   = "GRAD_TEST"
  !> \brief action value for Ensemble Kalman Filter analysis
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_ENKF   = "ENKF"
  !> \brief action value for Ensemble Kalman Smoother analysis
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_ENKS   = "ENKS"
  !> \brief action value for data assimilation
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_ASSIM   = "ASSIM"
  !> \brief action value for data assimilation, starting from the result of previous run, used the previous analysis as background
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_ASSIM2   = "ASSIM2"
  !> \brief action value for Implicit Particle Filter (IPF)
  !! The IPF consists in a data assimilation process followed by the generation of weighted samples
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_IPF   = "IPF"
  !> \brief Action value for direct simulation
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_DIRECT  = "DIRECT"
  !> \brief Action value for adjoint simulation
  !! The meaning of this variable depends on the context:
  !!   - in the solver, it means solving only the adjoint model,
  !!     the obs gap is supposed to be saved from previous computation
  !!   - in the minimization driver, it means computing the cost function (direct model) followed by the adjoint
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_ADJOINT = "ADJOINT"
  !> \brief Action value for the computation of the cost function (this is direct followed by cost)
  !! used only by the solver
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_COST = "COST"
  !> \brief action value for the computation of the gradient of the cost function (this is adjoint followed by g)
  !! used only by the solver
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_GRADIENT = "GRAD"
  !> \brief Action value for the computation of the cost function and the gradient simultaneously
  !!  used only by the solver that support this feature
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_COSTGRAD = "COST_GRAD"
  !> \brief Action value to tell the solver that m1qn3 is moving to the next iteration
  !! used by the solver to save some data (restart files ...)
  !<
  CHARACTER(LEN=*), PARAMETER :: RUN_NOTHING = "NOTHING"
  !> \brief action value for the generation of a control vector
  !! used only by the solver
  !<
  CHARACTER(LEN=*), PARAMETER :: MAKE_CTL = "CTL"
  !> \brief action value for the generation of a background control vector
  !! used only by the solver
  !<
  CHARACTER(LEN=*), PARAMETER :: MAKE_BG = "BG"
  !> \brief action value for changing variable
  !! used only by the solver. The change of variable is sometime necessary to accelerate the convergence of minimization algorithms. When such a change of variable is used, it is necessary to apply it to the analysed control variable.
  !<
  CHARACTER(LEN=*), PARAMETER :: CHANGE_VAR = "CHAVAR"
  !> \brief action value for the generation of a background control vector from pseudo inversion of observation
  !!
  !<
  CHARACTER(LEN=*), PARAMETER :: OBS_TO_BG = "OBS2BG"
  !> \brief action value for the generation of analysed model trajectory
  !! used only by the solver
  !<
  CHARACTER(LEN=*), PARAMETER :: MAKE_DMT  = "DMT" !< direct model trajectory
  CHARACTER(LEN=*), PARAMETER :: MAKE_ADMT = "ADMT"!< analysed model trajectory
  CHARACTER(LEN=*), PARAMETER :: MAKE_TDMT = "TDMT"!< true model trajectory
  !> \brief action value for the first order adjoint sensitivity
  !! used both by the solver and the simulator
  !<
  CHARACTER(LEN=*), PARAMETER :: FOA_SENSITIVITY = "FOA_SENSITIVITY"
  !> \brief action value for the second order adjoint sensitivity
  !! used both by the solver and the simulator
  !<
  CHARACTER(LEN=*), PARAMETER :: SOA_SENSITIVITY = "SOA_SENSITIVITY"
!------------------------------------------------------------------------------

  !> \brief Maximum size for the action variable
  !<
  INTEGER, PARAMETER :: IP_ACTION_LEN = 15

  !> \brief User defined type for exchanging parameters between programs
  !! Defines parameters to be exchanged between the solver and the minimization driver
  !<
  TYPE exchange_param
    CHARACTER(LEN=IP_ACTION_LEN) :: &
      aa_solver_action,&!< solver action
      aa_simul_action,& !< simul action
      prefix   !< prefix used for file name, the default value is equal to the value of simul action
    CHARACTER(LEN=ip_snl) :: aa_solver_path   = "nothing" !< solver path for extenal solver
    CHARACTER(LEN=ip_snl) :: restart_fileName = "rt_restart" !< restart file name for the solver
    CHARACTER(LEN=ip_snl) :: restart_candidate= "rt_restart_cand" !< candidate file for the restart
    CHARACTER(LEN=ip_snl) :: direct_fileName  = "rt_direct" !< direct solution file, used for the adjoint
    !solver control variables
    LOGICAL :: l_simul_diverged !< Says if the simulation has diverged or not, to be used with solvers that support this feature
    LOGICAL :: l_restart_from_previous!<restart th solver from the result of the previous iteration
    LOGICAL :: l_first_simul !< says if this is the first call to simulator
    LOGICAL :: l_run_from_ctl!< says ifthe solver should run from ctl
    LOGICAL :: l_save_pdata = .FALSE.!< says if the solver should save plotting data
    INTEGER :: i_iter = 0!< iteration ordinal number
    INTEGER :: i_nsimul_in_iter = 0!< number of simulation in the current iteration

!     !some model parameters
!     REAL(cp) :: r_mu
!     REAL(cp) :: r_sigma
!     REAL(cp) :: r_v0
!     REAL(cp) :: r_omega

     !regularization weighting parameters
     LOGICAl  :: l_useGD !< use Generalized diffusion projection?
     REAL(cp) :: r_wb    !< background weighting parameter
     REAL(cp) :: r_wGrad !< weighting parameter for gradient regularization
!
!     !informations on control parameters
!     REAL(cp):: r_sigmaB !< standard deviation of error in the background, used for synthetic background
!     LOGICAl :: l_amplitude !<
!     LOGICAl :: l_location  !<
!     LOGICAl :: l_sigma     !<

!     !observation parameters
!     INTEGER  :: i_obs_level !< observation level
!     INTEGER  :: i_nobs_x !< obs count in the space direction
!     INTEGER  :: i_nobs_t !< obs count in the time direction
!     REAL(cp) :: r_sigmaR

    !Implicit particle filter parameters
    INTEGER :: i_ipf_nparticle !< number of particles for the implicit particle filter
    INTEGER :: i_ipf_maxiter !< maximum number of iterations to compute the weight associated to a particle
    REAL(cp):: r_ipf_tol !< tolerance use to compute particle weight

!     !CS parameters
!     REAL(cp) :: r_max_coef !< max value for randomly generated coefficients
!     REAL(cp) :: r_nz_ratio !< ratio of nonzero coefficients for randomly andomly generated control vector
!     REAL(cp) :: r_mes_fact !< multiplicative factor between the number of measurements and the nomber of nonzero (nb_mes = rl_mes_fact*nb_nz)
!     REAL(cp) :: r_cs_reg   !< regularization parameter

    !control parameters
    INTEGER  :: i_nctl=0 !< size of the control vector!, PRIVATE
    INTEGER  :: i_npctl=0!< size of the Preconditionned control vector!, PRIVATE
    REAL(cp) :: r_cost   !< value of the cost function at ctl
    REAL(cp) :: r_costf  !< value of the tangent linear of the cost function at ctl (f for forward)
    REAL(cp) :: r_costb  !< adjoint of the cost function at ctl, used for convenience in the adjoint code
    INTEGER  :: i_ctl_lev !<level of ctl when using wavelet model
    LOGICAl  :: l_B_matrix=.FALSE. !< Says if B and Binv are given under the form of matrices, default if FALSE
    REAL(cp), DIMENSION(:), POINTER :: ra_dctl      => NULL()!< control vector
    REAL(cp), DIMENSION(:), POINTER :: ra_pctl      => NULL()!< preconditionned control vector
    REAL(cp), DIMENSION(:), POINTER :: ra_grad      => NULL()!< grad of the cost function at the control vector
    REAL(cp), DIMENSION(:), POINTER :: ra_b_ctl     => NULL()!< background control vector
    REAL(cp), DIMENSION(:), POINTER :: ra_sigma_ctl => NULL()!< standard deviation of errors in the b_ctl
    REAL(cp), DIMENSION(:), POINTER :: ra_Bm1       => NULL()!< inverse covariance of errors in the b_ctl
    !coordinates informations
    !> \brief number of physical dimensions of the control vector. It should always be greater than zero. If the CTL is a set of independant parameters, use 1.
    !<
    INTEGER :: i_ndim = 1
    !> \brief control variable coordinates, real coordinates in the computation domain; This is useful if the control vector is the discretization of a function (initial condition for example); The default values are integer 1..nctl
    !<
    REAL(dp), DIMENSION(:,:), POINTER :: ra_rcoord => NULL()
    !Matrix covariance and inverse
    REAL(cp), DIMENSION(:,:), POINTER :: ra_B_mat => NULL()   !< covariance matrix of errors in the b_ctl
    REAL(cp), DIMENSION(:,:), POINTER :: ra_Binv_mat => NULL()!< inverse covariance matrix of errors in the b_ctl

    !checking parameters
    REAL(cp) :: r_nothing!< use to avoid non used warning in some generic subroutines
    REAL(cp), DIMENSION(:), POINTER :: &
        ra_ctl_lb => NULL(),& !< lower bound of acceptable values for the control variable
        ra_ctl_ub => NULL()   !< upper bound of acceptable values for the control variable
    LOGICAl , DIMENSION(:), POINTER :: &
        la_check_ctl_lb => NULL(),& !< Check lower bound of acceptable values for the control variable
        la_check_ctl_ub => NULL()   !< Check upper bound of acceptable values for the control variable
  END TYPE exchange_param

  TYPE exchange_param_components
    CHARACTER(LEN=ip_cnl) ::& !, PRIVATE (removed for compatibility with old compilers)
      aa_title   = "aa_title"   ,&
      aa_history = "aa_history" ,&
      ep_bName   = "ep_bName"   ,&
      cost_bName = "cost_bName" ,&
      grad_bName = "grad_bName" ,&
      ctl_bName  = "ctl_bName"  ,&
      bctl_bName = "bctl_bName" ,&
      obs_bName  = "obs_bName"  ,&
      dmt_bName  = "dmt_bName"  ,&
      ims_bName  = "ims_bName"  ,&
      ogap_bName = "ogap_bName" ,&
      ctl_all_bName="ctl_all_bName",&
      ctl_bound_bName="ctl_bound_bName",&
      rt_dir     = "runtime_dir",& !< runtime directory, temporary runtime data are saved in this directory
      input_dir  = "input_dir" ,& !< input directory, input data are read from this directory
      output_dir = "output_dir" !< output directory, output data are saved in this directory

    CHARACTER(LEN=ip_cnl) :: &
      aa_solver_action = "aa_solver_action",&!< solver action
      aa_simul_action  = "aa_simul_action",&!< simul action
      prefix           = "prefix",&! simul action
      aa_solver_path   = "aa_solver_path",& !< solver path for extenal solver
      restart_fileName = "restart_fileName",& !< restart file name for the solver
      restart_candidate= "restart_candidate",& !< candidate file for the restart
      direct_fileName  = "direct_fileName" !< direct solution file, used for the adjoint
    !solver control variables
    CHARACTER(LEN=ip_lnl) ::&
        l_simul_diverged = "l_simul_diverged"&
      , l_restart_from_previous = "l_restart_from_previous"&
      , l_first_simul = "l_first_simul"&
      , l_run_from_ctl = "l_run_from_ctl"&
      , l_save_pdata   = "l_save_pdata"&
      , i_iter = "i_iter"&
      , i_nsimul_in_iter = "i_nsimul_in_iter"

!     !some model parameters
!     CHARACTER(LEN=ip_cnl) ::&
!       r_mu    = "r_mu",&
!       r_sigma = "r_sigma",&
!       r_v0    = "r_v0",&
!       r_omega = "r_omega"

!     !regularization weighting parameters
    CHARACTER(LEN=ip_cnl)  ::&
      l_useGD = "l_useGD",&
      r_wb    = "r_wb",&
      r_wGrad = "r_wGrad"
      !r_wGD   = "r_wGD",&
!
!     !informations on control parameters
!     CHARACTER(LEN=ip_cnl) ::&
!       r_sigmaB    = "r_sigmaB",&
!       l_amplitude = "l_amplitude",&
!       l_location = "l_location",&
!       l_sigma = "l_sigma"
!
!     !observation parameters
!     CHARACTER(LEN=ip_cnl) ::&
!       i_obs_level = "i_obs_level",&
!       i_nobs_x = "i_nobs_x",&
!       i_nobs_t = "i_nobs_t",&
!       r_sigmaR = "r_sigmaR"

    !informations on implicit particle filter
    CHARACTER(LEN=ip_cnl) ::&
      i_ipf_nparticle = "i_ipf_nparticle",&
      i_ipf_maxiter   = "i_ipf_maxiter",&
      r_ipf_tol       = "r_ipf_tol"

!     !CS parameters
!     CHARACTER(LEN=ip_cnl) ::&
!       r_max_coef = "r_max_coef",&
!       r_nz_ratio = "r_nz_ratio",&
!       r_mes_fact = "r_mes_fact",&
!       r_cs_reg = "r_cs_reg"

    !control parameters
    CHARACTER(LEN=ip_cnl) ::&
      i_nctl = "i_nctl",&
      i_npctl = "i_npctl",&
      r_cost = "r_cost",&
      i_ctl_lev = "i_ctl_lev",&
      l_B_matrix= "l_B_matrix"
    !variables related to the control vector
    !short name
    CHARACTER(LEN=ip_cnl) ::&
      ra_ctl     = "ra_ctl",&
      ra_pctl   = "ra_pctl",&
      ra_dctl    = "ra_dctl",&
      ra_grad    = "ra_grad",&
      ra_b_ctl   = "ra_b_ctl",&
      ra_sigma_ctl   = "ra_sigma_ctl",&
      ra_Bm1   = "ra_Bm1",&
      ra_B_mat   = "ra_B_mat",&
      ra_Binv_mat   = "ra_Binv_mat"
    !units
    CHARACTER(LEN=ip_cnl) ::&
      ra_ctlu    = "",&
      ra_pctlu  = "",&
      ra_dctlu   = "",&
      ra_gradu   = "N/A",&
      ra_b_ctlu  = "",&
      ra_sigma_ctlu  = "N/A",&
      ra_Bm1u  = "N/A",&
      ra_B_matu  = "N/A",&
      ra_Binv_matu  = "N/A"
    !long name
    CHARACTER(LEN=ip_lnl) ::&
      ra_ctlln   = "Control vector",&
      ra_pctlln = "Preconditionned control vector",&
      ra_dctlln  = "Increment of the control vector",&
      ra_gradln  = "Gradient of the cost function at the control vector",&
      ra_b_ctlln = "background estimation of the control vector",&
      ra_sigma_ctlln = "Standard deviation of the error in the control vector",&
      ra_Bm1ln = "Inverse standard Deviation of the error in the control vector",&
      ra_B_matln = "Matrix of error covariance associated with the background",&
      ra_Binv_matln = "Inverse matrix of error covariance associated with the background"

    !variables related to the coordinates
    !short name
    CHARACTER(LEN=ip_cnl) ::&
      i_ndim   = "i_ndim",&
      ra_rcoord = "ra_rcoord"
    !units
    CHARACTER(LEN=ip_cnl) ::&
      i_ndimu   = "N/A",&
      ra_rcoordu = "m"
    !long name
    CHARACTER(LEN=ip_lnl) ::&
      i_ndimln   = "Number of physical dimension of the problem",&
      ra_rcoordln = "Real coordinates of the control vector in the physical space"

    !checking parameters
    CHARACTER(LEN=ip_cnl) ::r_nothing = "r_nothing"
    CHARACTER(LEN=ip_cnl) ::&
      ra_ctl_lb   = "ra_ctl_lb",&
      ra_ctl_ub   = "ra_ctl_ub",&
      la_check_ctl_lb   = "la_check_ctl_lb",&
      la_check_ctl_ub   = "la_check_ctl_ub"
    CHARACTER(LEN=ip_cnl) ::&
      ra_ctl_lbu  = "N/A",&
      ra_ctl_ubu  = "N/A",&
      la_check_ctl_lbu  = "N/A",&
      la_check_ctl_ubu  = "N/A"
    CHARACTER(LEN=ip_lnl) ::&
      ra_ctl_lbln = "Lower bound of the control vector",&
      ra_ctl_ubln = "Upper bound of the control vector",&
      la_check_ctl_lbln = "Checking status of the lower bound of the control vector",&
      la_check_ctl_ubln = "Checking status of the upper bound of the control vector"
  END TYPE exchange_param_components

  !Data structure associated with the netcdf file containe exchange_param variable
  TYPE exchange_param_output
    CHARACTER(LEN=ip_snl) :: fileName, aa_title
    CHARACTER(LEN=ip_fnl) :: aa_history
      INTEGER :: ncid = -1

      !control parameters
      INTEGER :: i_nctl=0 !< size of the control vector
      INTEGER :: i_npctl=0 !< size of the changed control vector
      INTEGER :: r_costid = -1    !< value of the cost function at ctl
      INTEGER :: l_simul_divergedid = -1
      !attributes that defined the presence of optional fileds
      LOGICAL ::&
        l_B_matrix = .FALSE.!< Says if B and Binv are given under the form of matrices, default if FALSE
      !Coordinates
      INTEGER :: i_ndim = 1
      INTEGER :: ra_rcoordid
      !INTEGER :: i_ctl_lev !< level of ctl when using wavelet model
      INTEGER ::&
        ra_ctlid   = -1,&
        ra_pctlid = -1,&
        ra_dctlid  = -1,&
        ra_gradid  = -1,&
        ra_b_ctlid = -1,&
        ra_sigma_ctlid = -1,&
        ra_Bm1id = -1
      !Matrix operators
      INTEGER ::&
        ra_B_matid = -1,&
        ra_Binv_matid = -1

      !checking parameters
      !REAL(cp) :: r_nothing! use to avoid non used warning in some generic subroutines
      INTEGER :: &
        ra_ctl_lbid = -1,&
        ra_ctl_ubid = -1,&
        la_check_ctl_lbid = -1,&
        la_check_ctl_ubid = -1

    LOGICAL :: isOpened    = .FALSE.
  END TYPE exchange_param_output

  TYPE(exchange_param_components), PRIVATE, SAVE :: tm_epAtt
CONTAINS

  !
  !********************** Exchange parameter routines ************
  !

  !> \brief Initializes the netcdf file data strucontaining the exchange parameter
  !! @param [in, out] td_ncfile file data structure
  !! @param [in] ada_fileName name of the netcdf file
  !! @param [in] nctl size of the control variable
  !! @param [in] npctl size of the preconditionned control variable
  !<
  SUBROUTINE initEPOutput(td_ncfile, ada_fileName, nctl, npctl)
    type(exchange_param_output), INTENT(INOUT)    :: td_ncfile
    CHARACTER(LEN=ip_snl), INTENT(IN) :: ada_fileName
    INTEGER, INTENT(IN) :: nctl
    INTEGER, INTENT(IN) :: npctl

    td_ncfile%fileName = ada_fileName
    td_ncfile%i_nctl  = nctl
    td_ncfile%i_npctl  = npctl

  END SUBROUTINE initEPOutput

  !> \brief Closes the netcdf file containing the exchange parameter
  !! @param [in, out] td_ncfile file data structure
  !<
  SUBROUTINE ncEPclose(td_ncfile)
    type(exchange_param_output), INTENT(INOUT) :: td_ncfile

    CALL chkerr( nf90_close( td_ncfile%ncid), fname=td_ncfile%filename )
    td_ncfile%ncid = -1
    td_ncfile%isOpened = .FALSE.
  END SUBROUTINE ncEPclose

  !> \brief Opens the netcdf file containing the exchange parameter
  !! @param [in, out] td_ncfile file data structure
  !! @param [in] td_ep exchange parameter structure
  !<
  SUBROUTINE ncEPOpen(td_ncfile, td_ep)
    TYPE(exchange_param_output), INTENT(INOUT) :: td_ncfile
    TYPE(exchange_param), INTENT(INOUT), OPTIONAL :: td_ep
    !local variables
    INTEGER :: il_nDimensions, il_nVariables, il_nAttributes, il_unlimitedDimId, il_formatNum,&
      il_nctlid, il_ndimid, il_npctlid
    character(len = nf90_max_name) :: lc_name

    CALL debug(TRIM(td_ncfile%filename), 'In ncEPOpen, opening ', tag=dNETCDF)

    CALL chkerr(nf90_open(TRIM(td_ncfile%filename), NF90_NOCLOBBER, td_ncfile%ncid), fname=td_ncfile%filename)

    CALL chkerr(nf90_inquire(td_ncfile%ncid, il_nDimensions, il_nVariables, il_nAttributes, &
        il_unlimitedDimId, il_formatNum))
    CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_epAtt%i_nctl, il_nctlid ) )
    CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_nctlid, lc_name, td_ncfile%i_nctl ) )
    !newly added
    IF( nf90_inq_dimid( td_ncfile%ncid, tm_epAtt%i_npctl, il_npctlid )==NF90_NOERR )THEN
      CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_npctlid, lc_name, td_ncfile%i_npctl ) )
    ELSE
      IF( nf90_inq_dimid( td_ncfile%ncid, "i_nchanged", il_npctlid )==NF90_NOERR )THEN
        CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_npctlid, lc_name, td_ncfile%i_npctl ) )
      ELSE
        td_ncfile%i_npctl = td_ncfile%i_nctl
      END IF
    END IF

    IF( nf90_inq_dimid( td_ncfile%ncid, tm_epAtt%i_ndim, il_ndimid )/=NF90_NOERR )THEN
      CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_ndimid, lc_name, td_ncfile%i_ndim ) )
    ELSE
      td_ncfile%i_ndim = 1
    END IF

    !PRINT*, "getting attributres"
    CALL read_general_ep_att(td_ncfile)
    IF(PRESENT(td_ep)) CALL read_ep_att(td_ncfile, td_ep)!this should always be called after read_general_ep_att

    CALL debug(' In ncEPOpen, after read_ep_att', tag=dNETCDF)
    IF (PRESENT(td_ep)) td_ep%i_ndim = td_ncfile%i_ndim
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%r_cost          , td_ncfile%r_costid          ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%l_simul_diverged, td_ncfile%l_simul_divergedid) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_ctl          , td_ncfile%ra_ctlid          ) )
    !may not be present in old versions of files
    IF( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_pctl, td_ncfile%ra_pctlid) /= NF90_NOERR ) THEN
      IF( nf90_inq_varid( td_ncfile%ncid, "ra_c_ctl"      , td_ncfile%ra_pctlid) /= NF90_NOERR )&
        td_ncfile%ra_pctlid = -1
    END IF

    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_dctl         , td_ncfile%ra_dctlid         ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_grad         , td_ncfile%ra_gradid         ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_b_ctl        , td_ncfile%ra_b_ctlid        ) )
    !may not be present in old versions of files
    IF ( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_rcoord  , td_ncfile%ra_rcoordid   ) /= NF90_NOERR )&
      td_ncfile%ra_rcoordid = -1
    IF(td_ncfile%l_B_matrix)THEN
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_B_mat   , td_ncfile%ra_B_matid    ) )
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_Binv_mat, td_ncfile%ra_Binv_matid ) )
    ELSE
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_sigma_ctl    , td_ncfile%ra_sigma_ctlid    ) )
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_Bm1          , td_ncfile%ra_Bm1id          ) )
    END IF
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_ctl_lb       , td_ncfile%ra_ctl_lbid       ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%ra_ctl_ub       , td_ncfile%ra_ctl_ubid       ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%la_check_ctl_lb , td_ncfile%la_check_ctl_lbid ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_epAtt%la_check_ctl_ub , td_ncfile%la_check_ctl_ubid ) )

    td_ncfile%isOpened = .TRUE.
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncEPOpen

  !> \brief Reads the general attributes of the netcdf file
  !! @param [in, out] td_ncfile data structure associated with the exchange_param netcdf file
  !! eep that give general parameters for make_fileName is overwritten, use it cautiously
  !<
  SUBROUTINE read_general_ep_att(td_ncfile)
    type(exchange_param_output), INTENT(INOUT) :: td_ncfile
    !local variables
    INTEGER :: il_B_matrix

    CALL debug("  Entering read_ep_att", tag=dNETCDF)
    !reading attributes that are stored in td_ncfile data structure
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_title   , td_ncfile%aa_title   ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_history , td_ncfile%aa_history ) )
    !This attribute has been added later and may not be present in every file, default is .FALSE.
    IF( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_B_matrix , il_B_matrix )/= NF90_NOERR )&
      il_B_matrix = INT(.FALSE.)
    td_ncfile%l_B_matrix = int2l(il_B_matrix)
    !Reading attributes that are stored in eep data structure
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ep_bName   , eep%ep_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%cost_bName , eep%cost_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%grad_bName , eep%grad_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_bName  , eep%ctl_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%bctl_bName , eep%bctl_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%obs_bName  , eep%obs_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%dmt_bName  , eep%dmt_bName ) )
    IF( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ims_bName, eep%ims_bName )/= NF90_NOERR )&
      eep%ims_bName = "ims.nc"!the file may not contains this attribute
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ogap_bName , eep%ogap_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_all_bName  , eep%ctl_all_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_bound_bName, eep%ctl_bound_bName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%rt_dir     , eep%rt_dir ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%input_dir  , eep%input_dir ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%output_dir , eep%output_dir ) )
    CALL debug("  exiting read_ep_att ", tag=dNETCDF)
  END SUBROUTINE read_general_ep_att

  !> \brief Reads the attributes of the netcdf file and initialize (dynamic allocation) the exchange_param data structure
  !! @param [in, out] td_ncfile data structure associated with the exchange_param netcdf file
  !! @param [in] td_ep exchange parameter data structure
  !! td_ep is overwritten, and the global data structure, use it cautiously.
  !<
  SUBROUTINE read_ep_att( td_ncfile, td_ep )
    type(exchange_param_output), INTENT(INOUT) :: td_ncfile
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local variables
    INTEGER :: il_simul_diverged, il_from_previous, il_first_simul&
    , il_run_from_ctl, il_save_pdata, il_useGD
    !, il_amplitude, il_location, il_sigma

    CALL debug("  In read_ep_att, reading attributes", tag=dNETCDF)
    CALL set_B_matrix_status(td_ep, td_ncfile%l_B_matrix)
    CALL debug(130, "  In read_ep_att ", tag=dNETCDF)
    CALL set_ctlsize( td_ep, td_ncfile%i_nctl, td_ncfile%i_ndim, td_ncfile%i_npctl )

    CALL debug(200, "  In read_ep_att ", tag=dNETCDF)

    !Reading attributes that are stored in td_ep data structure
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_solver_action , td_ep%aa_solver_action ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_simul_action  , td_ep%aa_simul_action ) )
    !This variable may not be present in every file as it was introduced later
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%prefix, td_ep%prefix ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_solver_path   , td_ep%aa_solver_path ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%restart_fileName , td_ep%restart_fileName ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%restart_candidate, td_ep%restart_candidate ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%direct_fileName  , td_ep%direct_fileName ) )

    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_simul_diverged        , il_simul_diverged ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_restart_from_previous , il_from_previous ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_first_simul           , il_first_simul ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_run_from_ctl          , il_run_from_ctl ) )
    !The variable l_save_pdata does not necessarily exist in every file
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_save_pdata, il_save_pdata ))
    td_ep%l_simul_diverged        = int2l(il_simul_diverged)
    td_ep%l_restart_from_previous = int2l(il_from_previous)
    td_ep%l_first_simul  = int2l(il_first_simul)
    td_ep%l_run_from_ctl = int2l(il_run_from_ctl)
    td_ep%l_save_pdata   = int2l(il_save_pdata)
    !The variables i_iter and i_nsimul_in_iter do not necessarily exist in every file
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_iter  , td_ep%i_iter ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nsimul_in_iter  , td_ep%i_nsimul_in_iter ) )

!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_mu , td_ep%r_mu ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigma , td_ep%r_sigma ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_v0 , td_ep%r_v0 ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_omega , td_ep%r_omega ) )
    ! ! the regularization parameters where removed at some point
    ! ! there might be NetCDF error if reading those from some old files
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_useGD , il_useGD ) )
    td_ep%l_useGD = int2l(il_useGD)
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wb , td_ep%r_wb ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wGrad , td_ep%r_wGrad ) )
    !     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wGD , td_ep%r_wGD ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_amplitude , il_amplitude ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_location  , il_location ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_sigma     , il_sigma ) )
!     td_ep%l_amplitude    = int2l(il_amplitude)
!     td_ep%l_location     = int2l(il_location)
!     td_ep%l_sigma        = int2l(il_sigma)
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_obs_level, td_ep%i_obs_level) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nobs_x   , td_ep%i_nobs_x   ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nobs_t   , td_ep%i_nobs_t   ) )
!     !this attribute has been added later and may not be present in every file
!     IF( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigmaB, td_ep%r_sigmaB)/= NF90_NOERR )&
!       td_ep%r_sigmaB = REAL(1.0d10, cp)!very large value
!     !CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigmaB   , td_ep%r_sigmaB   ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigmaR   , td_ep%r_sigmaR   ) )

    !reading implicit filter parameters
    IF (td_ep%aa_solver_action==RUN_IPF) THEN
      CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ipf_nparticle, td_ep%i_ipf_nparticle ) )
      CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ipf_maxiter  , td_ep%i_ipf_maxiter   ) )
      CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_ipf_tol      , td_ep%r_ipf_tol       ) )
    END IF
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_max_coef , td_ep%r_max_coef ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_nz_ratio , td_ep%r_nz_ratio ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_mes_fact , td_ep%r_mes_fact ) )
!     CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_cs_reg   , td_ep%r_cs_reg   ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ctl_lev  , td_ep%i_ctl_lev  ) )
    CALL debug("  ... Done", tag=dNETCDF)
    !End of saving
  END SUBROUTINE read_ep_att

  SUBROUTINE ncEPCreate( td_ncfile, td_ep )
    type(exchange_param_output), INTENT(INOUT)      :: td_ncfile
    TYPE(exchange_param), INTENT(IN) :: td_ep
    !local variables
    INTEGER, DIMENSION(0) :: ila_empty
    INTEGER, DIMENSION(2) :: ila_dims2D
    INTEGER :: il_nctl, il_ndim, il_npctl

    CALL debug( TRIM(td_ncfile%filename), 'In ncEPCreate: creating ', tag=dNETCDF )
    td_ncfile%l_B_matrix = get_B_matrix_status(td_ep)
    td_ncfile%i_ndim = td_ep%i_ndim

    CALL chkerr( nf90_create( TRIM(td_ncfile%fileName), NF90_CLOBBER, td_ncfile%ncid ), fname=td_ncfile%filename )
    CALL chkerr( nf90_def_dim( td_ncfile%ncid, tm_epAtt%i_nctl, td_ncfile%i_nctl, il_nctl ) )
    CALL chkerr( nf90_def_dim( td_ncfile%ncid, tm_epAtt%i_ndim, td_ncfile%i_ndim, il_ndim ) )
    CALL chkerr( nf90_def_dim( td_ncfile%ncid, tm_epAtt%i_npctl, td_ncfile%i_npctl, il_npctl ) )
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_ctl  , NF90_DOUBLE, il_nctl, td_ncfile%ra_ctlid  ) )
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_dctl , NF90_DOUBLE, il_nctl, td_ncfile%ra_dctlid ) )
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_b_ctl, NF90_DOUBLE, il_nctl, td_ncfile%ra_b_ctlid) )
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_pctl, NF90_DOUBLE, il_npctl, td_ncfile%ra_pctlid ) )
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_grad , NF90_DOUBLE, il_npctl, td_ncfile%ra_gradid  ) )
    CALL debug(200, '... in ncEPCreate ', tag=dTRACE)
    ila_dims2D=(/il_ndim,il_nctl/)
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_rcoord, NF90_DOUBLE, ila_dims2D, td_ncfile%ra_rcoordid) )
    CALL debug(300, '... in ncEPCreate ', tag=dTRACE)
    IF(td_ncfile%l_B_matrix)THEN
      ila_dims2D=(/il_nctl,il_nctl/)
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_B_mat    , NF90_DOUBLE, ila_dims2D, td_ncfile%ra_B_matid    ) )
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_Binv_mat , NF90_DOUBLE, ila_dims2D, td_ncfile%ra_Binv_matid ) )
    ELSE
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_sigma_ctl, NF90_DOUBLE, il_nctl, td_ncfile%ra_sigma_ctlid) )
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_Bm1      , NF90_DOUBLE, il_nctl, td_ncfile%ra_Bm1id      ) )
    END IF
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_ctl_lb    , NF90_DOUBLE, il_nctl, td_ncfile%ra_ctl_lbid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%ra_ctl_ub    , NF90_DOUBLE, il_nctl, td_ncfile%ra_ctl_ubid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%la_check_ctl_lb , NF90_INT, il_nctl, td_ncfile%la_check_ctl_lbid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%la_check_ctl_ub , NF90_INT, il_nctl, td_ncfile%la_check_ctl_ubid))

    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%r_cost    , NF90_DOUBLE   , ila_empty, td_ncfile%r_costid  ))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_epAtt%l_simul_diverged, NF90_INT, ila_empty, td_ncfile%l_simul_divergedid  ))

    CALL debug(700, '... in ncEPCreate ', tag=dTRACE)
    CALL save_ep_att(td_ncfile, td_ep)
    CALL chkerr(nf90_enddef(td_ncfile%ncid))
    td_ncfile%isOpened = .TRUE.
    CALL debug('... ncEPCreate -> done', tag=dNETCDF)
  END SUBROUTINE ncEPCreate

  SUBROUTINE save_ep_att( td_ncfile, td_ep )
    !sauvegarde des attributs (param?ｿｽtres) de la trajectoire
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER :: il_simul_diverged, il_from_previous, il_first_simul&
    , il_run_from_ctl, il_save_pdata, il_B_matrix, il_useGD
    !, il_amplitude, il_location, il_sigma

    CALL debug("  Saving attributes", tag=dNETCDF)
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_title   , td_ncfile%aa_title   ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_history , td_ncfile%aa_history ) )
    il_B_matrix = INT(td_ncfile%l_B_matrix)
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_B_matrix        , il_B_matrix ) )

    !writting ...
    !variable attributes
    CALL chkerr( nf90_put_att( td_ncfile%ncid, td_ncfile%ra_ctlid, "_FillValue"  , 0.0_dp ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, td_ncfile%ra_dctlid, "_FillValue" , 0.0_dp ) )
    !Global attributes
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ep_bName   , eep%ep_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%cost_bName , eep%cost_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%grad_bName , eep%grad_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_bName  , eep%ctl_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%bctl_bName , eep%bctl_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%obs_bName  , eep%obs_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%dmt_bName  , eep%dmt_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ims_bName  , eep%ims_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ogap_bName , eep%ogap_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_all_bName  , eep%ctl_all_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%ctl_bound_bName, eep%ctl_bound_bName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%rt_dir     , eep%rt_dir ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%input_dir  , eep%input_dir ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%output_dir , eep%output_dir ) )

    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_solver_action , td_ep%aa_solver_action ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_simul_action  , td_ep%aa_simul_action ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%prefix           , td_ep%prefix          ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%aa_solver_path   , td_ep%aa_solver_path ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%restart_fileName , td_ep%restart_fileName ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%restart_candidate, td_ep%restart_candidate ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%direct_fileName  , td_ep%direct_fileName ) )

    il_simul_diverged = INT(td_ep%l_simul_diverged)
    il_from_previous  = INT(td_ep%l_restart_from_previous)
    il_first_simul    = INT(td_ep%l_first_simul)
    il_run_from_ctl   = INT(td_ep%l_run_from_ctl)
    il_save_pdata     = INT(td_ep%l_save_pdata)
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_simul_diverged        , il_simul_diverged ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_restart_from_previous , il_from_previous ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_first_simul           , il_first_simul ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_run_from_ctl          , il_run_from_ctl ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_save_pdata            , il_save_pdata   ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_iter  , td_ep%i_iter ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nsimul_in_iter, td_ep%i_nsimul_in_iter ) )

!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_mu , td_ep%r_mu ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigma , td_ep%r_sigma ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_v0 , td_ep%r_v0 ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_omega , td_ep%r_omega ) )

    il_useGD = INT(td_ep%l_useGD)
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_useGD , il_useGD ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wb , td_ep%r_wb ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wGrad , td_ep%r_wGrad ) )

!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_wGD , td_ep%r_wGD ) )
!     il_amplitude = INT(td_ep%l_amplitude)
!     il_location  = INT(td_ep%l_location)
!     il_sigma     = INT(td_ep%l_sigma)
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_amplitude , il_amplitude ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_location  , il_location ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%l_sigma     , il_sigma ) )
!     !
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_obs_level, td_ep%i_obs_level) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nobs_x   , td_ep%i_nobs_x   ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_nobs_t   , td_ep%i_nobs_t   ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigmaB   , td_ep%r_sigmaB   ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_sigmaR   , td_ep%r_sigmaR   ) )

    !reading implicit filter parameters
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ipf_nparticle, td_ep%i_ipf_nparticle ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ipf_maxiter  , td_ep%i_ipf_maxiter   ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_ipf_tol      , td_ep%r_ipf_tol       ) )

!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_max_coef , td_ep%r_max_coef ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_nz_ratio , td_ep%r_nz_ratio ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_mes_fact , td_ep%r_mes_fact ) )
!     CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%r_cs_reg   , td_ep%r_cs_reg   ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_epAtt%i_ctl_lev  , td_ep%i_ctl_lev  ) )
    CALL debug("  ...Done", tag=dNETCDF)
    !End of saving
  END SUBROUTINE save_ep_att

  !> \brief Gets the size of the preconditionned control vector
  !! @param [in] td_ep exchange parameter structure
  !<
  FUNCTION get_pctlSize(td_ep) RESULT(il_npctl)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER :: il_npctl

    il_npctl = td_ep%i_npctl
  END FUNCTION get_pctlSize

  !> \brief Gets the size of the control vector
  !! @param [in] td_ep exchange parameter structure
  !<
  FUNCTION get_ctlSize(td_ep) RESULT(il_nctl)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER :: il_nctl

    il_nctl = td_ep%i_nctl
  END FUNCTION get_ctlSize

  !> \brief (deprecated) Gets the input directory
  !! @param [in] td_ep exchange parameter structure
  !<
  FUNCTION get_input_dir(td_ep) RESULT(ala_dir)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    CHARACTER(LEN=ip_snl) :: ala_dir

    CALL nothing(td_ep%r_nothing)
    ala_dir = get_data_dir(INPUT_FILE)
  END FUNCTION get_input_dir

  !> \brief Returns the number of elements of the control vector that are out of bounds
  !! @param [in] td_ep exchange parameter
  !<
  FUNCTION get_nb_outOfBound(td_ep) RESULT(il_nb_outOfBound)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER :: il_nb_outOfBound

    il_nb_outOfBound = COUNT(&
      ( td_ep%la_check_ctl_lb.AND.(td_ep%ra_b_ctl+td_ep%ra_dctl <= td_ep%ra_ctl_lb) ).OR.&
      ( td_ep%la_check_ctl_ub.AND.(td_ep%ra_b_ctl+td_ep%ra_dctl >= td_ep%ra_ctl_ub) )&
    )
  END FUNCTION get_nb_outOfBound

  !> \brief Gets the fields l_B_matrix
  !! @param [in] td_ep exchange parameter
  !<
  FUNCTION get_B_matrix_status(td_ep)RESULT(ll_status)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    !local variables
    LOGICAL :: ll_status

    ll_status = td_ep%l_B_matrix
  END FUNCTION get_B_matrix_status

  !> \brief Sets the fields l_B_matrix
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] ld_status (logical) status of B
  !! .TRUE. is B is stored as a full matrix
  !! .FALSE. is B is a diagonal matrix and stored as a vector
  !<
  SUBROUTINE set_B_matrix_status(td_ep, ld_status)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    LOGICAL, INTENT(IN) :: ld_status
    INTEGER :: il_nctl

    CALL debug(100, "  In set_B_matrix_status ", tag=dNETCDF)
    il_nctl = td_ep%i_nctl
    CALL debug(110, "  In set_B_matrix_status ", tag=dNETCDF)
    IF(.NOT.(ld_status.EQV.td_ep%l_B_matrix))THEN
      IF(td_ep%i_nctl>0)THEN
        IF(td_ep%l_B_matrix)THEN
          !CALL debug(200, "  In set_B_matrix_status ", tag=dNETCDF)
          DEALLOCATE(td_ep%ra_B_mat, td_ep%ra_Binv_mat)
          !CALL debug(250, "  In set_B_matrix_status ", tag=dNETCDF)
          ALLOCATE( td_ep%ra_sigma_ctl(il_nctl), td_ep%ra_Bm1(il_nctl) )
        ELSE
          !CALL debug(300, "  In set_B_matrix_status ", tag=dNETCDF)
          DEALLOCATE(td_ep%ra_sigma_ctl, td_ep%ra_Bm1)
          !CALL debug(350, "  In set_B_matrix_status ", tag=dNETCDF)
          ALLOCATE( td_ep%ra_B_mat(il_nctl,il_nctl), td_ep%ra_Binv_mat(il_nctl,il_nctl) )
        END IF
      END IF
      td_ep%l_B_matrix = ld_status
    END IF
  ENDSUBROUTINE set_B_matrix_status

  !> \brief Gets the dimensionality of the control vector in the physical space; Useful when the CTL is a discretization of a function
  !! @param [in] td_ep exchange parameter
  !<
  FUNCTION get_ctl_dim(td_ep)RESULT(il_dim)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    !local variables
    INTEGER :: il_dim

    il_dim = td_ep%i_ndim
  END FUNCTION get_ctl_dim

  !> \brief Resizes the array variables, control parameters
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] nctl size of the control vector
  !! @param [in] ndim number of physical dimensions of the control vector. The default value 1 is assumed.
  !! @param [in] npctl size of the preconditionned CTL, when nctl is zero, npctl is automatically set to zero
  !! \details  this procedure set the boolean parameters to .FALSE., covariance to identity, gradien to -999 and other fields to 0. id_ndim must be greater than zero if id_nctl>0. If the control vector is a set independant parameters, id_ndim mus be 1.
  !<
  SUBROUTINE set_ctlsize( td_ep, nctl, ndim, npctl )
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    INTEGER, INTENT(IN) :: nctl
    INTEGER, INTENT(IN) :: ndim, npctl
    !local variables
    INTEGER :: il_npctl

    IF(nctl==0)THEN
      il_npctl = 0
    ELSE
      il_npctl = npctl
    END IF

    IF( ANY( (/ndim, nctl, il_npctl/)<0 ) )THEN
      CALL debug((/ndim, nctl, il_npctl/), "(/ndim, nctl, npctl/) must all be >= 0, provided:",tag=dALLWAYS)
      CALL debug("npctl is automatically set to zero if nctl is zero",tag=dALLWAYS)
      CALL stop_program( )
    END IF
    CALL debug((/ndim, nctl, il_npctl/), "In set_ctlsize; ===========, (/ndim, nctl, il_npctl/) = ",tag=dTRACE)

    !if the current size is different from the required size, reallocate
    !CALL debug(0, "In set_ctlsize ", tag=dTRACE )
    CALL set_pctlsize(td_ep, il_npctl)
    IF(td_ep%i_nctl/=nctl)THEN
      IF(td_ep%i_nctl>0)THEN
        DEALLOCATE(td_ep%ra_dctl        , td_ep%ra_b_ctl       &
                 , td_ep%ra_ctl_lb      , td_ep%ra_ctl_ub      &
                 , td_ep%la_check_ctl_lb, td_ep%la_check_ctl_ub&
                 , td_ep%ra_rcoord &
        )
        IF(td_ep%l_B_matrix)THEN
          DEALLOCATE(td_ep%ra_B_mat, td_ep%ra_Binv_mat)
        ELSE
          DEALLOCATE(td_ep%ra_sigma_ctl, td_ep%ra_Bm1)
        END IF
      END IF
      IF(nctl > 0) THEN
        ALLOCATE(td_ep%ra_dctl(nctl)  , td_ep%ra_b_ctl(nctl)       &
               , td_ep%ra_ctl_lb(nctl), td_ep%la_check_ctl_lb(nctl)&
               , td_ep%ra_ctl_ub(nctl), td_ep%la_check_ctl_ub(nctl)&
               , td_ep%ra_rcoord(ndim , nctl)&
        )
        td_ep%ra_rcoord = 0.0_cp
        IF(td_ep%l_B_matrix)THEN
          ALLOCATE( td_ep%ra_B_mat(nctl,nctl), td_ep%ra_Binv_mat(nctl,nctl) )
        ELSE
          ALLOCATE( td_ep%ra_sigma_ctl(nctl), td_ep%ra_Bm1(nctl) )
        END IF
      END IF
    END IF
    !set the new size and set the fields to the default values
    td_ep%i_nctl = nctl
    td_ep%i_ndim = ndim
    CALL set_default_ctl_param(td_ep)
  END SUBROUTINE set_ctlsize

  !> \brief Resizes the array variables related to the preconditionned CTL
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] npctl size of the change control variable
  !! \details this procedure set the gradient to -999 and the changed CTL to 0.
  !<
  SUBROUTINE set_pctlsize(td_ep, npctl)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    INTEGER, INTENT(IN) :: npctl

    IF( npctl<0 )THEN
      CALL debug( npctl, "In set_pctlsize: npctl must be >= 0, provided:",tag=dALLWAYS)
      CALL stop_program( )
    END IF
    CALL debug(npctl, "In set_pctlsize; ===========, npctl = ",tag=dTRACE)

    IF( (td_ep%i_npctl /= npctl).and.(td_ep%i_npctl>0) )THEN
        DEALLOCATE(td_ep%ra_pctl, td_ep%ra_grad)
    END IF

    IF(npctl>0)THEN
        ALLOCATE( td_ep%ra_pctl(npctl), td_ep%ra_grad(npctl) )
        td_ep%ra_pctl = 0.0_cp
        td_ep%ra_grad  = -999.0_cp
    END IF

    td_ep%i_npctl = npctl
  END SUBROUTINE set_pctlsize

  SUBROUTINE set_default_ctl_param(td_ep)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    !local variables
    INTEGER :: ibi

    IF( (td_ep%i_nctl>0) )THEN
      CALL debug(SIZE(td_ep%ra_dctl ), 'in set_default_ctl_param, SIZE(td_ep%ra_dctl ) = ', tag=dTRACE)
      td_ep%ra_dctl         = 0.0_cp
      td_ep%ra_b_ctl        = 0.0_cp
      td_ep%la_check_ctl_lb = .FALSE.
      td_ep%la_check_ctl_ub = .FALSE.
      td_ep%ra_ctl_lb       = 0.0_cp
      td_ep%ra_ctl_ub       = 0.0_cp

      IF(td_ep%l_B_matrix)THEN
        DO ibi = 1,td_ep%i_nctl
          td_ep%ra_B_mat(ibi, ibi)     = 1.0_cp
          td_ep%ra_Binv_mat(ibi, ibi)  = 1.0_cp
        END DO
      ELSE
        td_ep%ra_sigma_ctl = 1.0_cp
      END IF
    END IF
  END SUBROUTINE set_default_ctl_param

  !> \brief Read exchange parameters between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter, contains the name of the file as input
  !!     and the read params as output
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !! \details the file is supposed to be opened and the associated identifier store in the variable id_fileId;
  !<
  SUBROUTINE read_ep_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    type(exchange_param_output), INTENT(INOUT) :: td_ncfile

    CALL read_ep_att( td_ncfile, td_ep )
    CALL read_ctl_size_internal(td_ep, td_ncfile)

    !If one is going to read many blocs in the same file, it is recommended to rewind
    !reading other blocs
    CALL read_ctl_all_internal(td_ep, td_ncfile)
    CALL read_cost_internal(td_ep, td_ncfile)
    CALL read_grad_internal(td_ep, td_ncfile)
    !assigning file name informations
  END SUBROUTINE read_ep_internal

  !> \brief Write exchange parameters between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter, contains the name of the file as input
  !!  and the read params as output
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !! \details the file is supposed to be opened and the associated identifier store in the variable id_fileId;
  !<
  SUBROUTINE write_ep_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL write_ctl_all_internal(td_ep, td_ncfile)
    CALL write_cost_internal(td_ep, td_ncfile)
    CALL write_grad_internal(td_ep, td_ncfile)
  END SUBROUTINE write_ep_internal

  !> \brief Internal subroutine to read the simulation state
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !<
  SUBROUTINE read_simulStatus_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    !local variable
    INTEGER :: il_simul_diverged

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%l_simul_divergedid, il_simul_diverged ) )
    td_ep%l_simul_diverged = int2l(il_simul_diverged)
  END SUBROUTINE read_simulStatus_internal

  !> \brief Internal subroutine to read the ctl size, this subroutine allocates space for dynamic fields
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !<
  SUBROUTINE read_ctl_size_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    CALL debug(100, 'In read_ctl_size_internal', tag=dTRACE)
    CALL set_ctlsize(td_ep, td_ncfile%i_nctl, td_ncfile%i_ndim, td_ncfile%i_npctl)
  END SUBROUTINE read_ctl_size_internal

  !>\brief Internal subroutine to read the cost function between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_cost_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%r_costid, td_ep%r_cost ) )
  END SUBROUTINE read_cost_internal

  !>\brief Internal subroutine to read coordinates data
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile data structure containing the file id and the variable id
  !!\details the file is supposed to be openned; the pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_coord_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    IF( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_rcoordid, td_ep%ra_rcoord )/=NF90_NOERR )td_ep%ra_rcoord=0.0_cp
  END SUBROUTINE read_coord_internal

  !>\brief Internal subroutine to read the gradient of the cost function between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_grad_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_gradid, td_ep%ra_grad ) )
  END SUBROUTINE read_grad_internal

  !>\brief Internal subroutine to read ctl data
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !! this routine set the background and zeroe the increment
  !<
  SUBROUTINE read_ctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL read_dctl_internal(td_ep, td_ncfile)
    CALL read_bctl_internal(td_ep, td_ncfile)
    CALL read_pctl_internal(td_ep, td_ncfile)
  END SUBROUTINE read_ctl_internal

  !>\brief Internal subroutine to read ctl data (ctl increment or delta ctl)
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !! the true CTL is obtained by addind the background and the increment
  !<
  SUBROUTINE read_dctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_dctlid, td_ep%ra_dctl ) )
  END SUBROUTINE read_dctl_internal

  !>\brief Internal subroutine to read the preconditionned ctl data (preconditionned ctl increment)
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !! the true CTL is obtained by addind the background and the increment
  !<
  SUBROUTINE read_pctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    IF( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_pctlid, td_ep%ra_pctl )/=NF90_NOERR)&
      td_ep%ra_pctl = 0.0_cp
  END SUBROUTINE read_pctl_internal

!   !>\brief Internal subroutine to read the tangent linear of the preconditionned ctl data (preconditionned ctl increment)
!   !! @param [in, out] td_ep exchange parameter
!   !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
!   !!@details the file is supposed to be opened and the associated
!   !! identifier store in the variable id_fileId; The pointer fields
!   !! in td_ep are supposed to be associated and the allocatable fields allocated
!   !<
!   SUBROUTINE read_pctlf_internal(td_ep, td_ncfile)
!     TYPE(exchange_param), INTENT(INOUT)     :: td_ep
!     type(exchange_param_output), INTENT(IN) :: td_ncfile
!
!     IF( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_pctlfid, td_ep%ra_pctlf )/=NF90_NOERR)&
!       td_ep%ra_pctlf = 0.0_cp
!   END SUBROUTINE read_pctlf_internal

  !>\brief Internal subroutine to read backgroung ctl data and associated informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; the pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_bctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL set_B_matrix_status(td_ep, td_ncfile%l_B_matrix)
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_b_ctlid    , td_ep%ra_b_ctl     ) )
    IF(td_ncfile%l_B_matrix)THEN
      CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_B_matid, td_ep%ra_B_mat ) )
      CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_Binv_matid, td_ep%ra_Binv_mat ) )
    ELSE
      CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_sigma_ctlid, td_ep%ra_sigma_ctl ) )
      td_ep%ra_Bm1 = 0.0_cp
      WHERE( td_ep%ra_sigma_ctl > 0.0_cp ) td_ep%ra_Bm1 = 1.0_cp/td_ep%ra_sigma_ctl**2
    END IF
    CALL read_coord_internal(td_ep, td_ncfile)
  END SUBROUTINE read_bctl_internal

  !>\brief Internal subroutine to read ctl bounds and associated informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_ctl_bound_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    INTEGER, DIMENSION( td_ncfile%i_nctl ) :: ila_check
    INTEGER :: ibi

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_ctl_lbid, td_ep%ra_ctl_lb ) )
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%ra_ctl_ubid, td_ep%ra_ctl_ub ) )
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%la_check_ctl_lbid, ila_check ) )
    DO ibi = 1, td_ncfile%i_nctl
      td_ep%la_check_ctl_lb(ibi) = int2l( ila_check(ibi) )
    END DO
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%la_check_ctl_ubid, ila_check ) )
    DO ibi = 1, td_ncfile%i_nctl
      td_ep%la_check_ctl_ub(ibi) = int2l( ila_check(ibi) )
    END DO
  END SUBROUTINE read_ctl_bound_internal

  !>\brief Internal subroutine to read all the ctl informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE read_ctl_all_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(INOUT)     :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL read_dctl_internal(td_ep, td_ncfile)
    CALL read_bctl_internal(td_ep, td_ncfile)
    CALL read_pctl_internal(td_ep, td_ncfile)
    CALL read_ctl_bound_internal(td_ep, td_ncfile)
  END SUBROUTINE read_ctl_all_internal


  !> \brief Internal routine to write the simultaneous parameter
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  SUBROUTINE write_simulStatus_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    INTEGER :: il_simul_diverged

    il_simul_diverged = INT(td_ep%l_simul_diverged)
    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%l_simul_divergedid, il_simul_diverged ) )
  END SUBROUTINE write_simulStatus_internal

  !> \brief Internal routine to write the cost function between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_cost_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%r_costid, td_ep%r_cost ) )
  END SUBROUTINE write_cost_internal

  !> \brief Writes the coordinates
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !<
  SUBROUTINE write_coord_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_rcoordid, td_ep%ra_rcoord ) )
  END SUBROUTINE write_coord_internal

  !> \brief Writes the gradient of the cost function between the solver and the minimization driver
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !<
  SUBROUTINE write_grad_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_gradid, td_ep%ra_grad ) )
  END SUBROUTINE write_grad_internal

  !>\brief Internal subroutine to write ctl data
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !! \details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_ctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile
    !local variables
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: rla_ctl

    ALLOCATE( rla_ctl(td_ep%i_nctl) )
    rla_ctl = td_ep%ra_dctl + td_ep%ra_b_ctl
    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_ctlid, rla_ctl ) )
    DEALLOCATE(rla_ctl)
    CALL write_bctl_internal(td_ep, td_ncfile)
    CALL write_dctl_internal(td_ep, td_ncfile)
    CALL write_pctl_internal(td_ep, td_ncfile)
  END SUBROUTINE write_ctl_internal

  !>\brief Internal subroutine to write backgroung ctl data and associated informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_bctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_b_ctlid    , td_ep%ra_b_ctl     ) )
    IF(td_ncfile%l_B_matrix)THEN
      CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_B_matid   , td_ep%ra_B_mat    ) )
      CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_Binv_matid, td_ep%ra_Binv_mat ) )
    ELSE
      CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_sigma_ctlid, td_ep%ra_sigma_ctl ) )
    END IF
  END SUBROUTINE write_bctl_internal

  !>\brief Internal subroutine to write the increment of ctl data
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_dctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_dctlid, td_ep%ra_dctl ) )
  END SUBROUTINE write_dctl_internal

  !>\brief Internal subroutine to write the preconditionned increment of the  ctl data
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_pctl_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN)        :: td_ep
    type(exchange_param_output), INTENT(IN) :: td_ncfile

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_pctlid, td_ep%ra_pctl ) )
  END SUBROUTINE write_pctl_internal

!   !>\brief Internal subroutine to write the tangent linear preconditionned increment of the  ctl data
!   !! @param [in, out] td_ep exchange parameter
!   !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
!   !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
!   !<
!   SUBROUTINE write_pctlf_internal(td_ep, td_ncfile)
!     TYPE(exchange_param), INTENT(IN)        :: td_ep
!     type(exchange_param_output), INTENT(IN) :: td_ncfile
!
!     CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_pctlfid, td_ep%ra_pctlf ) )
!   END SUBROUTINE write_pctlf_internal

  !>\brief Internal subroutine to write ctl bounds and associated informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_ctl_bound_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    type(exchange_param_output), INTENT(IN)         :: td_ncfile
    INTEGER, DIMENSION( td_ncfile%i_nctl ) :: ila_check
    INTEGER :: ibi

    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_ctl_lbid, td_ep%ra_ctl_lb ) )
    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%ra_ctl_ubid, td_ep%ra_ctl_ub ) )
    DO ibi = 1, td_ncfile%i_nctl
      ila_check(ibi) = INT(td_ep%la_check_ctl_lb(ibi) )
    END DO
    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%la_check_ctl_lbid, ila_check ) )
    DO ibi = 1, td_ncfile%i_nctl
      ila_check(ibi) = INT(td_ep%la_check_ctl_ub(ibi) )
    END DO
    CALL chkerr( nf90_put_var( td_ncfile%ncid, td_ncfile%la_check_ctl_ubid, ila_check ) )
  END SUBROUTINE write_ctl_bound_internal

  !>\brief Internal subroutine to write all the ctl informations
  !! @param [in, out] td_ep exchange parameter
  !! @param [in] td_ncfile Data structure used to manipulate the exchange parameter file
  !!\details the file is supposed to be opened and the associated identifier store in the variable id_fileId; The pointer fields in td_ep are supposed to be associated and the allocatable fields allocated
  !<
  SUBROUTINE write_ctl_all_internal(td_ep, td_ncfile)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    type(exchange_param_output), INTENT(IN)         :: td_ncfile

    CALL write_dctl_internal(td_ep, td_ncfile)
    CALL write_bctl_internal(td_ep, td_ncfile)
    CALL write_pctl_internal(td_ep, td_ncfile)
    CALL write_ctl_bound_internal(td_ep, td_ncfile)
  END SUBROUTINE write_ctl_all_internal

  !> \brief Reads the ctl size from the given file and allocates space for dynamic fields
  !! @param [in] td_ep exchange parameter structure
  !! @param [in] ada_fName file name
  !<
  SUBROUTINE read_ctl_size(td_ep, ada_fName)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    type(exchange_param_output) :: tl_ncfile
    CHARACTER(LEN=ip_snl), INTENT(IN)   :: ada_fName

    CALL initEPOutput(tl_ncfile, ada_fName, 0, 0)
    CALL ncEPOpen( tl_ncfile )
    CALL read_ctl_size_internal(td_ep, tl_ncfile)
    CALL ncEPclose(tl_ncfile)
  END SUBROUTINE read_ctl_size

  !> \brief Reads exchange parameter data to file
  !! @param [in] td_ep exchange parameter structure
  !! @param [in] id_dType data type (CTL_DATA, )
  !! @param [in] id_status file status (input, or runtime)
  !! @param [in] from (optional) data type to be considered when using automatic filename
  !! \details the status information is used to choose the input directory
  !! When from is present, the file name is built using from instead of id_dType,
  !!  this allow the user to read few data from a file that contains more. for
  !!  exemple reading only the ctl from a file that contains all ep data. It is
  !!  mainly using when reading the ctl size, as it can be read from any file
  !<
  SUBROUTINE read_ep_data(td_ep, id_dType, id_status, from)
    TYPE(exchange_param), INTENT(INOUT) :: td_ep
    INTEGER, INTENT(IN) :: id_dType, id_status
    INTEGER, OPTIONAL, INTENT(IN) :: from
    !local variables
    type(exchange_param_output) :: tl_ncfile
    CHARACTER(LEN=ip_fnl) :: ala_fName
    INTEGER il_from

    IF ( PRESENT(from) ) THEN
      il_from = from
    ELSE
      SELECT CASE (id_dType)
        CASE (CTL_SIZE, COORD_DATA)
          CALL stop_program('In read_ep_data; parameter <from> is mandatory when reading ctl size or coordinates')
        CASE DEFAULT
          il_from = id_dType
      END SELECT
    END IF
    ala_fName = make_fileName(il_from, id_status)
    CALL debug(ala_fName, 'In read_ep_data, reading the file ** ', tag=dTRACE)
    CALL initEPOutput(tl_ncfile, ala_fName, 0, 0)
    !CALL debug('In read_ep_data 100 ', tag=dTRACE)
    CALL ncEPOpen( tl_ncfile )
    !CALL debug('In read_ep_data 200 ', tag=dTRACE)
    !processing simulation state when needed
    SELECT CASE(id_dType)
      CASE(EP_DATA, COST_DATA, GRAD_DATA)
        CALL read_simulStatus_internal(td_ep, tl_ncfile) !always read the status
      CASE DEFAULT
        !nothing
    END SELECT
    SELECT CASE(id_dType)
      CASE (CTL_SIZE)
        CALL read_ctl_size_internal(td_ep, tl_ncfile)
      CASE (COORD_DATA)
        CALL read_coord_internal(td_ep, tl_ncfile)
      CASE (CTL_DATA, ACTL_DATA, TCTL_DATA)
        CALL read_ctl_internal(td_ep, tl_ncfile)
      CASE (BCTL_DATA)
        CALL read_bctl_internal(td_ep, tl_ncfile)
      CASE (DCTL_DATA)
        CALL read_dctl_internal(td_ep, tl_ncfile)
      CASE (PCTL_DATA)
        CALL read_pctl_internal(td_ep, tl_ncfile)
      CASE (CTL_BOUND_DATA)
        CALL read_ctl_bound_internal(td_ep, tl_ncfile)
      CASE (CTL_ALL_DATA)
        CALL read_ctl_internal(td_ep, tl_ncfile)
      CASE (EP_DATA)
        CALL read_ep_internal(td_ep, tl_ncfile)
      CASE (COST_DATA)
        CALL read_cost_internal(td_ep, tl_ncfile)
      CASE (GRAD_DATA)
        CALL read_grad_internal(td_ep, tl_ncfile)
      CASE DEFAULT
        CALL stop_program(id_dType, 'In read_ep_data; bad data type : ')
    END SELECT
    CALL ncEPclose(tl_ncfile)
  END SUBROUTINE read_ep_data

  !> \brief Writes exchange parameter data to file
  !! @param [in] td_ep exchange parameter structure
  !! @param [in] id_dType data type (CTL_DATA, )
  !! @param [in] id_status file status (output, or runtime)
  !! @param [in] ada_text text description of the action that created the
  !!   information to be saved
  !! @param [in] as (optional) as type under which the data are actually saved
  !! @param [in] prefix (optional) prefix of the filename
  !! \details the status information is used to choose the output directory
  !! the size of the control vector is automatically saved
  !! the parameter as affect only the file name, it makes it possible to save
  !! data in a file with the name associated with as instead of id_dType. For
  !! exemple, CTL_ALL_DATA, EP_DATA can be saved as CTL or BCTL, this is usefull
  !! for restart. It is up to the user to use this parameter only if necessary.
  !<
  SUBROUTINE write_ep_data(td_ep, id_dType, id_status, ada_text, as, prefix)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER, INTENT(IN) :: id_dType, id_status
    CHARACTER(LEN=*), INTENT(IN) :: ada_text
    INTEGER, INTENT(IN), OPTIONAL :: as
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: prefix
    !local variables
    type(exchange_param_output) :: tl_ncfile
    CHARACTER(LEN=ip_fnl) :: ala_fName
    REAL(cp), DIMENSION(:), POINTER:: rla_plot_data
    INTEGER :: il_as
    LOGICAL :: ll_plot_allocated

    !initializing dynamic space
    rla_plot_data => NULL()
    ll_plot_allocated = .FALSE.
    IF( PRESENT(as) ) THEN
      il_as = as
    ELSE
      il_as = id_dType
    END IF
    ala_fName = make_fileName(il_as, id_status, prefix=prefix)
    CALL initEPOutput(tl_ncfile, ala_fName, td_ep%i_nctl, td_ep%i_npctl)
    tl_ncfile%aa_title = ""
    tl_ncfile%aa_history = make_ncHistory(ada_text)
    CALL ncEPCreate( tl_ncfile, td_ep )

    !processing simulation state when needed
    CALL write_simulStatus_internal(td_ep, tl_ncfile) !always write the status
    CALL write_coord_internal(td_ep, tl_ncfile) !Write the coordinates
    SELECT CASE(id_dType)
      CASE (CTL_DATA, ACTL_DATA, TCTL_DATA)
        CALL write_ctl_internal(td_ep, tl_ncfile)
        ALLOCATE( rla_plot_data(td_ep%i_nctl) )
        ll_plot_allocated = .TRUE.
        rla_plot_data = td_ep%ra_dctl + td_ep%ra_b_ctl
      CASE (BCTL_DATA)
        CALL write_bctl_internal(td_ep, tl_ncfile)
        rla_plot_data => td_ep%ra_b_ctl
      CASE (DCTL_DATA)
        CALL write_dctl_internal(td_ep, tl_ncfile)
        rla_plot_data => td_ep%ra_dctl
      CASE (PCTL_DATA)
        CALL write_pctl_internal(td_ep, tl_ncfile)
        rla_plot_data => td_ep%ra_pctl
      CASE (CTL_BOUND_DATA)
        CALL write_ctl_bound_internal(td_ep, tl_ncfile)
      CASE (CTL_ALL_DATA)
        CALL write_ctl_all_internal(td_ep, tl_ncfile)
      CASE (EP_DATA)
        CALL write_ep_internal(td_ep, tl_ncfile)
      CASE (COST_DATA)
        CALL write_cost_internal(td_ep, tl_ncfile)
      CASE (GRAD_DATA)
        CALL write_grad_internal(td_ep, tl_ncfile)
        rla_plot_data => td_ep%ra_grad
      CASE DEFAULT
        CALL debug(id_dType, 'In write_ep_data; bad data type or no adequate routine : ',tag=dALLWAYS)
        CALL stop_program( '  Stopping the program ...')
    END SELECT
    CALL ncEPclose(tl_ncfile)

    IF( ASSOCIATED(rla_plot_data) )THEN
      CALL write_vector_for_plot(ala_fName, rla_plot_data, td_ep%ra_rcoord)
      IF ( ll_plot_allocated ) THEN
        DEALLOCATE( rla_plot_data )
      END IF
    END IF

  END SUBROUTINE write_ep_data

  !> @brief Print the control vector and its bounds
  !! @param [in] td_ep exchange parameter
  !<
  SUBROUTINE print_ctl_bound(td_ep)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    INTEGER :: ibi

    WRITE(*,*)"print_ctl_bound: control vector and its bounds"
    WRITE(*,*)"  idx left_bound control_val right_bound"
    WRITE(*,*)"--------------------------------------------"
    DO ibi = 1, SIZE(td_ep%ra_dctl)
       WRITE(*, FMT=IFORMAT, ADVANCE='NO') ibi
       IF( td_ep%la_check_ctl_lb(ibi) )THEN
          WRITE(*, FMT=RFORMAT, ADVANCE='NO') td_ep%ra_ctl_lb(ibi)
       ELSE
          WRITE(*, FMT=AFORMAT, ADVANCE='NO') '     -      '
       END IF
       WRITE(*, FMT=RFORMAT, ADVANCE='NO') td_ep%ra_b_ctl(ibi) + td_ep%ra_dctl(ibi)
       IF( td_ep%la_check_ctl_ub(ibi) )THEN
          WRITE(*, FMT=RFORMAT) td_ep%ra_ctl_ub(ibi)
       ELSE
          WRITE(*, FMT=AFORMAT) '     -      '
       END IF
    END DO
    WRITE(*,*)"--------------------------------------------"
  END SUBROUTINE print_ctl_bound

  !> @brief print exchange parameters for diagnostics
  !! @param [in] td_ep exchange parameters
  !<
  SUBROUTINE print_ep(td_ep)
    TYPE(exchange_param), INTENT(IN) :: td_ep
    CHARACTER(LEN=80)   :: array_rformat
    array_rformat = "(A,"//TRIM( NUM2STR(td_ep%i_nctl) )//RFORMAT//")"

    CALL debug('com_tools::print_ep : exchange_param-----------------------')
    CALL debug(td_ep%aa_solver_action     , '  Action to be taken (action)    = ',tag=dALLWAYS)
    CALL debug(td_ep%aa_solver_path, '  Path to the external solver    = ',tag=dALLWAYS)
    CALL debug(eep%input_dir     , '  Input directory .............. = ',tag=dALLWAYS)
    CALL debug(eep%output_dir    , '  Output directory ............  = ',tag=dALLWAYS)
    CALL debug(eep%dmt_bName     , '  Direct model output file name. = ',tag=dALLWAYS)
    CALL debug(eep%ims_bName     , '  initial model state file name. = ',tag=dALLWAYS)
    CALL debug(eep%obs_bName     , '  Observation file name......... = ',tag=dALLWAYS)
    CALL debug(eep%ogap_bName    , '  Obsgap file Name ............. = ',tag=dALLWAYS)
    CALL debug(eep%ctl_bName     , '  CTL file Name ................ = ',tag=dALLWAYS)
    CALL debug(eep%bctl_bName    , '  background ctl file name...... = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_mu          , '  Diffusion coefficient (mu)     = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_sigma       , '  Bell shape sigma      (sigma)  = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_v0          , '  Vvelocity parameter    (v0)    = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_omega       , '  Velocity oscillation  (omega)  = ',tag=dALLWAYS)

    CALL debug(td_ep%l_useGD       , '  Use GD projection?             = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_wGD         , '  GD weight                      = ',tag=dALLWAYS)
    CALL debug(td_ep%r_wb          , '  Weighting param for the bg     = ',tag=dALLWAYS)
    CALL debug(td_ep%r_wGrad       , '  Weighting param for grad regul = ',tag=dALLWAYS)

!     CALL debug(td_ep%r_sigmaB      , '  bg standard deviation(sigmaB)  = ',tag=dALLWAYS)
!     CALL debug(td_ep%l_amplitude   , '  Bell shape amplitude (amplitu) = ',tag=dALLWAYS)
!     CALL debug(td_ep%l_location    , '  Bell shape location (location) = ',tag=dALLWAYS)
!     CALL debug(td_ep%l_sigma       , '  B. shape sigma control?        = ',tag=dALLWAYS)
    !parameters of the implicit particle filter
    CALL debug(td_ep%i_ipf_nparticle, '  Number of particles in the IPF = ',tag=dALLWAYS)
    CALL debug(td_ep%i_ipf_maxiter  , '  Max iteration for IPF weight   = ',tag=dALLWAYS)
    CALL debug(td_ep%r_ipf_tol      , '  Tolerance for IPF weight       = ',tag=dALLWAYS)
    !
!     CALL debug(td_ep%i_obs_level   , '  Observation level              = ',tag=dALLWAYS)
!     CALL debug(td_ep%i_nobs_x      , '  Obs count along x direction    = ',tag=dALLWAYS)
!     CALL debug(td_ep%i_nobs_t      , '  Obs count along time direction = ',tag=dALLWAYS)
!     CALL debug(td_ep%r_sigmaR      , '  Obs standard deviation(sigmaR) = ',tag=dALLWAYS)

    CALL debug(td_ep%i_nctl        , '  Size of the control vector     = ',tag=dALLWAYS)
    CALL debug(td_ep%r_cost        , '  Value of the cost function     = ',tag=dALLWAYS)
    CALL debug(td_ep%ra_dctl       , '  delta CTL                      = ',tag=dALLWAYS)
    CALL debug(td_ep%ra_grad       , '  Gradient of the cost fnc(grad) = ',tag=dALLWAYS)
    CALL debug(td_ep%ra_b_ctl      , '  Background of C.vector (b_ctl) = ',tag=dALLWAYS)
    CALL debug(td_ep%ra_sigma_ctl  , '  ctl STD            (sigma_ctl) = ',tag=dALLWAYS)


    CALL debug(td_ep%l_first_simul , ' First simulation?               = ',tag=dALLWAYS)
    CALL debug(td_ep%l_run_from_ctl,' Run from CTL?                   = ',tag=dALLWAYS)
    CALL debug(td_ep%l_save_pdata  ,' Save date for plotting purpose? = ',tag=dALLWAYS)
    CALL debug(td_ep%i_iter        ,' current iteration              = ',tag=dALLWAYS)
    CALL debug(td_ep%i_nsimul_in_iter,' nsimul in the current iteration= ',tag=dALLWAYS)
    CALL debug(td_ep%l_restart_from_previous , ' Restart from previous simulation?              = ',tag=dALLWAYS)
    CALL debug(td_ep%l_simul_diverged , ' Did the simulation diverged?              = ',tag=dALLWAYS)
  END SUBROUTINE print_ep

  !
  !********************************** observation !!!!!!!!!!!!!!!!!!!!!!!!!!
  !

  !> \brief read observations from file
  !! @param [in, out] td_os observation structure,
  !! contains the name of the file and the dimension (type) of the problem as input and the read params as output
  !<
  SUBROUTINE read_obs(td_os)
    TYPE(obs_structure), INTENT(INOUT) :: td_os
    type(ObsOutput) :: tl_ncfile

    CALL initObsOutput( tl_ncfile, td_os%obs_fName )
    CALL ncObsOpen( tl_ncfile )
    CALL ncObsRead( tl_ncfile, td_os, rec=1 )

    td_os%i_obs_level  = tl_ncfile%i_obs_level
    td_os%ia_M_vector  = tl_ncfile%ia_M_vector
    td_os%ia_nobs = tl_ncfile%ia_nobs
    td_os%i_nobs  = tl_ncfile%i_nobs
    td_os%i_ndim  = tl_ncfile%i_ndim
    td_os%ra_dx   = tl_ncfile%ra_dx

    td_os%ra_Rm1 = 0.0_cp
    WHERE( td_os%ra_sigma > 0.0_cp ) td_os%ra_Rm1 = 1/td_os%ra_sigma**2
    CALL ncObsClose( tl_ncfile )
    !CALL debug('In read_obs: end of reading')
  END SUBROUTINE read_obs

  !> \brief read obsgap from file
  !! @param [in, out] td_os observation structure,
  !! contains the name of the file and the dimension (type) of the problem as input
  !! and the read params as output
  !<
  SUBROUTINE read_obsgap(td_os)
    TYPE(obs_structure), INTENT(INOUT) :: td_os
    type(ObsOutput) :: tl_ncfile
    INTEGER :: il_ogap

    CALL initObsOutput( tl_ncfile, td_os%ogap_fName )
    CALL ncObsOpen( tl_ncfile )
    CALL ncObsRead( tl_ncfile, td_os, rec=1, ogap=il_ogap )
    td_os%i_obs_level = tl_ncfile%i_obs_level
    td_os%ia_M_vector = tl_ncfile%ia_M_vector
    td_os%ia_nobs = tl_ncfile%ia_nobs
    td_os%ra_dx   = tl_ncfile%ra_dx
    td_os%i_nobs  = tl_ncfile%i_nobs
    td_os%i_ndim  = tl_ncfile%i_ndim

    CALL ncObsClose( tl_ncfile )
  END SUBROUTINE read_obsgap

!   !> \brief write observations in the friendly-plot format
!   !! @param [in] td_os, observation structure
!   !<
!   SUBROUTINE write_obs_for_plot(td_os)
!     TYPE(obs_structure), INTENT(INOUT) :: td_os
!     !local variables
!     TYPE(dyn_rVector), DIMENSION(SIZE(td_os%ia_icoord,1)+1) :: tla_save
!     INTEGER :: ibi, ndim
!
!     ndim = SIZE(td_os%ia_icoord,1)
!     DO ibi = 1, ndim
!       IF(td_os%l_rcoord)THEN
!         tla_save(ibi) = td_os%ra_rcoord(ibi, :)
!       ELSE
!         tla_save(ibi) = REAL(td_os%ia_icoord(ibi, :), cp)
!       END IF
!     END DO
!     tla_save(ndim+1) = td_os%ra_obs
!
!     CALL save_trj( tla_save, TRIM( make_pFileName(td_os%obs_fName) ), ld_column = .TRUE. )
!     CALL reset_drv_array( tla_save )
!   END SUBROUTINE write_obs_for_plot

  !> \brief write observations to file
  !! @param [in, out] td_os observation structure
  !! @param [in] ada_title stringg used as title (netcdf global parameter)
  !<
  SUBROUTINE write_obs(td_os, ada_title)
    TYPE(obs_structure), INTENT(INOUT) :: td_os
    type(ObsOutput) :: tl_ncfile
    CHARACTER(LEN=*), OPTIONAL :: ada_title
    CHARACTER(LEN=ip_snl)      :: ala_title

    IF ( PRESENT(ada_title) )THEN
      ala_title = TRIM(ada_title)
    ELSE
      ala_title = 'Observation (use the second param in call to write_obs for precise title)'
    END IF

    td_os%i_nobs  = SIZE(td_os%ra_obs)!important
    CALL initObsOutput(&
        tl_ncfile, td_os%obs_fName, nobs=td_os%i_nobs, &
        ndim=td_os%i_ndim, title = ala_title, &
        hist="COM_TOOLS::write_obs"&
    )
    !tl_ncfile%aa_title   = ala_title
    !tl_ncfile%aa_history = make_ncHistory("COM_TOOLS::write_obs")
    tl_ncfile%i_obs_level = td_os%i_obs_level
    tl_ncfile%ia_M_vector = td_os%ia_M_vector
    tl_ncfile%ia_nobs = td_os%ia_nobs
    tl_ncfile%ra_dx = td_os%ra_dx
    CALL ncObsCreate(tl_ncfile)
    CALL ncObsWrite(tl_ncfile, td_os)
    CALL ncObsClose(tl_ncfile)
    !writing obs to the friendly plotting format
    !CALL write_obs_for_plot(td_os)
    IF(td_os%l_rcoord)THEN
      CALL write_vector_for_plot(td_os%obs_fName, td_os%ra_obs, td_os%ra_rcoord)
    ELSE
      CALL write_vector_for_plot(td_os%obs_fName, td_os%ra_obs, REAL(td_os%ia_icoord, cp) )
    END IF
  END SUBROUTINE write_obs

  !> \brief write observations to file
  !! @param [in, out] td_os observation structure
  !! @param [in] ada_title stringg used as title (netcdf global parameter)
  !<
  SUBROUTINE write_obsgap(td_os, ada_title)
    TYPE(obs_structure), INTENT(INOUT) :: td_os
    type(ObsOutput) :: tl_ncfile
    CHARACTER(LEN=*), OPTIONAL :: ada_title
    CHARACTER(LEN=ip_snl)      :: ala_title
    INTEGER :: il_ogap

    IF ( PRESENT(ada_title) )THEN
      ala_title = TRIM(ada_title)
    ELSE
      ala_title = "Obsgap"
    END IF
    td_os%i_nobs  = SIZE(td_os%ra_obs)!important
    CALL initObsOutput(&
        tl_ncfile, td_os%ogap_fName, nobs=td_os%i_nobs, &
        ndim=td_os%i_ndim, title=ala_title,&
        hist="COM_TOOLS::write_obsgap"&
    )
    !tl_ncfile%aa_title   = ala_title
    !tl_ncfile%aa_history = "COM_TOOLS::write_obsgap"
    tl_ncfile%i_obs_level = td_os%i_obs_level
    tl_ncfile%ia_M_vector = td_os%ia_M_vector
    tl_ncfile%ia_nobs = td_os%ia_nobs
    tl_ncfile%ra_dx = td_os%ra_dx
    CALL ncObsCreate(tl_ncfile)
    CALL ncObsWrite(tl_ncfile, td_os, ogap=il_ogap)
    CALL ncObsClose(tl_ncfile)

    IF(td_os%l_rcoord)THEN
      CALL write_vector_for_plot(td_os%ogap_fName, td_os%ra_obsgap, td_os%ra_rcoord)
    ELSE
      CALL write_vector_for_plot(td_os%ogap_fName, td_os%ra_obsgap, REAL(td_os%ia_icoord, cp) )
    END IF
  END SUBROUTINE write_obsgap

END MODULE com_tools
