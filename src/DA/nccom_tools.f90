!>file nccom_tools.f90
!!\brief netcdf module for communication tools
!<
MODULE nccom_tools
  USE general_constant
  USE general_tools
  USE debug_tools
  USE com_tools
IMPLICIT NONE

  !> \brief Attributes of exchange parameters used for netcdf
  !<
  TYPE exchange_param_attributes
    CHARACTER(LEN=ip_cnl), PRIVATE ::&
      ep_fName   = "ep_fName",&! "Exchange parameter file name"   ,&
      cost_fName = "cost_fName",&! "Cost function file name" ,&
      grad_fName = "grad_fName",&! "Gradient file name" ,&
      ctl_fName  = "ctl_fName",&! "Control vector file name"  ,&
      bctl_fName = "bctl_fName",&! "Background control vector file name",&
      obs_fName  = "obs_fName",&! "Observation file name"  ,&
      dmt_fName  = "dmt_fName",&! "Direct model trajectory file name"  ,&
      ogap_fName = "ogap_fName",&! "Obs gap file namme",&
      ctl_all_fName = "ctl_all_fName",&! "Name of the file containing all the information on the control vector",&
      ctl_bound_fName = "ctl_bound_fName",&! "Name of the file containing the bounds of the control vector",&
      rt_dir     = "rt_dir",&! "Directory where temporary runtime data are saved",&
      input_dir  = "input_dir",&! "Directory containing the input data are read from",&
      output_dir = "output_dir",&! "directory where output data are saved",&
      aa_solver_action = "aa_solver_action",&! "Action required from the solver",& 
      aa_simul_action  = "aa_simul_action"! "Action required from the simulator"
        
    CHARACTER(LEN=ip_cnl) ::&
      aa_solver_path   = "aa_solver_path",&!"Path to the extenal solver",&
      restart_fileName = "restart_fileName",&!"Restart file name for the solver",&
      restart_candidate= "restart_candidate",&!"Candidate file for the restart",&
      direct_fileName  = "direct_fileName",&!"File containing the solution of the direct model",&
      l_simul_diverged = "l_simul_diverged",&!"Says if the simulation has diverged or not",&
      l_restart_from_previous = "l_restart_from_previous",&!"Restart th solver from the previous result",&
      l_first_simul    = "l_first_simul",&!"Says if this is the first call to simulator",&
      l_run_from_ctl   = "l_run_from_ctl",&!"Says ifthe solver should run from ctl"
    !some model parameters
    CHARACTER(LEN=ip_cnl) ::&
      r_mu = "r_mu",&
      r_sigma = "r_sigma",&
      r_v0 = "r_v0",&
      r_omega = "r_omega"
    
    !regularization weighting parameters
    CHARACTER(LEN=ip_cnl) ::&
      l_useGD = "l_useGD",&!"use Generalized diffusion projection?",&
      r_wGD = "r_wGD",&!"Weighting parameter for the GD term",&
      r_wb  = "r_wb",&!"Weighting parameter for the background",&
      r_wGrad = "r_wGrad"!"Weighting parameter for gradient regularization"

   !informations on control parameters
    CHARACTER(LEN=ip_cnl) ::&
      l_amplitude = "l_amplitude",&!"Control the amplitude",&
      l_location  = "l_location",&!"Control the location",&
      l_sigma  = "l_sigma"!"Control the sigma"

    !observation parameters
    CHARACTER(LEN=ip_cnl) ::&
      i_obs_level = "i_obs_level",&!"observation level",&
      i_nobs_x = "i_nobs_x",&!"obs count in the space direction",&
      i_nobs_t = "i_nobs_t",&!"obs count in the time direction",&
      r_sigmaR = "r_sigmaR"!"Standard deviation of the observation error"
    
    !CS parameters
    CHARACTER(LEN=ip_cnl) ::&
      r_max_coef = "r_max_coef",&!"Max value for randomly generated coefficients",&
      r_nz_ratio = "r_nz_ratio",&!"ratio of nonzero coefficients",&
      r_mes_fact = "r_mes_fact",&!"multiplicative factor between measurements count and the nonzero count",&
      r_cs_reg = "r_cs_reg"!"regularization parameter fo CS"

    !control parameters
    CHARACTER(LEN=ip_cnl), PRIVATE  :: i_nctl="i_nctl"!"Size of the control vector"
    CHARACTER(LEN=ip_cnl) :: r_cost = "r_cost"!"value of the cost function at ctl"
    CHARACTER(LEN=ip_cnl) :: i_ctl_lev = "i_ctl_lev"!"level of ctl when using wavelet model"
     
    CHARACTER(LEN=ip_cnl) :: ra_ctl   = "ra_ctl",&
    CHARACTER(LEN=ip_cnl) :: ra_ctlu  = "",&
    CHARACTER(LEN=ip_lnl) :: ra_ctlln = "Control vector"
    CHARACTER(LEN=ip_cnl) :: ra_grad   = "ra_grad",&
    CHARACTER(LEN=ip_cnl) :: ra_gradu  = "",&
    CHARACTER(LEN=ip_lnl) :: ra_gradln = "Gradient of the cost function at the control vector" 
    CHARACTER(LEN=ip_cnl) :: ra_b_ctl   = "ra_b_ctl",&
    CHARACTER(LEN=ip_cnl) :: ra_b_ctlu  = "N/A",&
    CHARACTER(LEN=ip_lnl) :: ra_b_ctlln = "background approximation of the control vector"!> 
    CHARACTER(LEN=ip_cnl) :: ra_sigma_ctl   = "ra_sigma_ctl",&
    CHARACTER(LEN=ip_cnl) :: ra_sigma_ctlu  = "N/A",&
    CHARACTER(LEN=ip_lnl) :: ra_sigma_ctlln = "standard deviation of errors in the b_ctl"!> 
    CHARACTER(LEN=ip_cnl) :: ra_Bm1       = "ra_Bm1",&
    CHARACTER(LEN=ip_cnl) :: ra_Bm1u      = "N/A",&
    CHARACTER(LEN=ip_lnl) :: ra_Bm1ln     = "(inverse matrix) covariance of errors in the b_ctl" 

    !checking parameters
    CHARACTER(LEN=ip_cnl) :: r_nothing! use to avoid non used warning in some generic subroutines
    CHARACTER(LEN=ip_cnl) :: ra_ctl_lb   = "ra_ctl_lb"
    CHARACTER(LEN=ip_cnl) :: ra_ctl_lbu  = "N/A"
    CHARACTER(LEN=ip_lnl) :: ra_ctl_lbln = "lower bound of acceptable values for the control variable"
    CHARACTER(LEN=ip_cnl) :: ra_ctl_ub   ="ra_ctl_ub"
    CHARACTER(LEN=ip_cnl) :: ra_ctl_ubu  ="N/A"
    CHARACTER(LEN=ip_lnl) :: ra_ctl_ubln ="upper bound of acceptable values for the control variable"
    CHARACTER(LEN=ip_cnl) :: la_check_ctl_lb   = "la_check_ctl_lb"
    CHARACTER(LEN=ip_cnl) :: la_check_ctl_lbu  = "N/A"
    CHARACTER(LEN=ip_lnl) :: la_check_ctl_lbln = "Check lower bound of acceptable values for the control variable"
    CHARACTER(LEN=ip_cnl) :: la_check_ctl_ub   = "la_check_ctl_ub"
    CHARACTER(LEN=ip_cnl) :: la_check_ctl_ubu  = "N/A",&
    CHARACTER(LEN=ip_lnl) :: la_check_ctl_ubln = "Check upper bound of acceptable values for the control variable"
    
		!cordinate variables
		CHARACTER(LEN=ip_cnl) :: x	= "x"
		CHARACTER(LEN=ip_cnl) :: xu	= "m"
		CHARACTER(LEN=ip_lnl) :: xln= "Longitude"
		CHARACTER(LEN=ip_cnl) :: y	= "y"
		CHARACTER(LEN=ip_cnl) :: yu	= "m"
		CHARACTER(LEN=ip_lnl) :: yln= "Latitude" 
		CHARACTER(LEN=ip_cnl) :: z	= "y"
		CHARACTER(LEN=ip_cnl) :: zu	= "m"
		CHARACTER(LEN=ip_lnl) :: zln= "Altitude"   
  END TYPE exchange_param_attributes

  TYPE obsOutput
    CHARACTER(LEN=ip_snl) :: filename
    INTEGER :: i_max_nobs = -1
    INTEGER :: ncid = -1
    INTEGER :: obs_id	= -1
    INTEGER :: sigma_id	= -1
    INTEGER :: Rm1_id	= -1
    INTEGER :: obsgap_id  = -1! concentration id
    INTEGER :: M_vactor_id = -1
    INTEGER :: icoord_id = -1
    INTEGER :: rcoord_id = -1
    INTEGER	:: timeStep_id = -1
    INTEGER	:: date_id = -1
    !champs pour la sauvegarde des attributs essentiels de la trajectoire
    INTEGER :: i_sub = 1 !fr�quence de sauvegarde
    INTEGER :: nb_record = -1 !number of record in the file
    
    REAL(KIND=dp)  :: delta_t = -1.0_dp! pas de temp
    REAL(KIND=dp)  :: delta_x = -1.0_dp! pas en x
    REAL(KIND=dp)  :: delta_y = -1.0_dp! pas en y
    LOGICAL	:: isOpened = .FALSE.
  END TYPE obsOutput
  
  !> \brief User defined type for attributes of observations
  !<
  TYPE obs_structure_attributes
    CHARACTER(LEN=ip_cnl) ::&
      obs_fName  = "obs_fName",&
      ogap_fName = "ogap_fName",&
      i_max_nobs = "i_max_nobs",&
      i_nb_dim = "i_nb_dim",&
      i_nobs   = "i_nobs",&
      l_rcoord = "l_rcoord",&
      l_icoord = "l_icoord",&
      l_icoord_allocated = "l_icoord_allocated",&
      l_rcoord_allocated = "l_rcoord_allocated"
      
    CHARACTER(LEN=ip_cnl) ::&
      ra_obs    = "ra_obs",&
      ra_obsu   = "N/A",&
      ra_obsln  = "observation data"
      
    CHARACTER(LEN=ip_cnl) ::&
      ra_sigma   = "ra_sigma",&
      ra_sigmau  = "N/A",&
      ra_sigmaln = "standard deviation of observation error"
      
    CHARACTER(LEN=ip_cnl) ::&
      ra_Rm1    = "ra_Rm1",&
      ra_Rm1u   = "N/A",&
      ra_Rm1ln  = "inverse covariance (diag matrix)"
      
    CHARACTER(LEN=ip_cnl) ::&
      ra_obsgap   = "ra_obsgap",&
      ra_obsgapu  = "N/A",&
      ra_obsgapln = "Diference between obs and model output"

    !wavelet related variables
    CHARACTER(LEN=ip_cnl) :: i_obs_level = "i_obs_level"
    CHARACTER(LEN=ip_cnl) ::&
      ia_M_vector = "ia_M_vector",&
      ia_M_vectoru = "N/A",&
      ia_M_vectorln = "Size of the wavelet grid at the coarsest scale"
      
    CHARACTER(LEN=ip_cnl) ::&
      ia_icoord = "ia_icoord",&
      ia_icoordu = "N/A",&
      ia_icoordln = "Integer coordinates, indices in the discretization grid"
      
    CHARACTER(LEN=ip_cnl) ::&
      ra_rcoord   = "ra_rcoord",&
      ra_rcoordu  = "N/A",&
      ra_rcoordln = "real coordinates in the computation domain"
  END TYPE obs_structure_attributes
  
CONTAINS

  
END MODULE nccom_tools