!>file ctl_tools.f90
!!\brief CTL tools tools
!<
MODULE ctl_tools
  USE general_constant
  USE general_tools
  USE debug_tools
IMPLICIT NONE


  !> \brief User defined type for control vector parameters
  !! Defines parameters that are used to define or manage the control vector
  !<
  TYPE bctl_param
    CHARACTER(LEN=ip_snl) ::&
      aa_title,& !> title of the output file
      aa_hist    !> history string for the output file
    REAL(cp) :: r_sigmaB
  END TYPE bctl_param

CONTAINS

  !> \brief load ctl parameters from a namelist file
  !! \param[in,out] td_bgp data structure to be loaded
  !! \param[in] ada_namelist name of the namelist file
  !! \param[in] ada_varName name of the variable associated with td_bgp
  !!
  !<
  SUBROUTINE load_bctlp(td_bgp, ada_namelist, ada_varName)
    INTEGER, PARAMETER:: ip_numnam = 68
    TYPE(bctl_param), INTENT(IN OUT) :: td_bgp
    CHARACTER(len=*), INTENT(IN) :: ada_namelist, ada_varName
    !local variables
    REAL(cp) :: rl_sigmaB
    CHARACTER(len=ip_snl) :: varName, title, history, ala_varName
    NAMELIST/NAM_BCTL_PARAM/&
      varName  ,&
      title    ,&
      history  ,&
      rl_sigmaB

    OPEN(ip_numnam,FILE=ada_namelist,FORM='FORMATTED',STATUS='OLD')
    varName = "{#}@"
    CALL debug(100, "In load_bctlp --- ", tag=dALLWAYS)
    ala_varName = TRIM(ada_varName)
    CALL uppercase(ala_varName)
    CALL debug(200, "In load_bctlp --- ", tag=dALLWAYS)
    DO WHILE(varName.NE.ala_varName)
      READ(ip_numnam, NAM_BCTL_PARAM)!reading the block
      CALL uppercase(varName)
      CALL debug(varName, "In load_bctlp, found varName = ", tag=dALLWAYS)
    END DO
    CLOSE(ip_numnam)
    td_bgp%aa_title = title
    td_bgp%aa_hist  = history
    td_bgp%r_sigmaB = rl_sigmaB

  END SUBROUTINE load_bctlp

  !> \brief print ctl parameters
  !! \param[in] td_bgp data structure for ctl parameters
  !!
  !<
  SUBROUTINE print_bctlp(td_bgp)
    TYPE(bctl_param), INTENT(IN) :: td_bgp

    CALL debug('', 'printing convObs_param---------------------------------')
    CALL debug(td_bgp%aa_title   , '  aa_title = ')
    CALL debug(td_bgp%aa_hist    , '  aa_hist  = ')
    CALL debug(td_bgp%r_sigmaB   , '  r_sigmaB = ')
    CALL debug('', '.......................................................')
  END SUBROUTINE print_bctlp

END MODULE ctl_tools