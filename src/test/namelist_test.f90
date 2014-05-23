PROGRAM namelist_test
  USE general_constant
  USE com_tools
  USE debug_tools
IMPLICIT NONE
  INTEGER, PARAMETER :: IP_MAX = 5,&
    ip_numnam = 63
  INTEGER :: il_scalar, il_size, ibi
  INTEGER, DIMENSION(IP_MAX) :: ila_array
  REAL(KIND=cp) :: rl_scalar
  REAL(KIND=cp), DIMENSION(IP_MAX) :: rla_array
  CHARACTER(LEN=ip_fnl) :: ala_history, ala_fName, ala_test


  NAMELIST/NAM_history/ala_history, ala_test
  NAMELIST/NAM_size/il_size
  NAMELIST/NAM_data/il_scalar, rl_scalar, ila_array, rla_array

  il_size = IP_MAX
  DO ibi = 1, IP_MAX
    ila_array(ibi) = ibi
    rla_array(ibi) = REAL(ibi, cp)
  END DO
  ala_fName = 'runtime_namelist.dat'
  ala_test = 'test numero 2'
  !CALL save_namelist(ala_fName)
  !CALL debug('', 'Saved values===============================================')
  !CALL print_namelist()
  CALL load_namelist(ala_fName)
  CALL debug('', 'Loaded values==============================================')
  CALL print_namelist()
CONTAINS

  SUBROUTINE save_namelist(ada_fName)
    CHARACTER(LEN=ip_fnl), INTENT(IN) :: ada_fName

    !initializing history
    ala_history = make_history('namelist_test')
		CALL debug(ala_history, 'In save_namelist, ala_history = ')

    OPEN(ip_numnam, FILE=ada_fName, FORM=FILE_FORMATTED, STATUS=FILE_REPLACE)
    WRITE(ip_numnam, NAM_history)! writing the history block
    WRITE(ip_numnam, NAM_size)!writing the size block
    WRITE(ip_numnam, NAM_data)!writing the size block
    CLOSE(ip_numnam)
  END SUBROUTINE save_namelist

  SUBROUTINE load_namelist(ada_fName)
    CHARACTER(LEN=ip_fnl), INTENT(IN) :: ada_fName

    OPEN(ip_numnam, FILE=ada_fName, FORM=FILE_FORMATTED, STATUS=FILE_OLD)
    READ(ip_numnam, NAM_history)! writing the history block
    READ(ip_numnam, NAM_size)!writing the size block
    READ(ip_numnam, NAM_data)!writing the size block
    CLOSE(ip_numnam)
  END SUBROUTINE load_namelist

  SUBROUTINE print_namelist()
    CALL debug('', 'In print_namelist, namelist values ++++++++++++++++++++++++++++')
    CALL debug(il_size, 'il_size = ')
    CALL debug(il_scalar, 'il_scalar = ')
    CALL debug(rl_scalar, 'rl_scalar = ')
    CALL debug(ila_array, 'ila_array = ')
    CALL debug(rla_array, 'rla_array = ')
    CALL debug(ala_test, 'ala_test = ')
    CALL debug(ala_history, 'ala_history = ')
    CALL debug('', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  END SUBROUTINE print_namelist

END PROGRAM namelist_test