!> \file checkpoint.f90
!! @brief Checpointing module routines for BALAISE, used mainly for generalised diffusion projection
!! @author Innocent Souopgui
!! @version 1.0
MODULE checkpoint
  USE debug_tools
IMPLICIT NONE

  !>
  LOGICAL, PRIVATE :: lm_checkpoint_on_disk = .TRUE.
  !>@brief File Id used in this module for checkpointing
  INTEGER, PRIVATE, PARAMETER :: im_nunit = 10

  !> @brief User define type for checkpoint
  TYPE checkpoint_type
    INTEGER :: i_unit = -1 !>
    INTEGER :: i_nRecord = -1!> maximum number of record to be saved
    INTEGER :: i_varSize = -1!> number of element in a single variable
    INTEGER :: i_recSize = -1!> Size of each record in bytes
    INTEGER :: i_nVar = -1!> number of variables in each record
    CHARACTER(len=ip_snl):: aa_fName
    !
  END TYPE checkpoint_type

  TYPE(checkpoint_type), PRIVATE, SAVE, DIMENSION(im_nunit):: tm_chkp

  !> @brief interface for opening checkpoint unit
  INTERFACE chkp_open
    MODULE PROCEDURE chkp_open_for_write
    MODULE PROCEDURE chkp_open_for_read
  END INTERFACE chkp_open

  !> @brief interface for writing/saving checkpoint data
  INTERFACE chkp_save
    MODULE PROCEDURE chkp_save_1v_1d!one variable on 1d
    MODULE PROCEDURE chkp_save_1v_2d!one variable on 2d
!     MODULE PROCEDURE chkp_save_2v_2d!two variables on 2d
  END INTERFACE chkp_save

  !> @brief interface for reading/restauring checkpoint data
  INTERFACE chkp_restaure
    MODULE PROCEDURE chkp_restaure_1v_1d!one variable on 1d
    MODULE PROCEDURE chkp_restaure_1v_2d !one variable on 2d
!     MODULE PROCEDURE chkp_restaure_2v_2d !two variables on 2d
  END INTERFACE chkp_restaure

CONTAINS

  !> @brief Check the id of the checkpointing unit
  !! @param [in] id_chkpId checkpoint identifier, this must be an integer between 1 and im_nunit
  !! @details  Causes the program to stop if the id_chkpId is out of bounds
  !<
  SUBROUTINE chkp_check_id(id_chkpId)
    INTEGER, INTENT(IN) :: id_chkpId

    IF (id_chkpId > im_nunit) THEN
      CALL debug('', 'In chkp_prepare_param, the id of the checkpointing unit is greater than the maximum allowed')
      CALL debug(id_chkpId, '  The id of the checkpointing unit is:')
      CALL debug(im_nunit, '  The maximum allowed for the id of a checkpointing unit is:')
      CALL stop_program('', 'The maximum is set by the parameter variable <im_nunit> in <gd_regul.f90>')
    ELSE IF(id_chkpId <= 0) THEN
      CALL debug('', 'In chkp_prepare_param, the id of the checkpointing unit must be a positive integer')
      CALL stop_program(id_chkpId, '  The id of the checkpointing unit is:')
    END IF
  END SUBROUTINE chkp_check_id

  !> @brief Prepares and checks parameters to open a checkpointing unit
  !! @param [in] id_chkpId checkpoint identifier, this must be an integer between 1 and im_nunit
  !! @param [in] id_nVar Number of variables to be saved in the checkpointing unit
  !! @param [in] id_varSize Size of each variable
  !! @param [in] nrec (optional) maximum number of record in the checkpointing unit, this is mandatory only if checkpoint is saved in memory
  !! @details  If one needs more checkpoint units, increase the value of the variable im_nunit
  !! Causes the program to stop if the id_chkpId is out of bounds
  !<
  SUBROUTINE chkp_prepare_param(id_chkpId, id_nVar, id_varSize, nrec)
    INTEGER, INTENT(IN) :: id_chkpId, id_nVar, id_varSize
    INTEGER, INTENT(IN), OPTIONAL :: nrec

    CALL chkp_check_id(id_chkpId)
    tm_chkp(id_chkpId)%i_nVar    = id_nVar
    tm_chkp(id_chkpId)%i_varSize = id_varSize
    tm_chkp(id_chkpId)%i_recSize = cp*id_nVar*id_varSize
    IF (lm_checkpoint_on_disk) THEN
      tm_chkp(id_chkpId)%i_unit  = 1300 + id_chkpId !1300 is the base for checkpointing unit, can be changed to any other positive integer. Just make sure that ther will be no conflict of file unit
      WRITE(tm_chkp(id_chkpId)%aa_fName, FMT='(A,I5.5,A)') 'runtime_gd_chkp_', id_chkpId, '.tmp'
    ELSE
      IF( PRESENT(nrec) ) THEN
        tm_chkp(id_chkpId)%i_nRecord = nrec
      ELSE
        tm_chkp(id_chkpId)%i_nRecord = -1
      END IF
    END IF
  END SUBROUTINE chkp_prepare_param

  !> @brief Open a checkpoint unit for write
  !! @param [in] id_chkpId checkpoint identifier, this must be an integer between 1 and im_nunit
  !! @param [in] id_nVar Number of variables to be saved in the checkpointing unit
  !! @param [in] id_varSize Size of each variable
  !! @param [in] nrec (optional) maximum number of record in the checkpointing unit, this is mandatory only if checkpoint is saved in memory
  !! @details  If one needs more checkpoint units, increase the value of the variable im_nunit
  !! This routine allocates space for memory checkpointing if necessary
  !! See checkpoint module for details on im_nunit
  !<
  SUBROUTINE chkp_open_for_write(id_chkpId, id_nVar, id_varSize, nrec)
    INTEGER, INTENT(IN) :: id_chkpId, id_nVar, id_varSize
    INTEGER, INTENT(IN), OPTIONAL :: nrec
    INTEGER :: il_iostat

    CALL chkp_prepare_param(id_chkpId, id_nVar, id_varSize, nrec)
    IF (lm_checkpoint_on_disk) THEN
      OPEN(UNIT=tm_chkp(id_chkpId)%i_unit, FILE=tm_chkp(id_chkpId)%aa_fName, ACTION="write", ACCESS='direct',&
           RECL=tm_chkp(id_chkpId)%i_recSize, STATUS='replace', FORM='unformatted', IOSTAT=il_iostat&
          )
      IF(il_iostat>0)THEN
        CALL stop_program(il_iostat, 'In chkp_open_for_write, error creating checkpointing file; IOSTAT = ')
      END IF
    ELSE
      CALL stop_program('', 'In chkp_open_for_write: checkpointing in memory is not yet available.')
    END IF
  END SUBROUTINE chkp_open_for_write

  !> @brief Open a checkpoint unit for read
  !! @param [in] id_chkpId checkpoint identifier, this must be an integer between 1 and im_nunit
  !! @details  This routine assumes the file or memory space already exists
  !<
  SUBROUTINE chkp_open_for_read(id_chkpId)
    INTEGER, INTENT(IN) :: id_chkpId
    INTEGER :: il_iostat

    CALL chkp_check_id(id_chkpId)
    IF (lm_checkpoint_on_disk) THEN
      OPEN( UNIT=tm_chkp(id_chkpId)%i_unit, FILE=tm_chkp(id_chkpId)%aa_fName, ACTION="read", ACCESS='direct',&
            RECL=tm_chkp(id_chkpId)%i_recSize, STATUS='old', FORM='unformatted', IOSTAT=il_iostat&
          )
      IF(il_iostat>0)THEN
        CALL debug(il_iostat, 'In chkp_open_for_read, error opening checkpointing file; IOSTAT = ')
        CALL debug(id_chkpId, 'Make sure that the file has been created, Checkpointing id = ')
        CALL stop_program(tm_chkp(id_chkpId)%aa_fName, '  File Name = ')
      END IF
    ELSE
      CALL stop_program('', 'In chkp_open_for_read: checkpointing in memory is not yet available.')
    END IF
  END SUBROUTINE chkp_open_for_read

  !> @brief Close a checkpoint unit
  !! @param [in] id_chkpId checkpoint identifier, this must be an integer between 1 and im_nunit
  !! @details  see checkpoint module for details
  !<
  SUBROUTINE chkp_close(id_chkpId)
    INTEGER, INTENT(IN) :: id_chkpId
    INTEGER :: il_iostat

    IF (lm_checkpoint_on_disk) THEN
      CLOSE( UNIT=tm_chkp(id_chkpId)%i_unit, IOSTAT=il_iostat )
      IF(il_iostat>0)THEN
        CALL debug(il_iostat, 'In chkp_close, error closing checkpointing file; IOSTAT = ')
        CALL stop_program(id_chkpId, 'Make sure that the file has been created, Checkpointing id = ')
      END IF
    ELSE
      CALL stop_program('', 'In chkp_close: checkpointing in memory is not yet available.')
    END IF
  END SUBROUTINE chkp_close

  !> @brief save vector (for adjoint calculations)
  !! @param [in] id_chkpId identifier of the unit where to save
  !! @param [in] rda_u intermediate result to be saved
  !! @param [in] id_rec record number to save the intermediate result
  !! \todo The project is to have two options for saving intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_save_vec(id_chkpId, rda_u, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: rda_u

    IF (lm_checkpoint_on_disk) THEN
      WRITE(UNIT=tm_chkp(id_chkpId)%i_unit, REC=id_rec) rda_u
    ELSE
      CALL stop_program('', 'In chkp_save_vec:  the checkpoint in memory is not yet operational')
    END IF
  END SUBROUTINE chkp_save_vec

  !> @brief restaures intermediate results (one variable, one dimension) for adjoint calculations
  !! @param [in] id_chkpId identifier of the unit where the result had been saved
  !! @param [in, out] rda_vec contained the restaure intermediate results
  !! @param [in] id_rec record number to restaure
  !! \todo The project is to have two options for restauring intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_restaure_vec(id_chkpId, rda_vec, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)   :: rda_vec

    IF (lm_checkpoint_on_disk) THEN
      READ(UNIT=tm_chkp(id_chkpId)%i_unit, REC=id_rec) rda_vec
    ELSE
      CALL stop_program('', 'In chkp_restaure_vec:  the checkpoint in memory is not yet operational')
    END IF
  END SUBROUTINE chkp_restaure_vec

  !> @brief save intermediate results of GD projection (one variable, one dimension) for adjoint calculations
  !! @param [in] id_chkpId identifier of the unit where to save
  !! @param [in] rda_u intermediate result to be saved
  !! @param [in] id_rec record number to save the intermediate result
  !! \todo The project is to have two options for saving intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_save_1v_1d(id_chkpId, rda_u, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: rda_u

    CALL chkp_save_vec(id_chkpId, rda_u, id_rec)
  END SUBROUTINE chkp_save_1v_1d

  !> @brief save intermediate results of GD projection (one variable, two dimensions) for adjoint calculations
  !! @param [in] id_chkpId identifier of the unit where to save
  !! @param [in] rda_u intermediate result to be saved
  !! @param [in] id_rec record number to save the intermediate result
  !! \todo The project is to have two options for saving intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_save_1v_2d(id_chkpId, rda_u, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN)   :: rda_u
    REAL(KIND=dp), DIMENSION( PRODUCT(SHAPE(rda_u)) )   :: rla_vec

    rla_vec = RESHAPE(rda_u, SHAPE(rla_vec) )
    CALL chkp_save_vec(id_chkpId, rla_vec, id_rec)
  END SUBROUTINE chkp_save_1v_2d

  !> @brief restaures intermediate results (one variable, one dimension) for adjoint calculations
  !! @param [in] id_chkpId identifier of the unit where the result had been saved
  !! @param [in, out] rda_u contained the restaure intermediate results
  !! @param [in] id_rec record number to restaure
  !! \todo The project is to have two options for restauring intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_restaure_1v_1d(id_chkpId, rda_u, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)   :: rda_u

    CALL chkp_restaure_vec(id_chkpId, rda_u, id_rec)
  END SUBROUTINE chkp_restaure_1v_1d

  !> @brief restaures intermediate results (one variable, two dimensions) for adjoint calculations
  !! @param [in] id_chkpId identifier of the unit where the result had been saved
  !! @param [in, out] rda_u contained the restaure intermediate results
  !! @param [in] id_rec record number to restaure
  !! \todo The project is to have two options for restauring intermediate results : 1st option to file and second option in memory whith total transparency for the user
  !! @details  The value of id_chkpId must be between 1 and im_nunit, see checkpoint module for details
  !<
  SUBROUTINE chkp_restaure_1v_2d(id_chkpId, rda_u, id_rec)
    INTEGER, INTENT(IN) :: id_chkpId, id_rec
    REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT)   :: rda_u
    REAL(KIND=dp), DIMENSION( PRODUCT(SHAPE(rda_u)) )   :: rla_vec

    CALL chkp_restaure_vec(id_chkpId, rla_vec, id_rec)
    rda_u = RESHAPE(rla_vec, SHAPE(rda_u) )
  END SUBROUTINE chkp_restaure_1v_2d

END MODULE checkpoint