!>
!!Routines to write and read matrix to/from file. These routines work in column order
!<
MODULE io_tools
	USE general_constant

IMPLICIT NONE

CONTAINS

  SUBROUTINE write_matrix(fileName, rda_A)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN) :: rda_A
    CHARACTER*(ip_snl)         , INTENT(IN) :: fileName
    INTEGER ib_i, ib_j, il_ios
    INTEGER, PARAMETER :: ip_fid = 41

    ip_fid = 41
    OPEN( UNIT = ip_fid, FILE = fileName, STATUS = 'REPLACE', FORM = 'FORMATTED', IOSTAT = il_ios)
    IF (il_ios /= 0) then
      WRITE(* , *) 'Error creating file', fileName
      STOP;
    END IF
    WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 1)
    WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 2)
    DO ib_j = LBOUND(rda_A,2), UBOUND(rda_A,2)
    	DO ib_i = LBOUND(rda_A,1), UBOUND(rda_A,1)
        WRITE(unit=ip_fid, fmt='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid end of line
      END DO
      WRITE(unit=ip_fid, fmt=*)''!
    END DO
    CLOSE(ip_fid )
  END SUBROUTINE write_matrix

  SUBROUTINE read_matrix(fileName, rda_A)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT) :: rda_A
    CHARACTER*(80)         , INTENT(IN)  :: fileName
    INTEGER ib_i, ib_j, il_nbRow, il_nbCol, il_ios
    INTEGER , PARAMETER :: ip_fid =615

    OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
    IF (il_ios /= 0) then
      WRITE(* , *) 'Error opening file', fileName
      STOP;
    END IF
    !stop
    READ(UNIT=ip_fid, FMT=IFMT) il_nbRow
    READ(UNIT=ip_fid, FMT=IFMT) il_nbCol
    PRINT*, "In read_matrix, il_nbRow = ", il_nbRow, "; il_nbCol = ", il_nbCol
    DO ib_j = 1, il_nbCol
    	DO ib_i = 1, il_nbRow
        READ(UNIT=ip_fid, FMT='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid going to the next line
      END DO
      READ(unit=ip_fid, fmt=*) !just to skip the end of line
    END DO
    CLOSE(ip_fid )
  END SUBROUTINE read_matrix

	!!read the head of the file: nbRow and nbCol
  SUBROUTINE readInfo(fileName, id_nbRow, id_nbCol)
    INTEGER , PARAMETER :: ip_fid =616
    CHARACTER*(80)         , INTENT(IN)  :: fileName
    INTEGER, INTENT(OUT) :: id_nbRow, id_nbCol
    INTEGER :: il_ios

    OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
    IF (il_ios /= 0) then
      WRITE(* , *) 'In readInfo : Error opening file ', fileName
      STOP;
    END IF
    READ(UNIT=ip_fid, FMT=IFMT) id_nbRow
    READ(UNIT=ip_fid, FMT=IFMT) id_nbCol
    PRINT*, "In readInfo, id_nbRow = ", id_nbRow, "; id_nbCol = ", id_nbCol
    CLOSE(ip_fid )
  END SUBROUTINE readInfo

END MODULE io_tools
