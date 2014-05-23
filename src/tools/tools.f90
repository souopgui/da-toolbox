!> \file tools.f90
!! \brief general tools definition
!! @author Innocent Souopgui
!! \details
!! 1D dynamic array, 2D dynamic array. The 2D Dynamic arrays look like java 2D, it is a 1D array of 1D array.
!<
MODULE tools
    USE general_constant
IMPLICIT NONE

     INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE assign_drv_vector!> assign 1D array to drv
        MODULE PROCEDURE assign_div_vector!> assign 1D array to div
        MODULE PROCEDURE assign_dsv_vector!> assign 1D array to dcv
        MODULE PROCEDURE assign_2DdrA_vector!> assign 1D array to d2Da
        MODULE PROCEDURE assign_2DdrA_2DA!> assign 2D array to d2Da
     END INTERFACE

    !> \brief sort routines, selection
    !<
    INTERFACE select_sort
      MODULE PROCEDURE select_sort_int
      !MODULE PROCEDURE select_sort_sp_real
      !MODULE PROCEDURE select_sort_dp_real
    END INTERFACE
    !> \brief switch routines
    !<
    INTERFACE switch
      MODULE PROCEDURE switch_int
      !MODULE PROCEDURE switch_sp_real
      !MODULE PROCEDURE switch_dp_real
    END INTERFACE

!      INTERFACE MINVAL
!         MODULE PROCEDURE d2DrA_minval!> minval in dynamic 2D real vector
!         MODULE PROCEDURE drv_minval!> minval in dynamic 1D real vector
!         MODULE PROCEDURE drv_array_minval!> minval of an array of dynamic 1D real vector
!      END INTERFACE MINVAL
!
!      INTERFACE MAXVAL
!         MODULE PROCEDURE d2DrA_maxval!> maxval in dynamic 2D real vector
!         MODULE PROCEDURE drv_maxval!> maxval in dynamic 1D real vector
!         MODULE PROCEDURE drv_array_maxval!> minval of an array of dynamic 1D real vector
!      END INTERFACE MAXVAL

     !INTERFACE NULLIFY
     !   MODULE PROCEDURE drv_nullify!> nullify dynamic real vector
     !   MODULE PROCEDURE dov_nullify!> nullify dynamic observation vector
     !END INTERFACE NULLIFY

    INTEGER, PARAMETER :: OBS_FIELD = 0!> constant for data field of dynamic observation vector
    INTEGER, PARAMETER :: OBS_CMAT  = 1!> constant for covariance matrix field of dynamic observation vector
    INTEGER, PARAMETER :: OBS_XIDX  = 2!> constant for H field of dynamic observation vector

    !> \brief dynamic real vector
    TYPE dyn_rVector
        REAL(KIND=cp), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_rVector

    !> \brief dynamic integer vector
    TYPE dyn_iVector
        INTEGER, DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_iVector

    !> \brief character vector
    TYPE dyn_sVector
        CHARACTER(LEN=ip_snl), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_sVector

    !> \brief 2D dynamic real array (array of arrays), java style 2D array
    TYPE dyn_2DrArray
        TYPE(dyn_rVector), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_2DrArray

    !> \brief 2D dynamic integer array (array of arrays), java style 2D array
    TYPE dyn_2DiArray
        TYPE(dyn_iVector), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_2DiArray

    !> \brief observation vector at a given date
    !!
    !<
    TYPE dyn_oVector
        INTEGER :: ts !> time step
        LOGICAL, PRIVATE :: H_allocated = .FALSE. !> is memory space allocated for idx
        LOGICAL, PRIVATE :: Rm1_allocated = .FALSE. !> is memory space allocated for Rm1
        REAL(KIND=cp), DIMENSION(:), POINTER :: dat  => NULL() !> observation data
        INTEGER      , DIMENSION(:), POINTER :: H => NULL()!> observation operator, indexes of obs data in the spatial grid
        !> inverse of the observation diagonal covariance matrix, can be extended to full matrix with covariates components
        !<
        REAL(KIND=cp), DIMENSION(:), POINTER :: Rm1    => NULL()
    END TYPE dyn_oVector
CONTAINS
    SUBROUTINE select_sort_int(ida_x)
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ida_x
      INTEGER :: ibi, ibj, il_min_idx, il_size

      il_size = SIZE(ida_x)
      DO ibi = 1, il_size-1
        il_min_idx = ibi
        DO ibj = ibi+1, il_size
          IF( ida_x(ibj)<ida_x(il_min_idx) ) il_min_idx = ibj
        END DO
        IF(ibi /= il_min_idx) THEN
          CALL switch( ida_x(ibi), ida_x(il_min_idx) )
        END IF
      END DO
    END SUBROUTINE select_sort_int

    SUBROUTINE switch_int(id_a, id_b)
      INTEGER, INTENT(INOUT) :: id_a, id_b
      !local variables
      INTEGER :: il_tmp

      il_tmp = id_a
      id_a = id_b
      id_b = il_tmp
    END SUBROUTINE switch_int

    !> \brief assign observation data to a dynamic observation vector
    !! \param[in, out] td_dra dynamic observation vector (structure)
    !! \param[in] rda_data observation data
    !! \param[in] ida_H indexes of observation data in spatial grid
    !! \param[in] rda_Rm1 diagonal matrix, inverse of the observation's covariance error (can be extended to 2D)
    !! \details ida_idx and rda_Rm1 are optional since they can be identic for all observation date in time
    !<
    SUBROUTINE assign_dov_vector(td_dov, rda_data, rda_Rm1, ida_H)
        TYPE(dyn_oVector), INTENT(INOUT)        :: td_dov
        REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: rda_Rm1
        INTEGER,       DIMENSION(:), OPTIONAL, INTENT(IN) :: ida_H
        LOGICAL :: ll_idx, ll_Rm1

        IF( PRESENT(ida_H) )THEN
            ll_idx = .TRUE.
        ELSE
            ll_idx = .FALSE.
        END IF
        IF( PRESENT(rda_Rm1) )THEN
            ll_Rm1 = .TRUE.
        ELSE
            ll_Rm1 = .FALSE.
        END IF
        CALL reset_dov(td_dov, SIZE(rda_data), ll_idx, ll_Rm1)
        td_dov%dat = rda_data
        IF(ll_idx) td_dov%H = ida_H
        IF(ll_Rm1) td_dov%Rm1 = rda_Rm1
    END SUBROUTINE assign_dov_vector

    !> \brief nullify pointer component of dynamic observation vector
    !! \param[in, out] td_dov observation vector to be nullify
    !<
    SUBROUTINE dov_nullify(td_dov)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov

        td_dov%dat => NULL()
        td_dov%Rm1 => NULL()
        td_dov%Rm1_allocated = .FALSE.
        td_dov%H => NULL()
        td_dov%H_allocated = .FALSE.
    END SUBROUTINE dov_nullify

    !< \brief Sample observation from state variable
    !! \param[in] rda_x state variable
    !! \param[in, out] rda_yo observation variable
    !! \param[in] rda_e additive error (same shape as rda_yo) to corrupt obs
    !>
    SUBROUTINE sample_1D_direct_obs(rda_x, rda_idx, rda_yo, rda_e)
        REAL(KIND=cp), DIMENSION(:), INTENT(IN)           :: rda_x
        INTEGER, DIMENSION(:), INTENT(IN)                 :: rda_idx
        REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)        :: rda_yo
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: rda_e
        INTEGER :: ibi

        DO ibi = 1, SIZE(rda_yo)
            rda_yo(ibi) = rda_x( rda_idx(ibi) )
        END DO
        IF (PRESENT(rda_e)) rda_yo = rda_yo + rda_e
    END SUBROUTINE sample_1D_direct_obs

    !>resize the observation vector structure
    !! \param[in, out] td_dov observation vector to be resized
    !! \param[in] id_size new size of the observation vector
    !! \param[in] ld_Rm1 says if the field R (covariance) need to be allocated, default is .FALSE.
    !! \param[in] ld_H says if the field H (projection obs operator) need to be allocated, default is .FALSE.
    !<
    SUBROUTINE reset_dov(td_dov, id_size, ld_Rm1, ld_H)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov
        INTEGER, INTENT(IN) :: id_size
        LOGICAL, INTENT(IN), OPTIONAL :: ld_Rm1, ld_H
        LOGICAL :: ll_Rm1, ll_H

        IF ( PRESENT(ld_Rm1) ) THEN
            ll_Rm1 = ld_Rm1
        ELSE
            ll_Rm1 = .FALSE.
        END IF
        IF ( PRESENT(ld_H) ) THEN
            ll_H = ld_H
        ELSE
            ll_H = .FALSE.
        END IF

        !dat field
        CALL reset_dov_field(td_dov, id_size, OBS_FIELD)
        !R field
        IF(ll_Rm1)THEN
            CALL reset_dov_field(td_dov, id_size, OBS_CMAT)
        ELSE
            CALL reset_dov_field(td_dov, 0, OBS_CMAT)!resize to zero
        END IF
        !H field
        IF(ll_H)THEN
            CALL reset_dov_field(td_dov, id_size, OBS_XIDX)
        ELSE
            CALL reset_dov_field(td_dov, 0, OBS_XIDX)!resize to zero
        END IF
    END SUBROUTINE reset_dov

    !>resize a field of the observation vector structure
    !! \param[in, out] td_dov observation vector to be resized
    !! \param[in] id_size new size of the observation vector
    !! \param[in] id_field identificator for the field to be resized. accepted values : OBS_FIELD (for data field), OBS_CMAT (for covariance matrix), OBS_XIDX (for x indices)
    !!
    !<
    SUBROUTINE reset_dov_field(td_dov, id_size, id_field)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov
        INTEGER, INTENT(IN) :: id_size, id_field
        INTEGER :: il_cSize

        !dat field
        IF( ASSOCIATED(td_dov%dat) )THEN
            il_cSize = SIZE(td_dov%dat)
        ELSE
            il_cSize = 0
        END IF

        SELECT CASE (id_field)
            CASE (OBS_FIELD)
                IF (il_cSize /= id_size) THEN !current size different from required size
                    IF(il_cSize > 0) THEN !current size greater than zero
                        DEALLOCATE(td_dov%dat)
                        NULLIFY(td_dov%dat)
                    END IF
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%dat(id_size) )
                        td_dov%dat = 0.0_cp
                    END IF
                END IF
            CASE (OBS_CMAT)
                IF(td_dov%Rm1_allocated)THEN!R has been allocated
                    IF (il_cSize /= id_size) THEN!current size different from required size
                        IF(il_cSize > 0) THEN!current size greater than zero
                            DEALLOCATE(td_dov%Rm1)
                            NULLIFY(td_dov%Rm1)
                            td_dov%Rm1_allocated = .FALSE.
                        END IF
                        IF(id_size>0)THEN!required size greater than zero
                            ALLOCATE( td_dov%Rm1(id_size) )
                            td_dov%Rm1 = 0.0_cp
                            td_dov%Rm1_allocated = .TRUE.
                        END IF
                    END IF
                ELSE!R has not been allocated
                    NULLIFY(td_dov%Rm1)
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%Rm1(id_size) )
                        td_dov%Rm1 = 0.0_cp
                        td_dov%Rm1_allocated = .TRUE.
                    END IF
                END IF
            CASE (OBS_XIDX)
                IF(td_dov%H_allocated)THEN!H has been allocated
                    IF (il_cSize /= id_size) THEN!current size different from required size
                        IF(il_cSize > 0) THEN!current size greater than zero
                            DEALLOCATE(td_dov%H)
                            NULLIFY(td_dov%H)
                            td_dov%H_allocated = .FALSE.
                        END IF
                        IF(id_size>0)THEN!required size greater than zero
                            ALLOCATE( td_dov%H(id_size) )
                            td_dov%H = 0
                            td_dov%H_allocated = .TRUE.
                        END IF
                    END IF
                ELSE!H has not been allocated
                    NULLIFY(td_dov%H)
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%H(id_size) )
                        td_dov%H = 0
                        td_dov%H_allocated = .TRUE.
                    END IF
                END IF
        END SELECT
    END SUBROUTINE reset_dov_field

    SUBROUTINE assign_drv_vector(td_dra, rd_data)
        TYPE(dyn_rVector), INTENT(INOUT) :: td_dra
        REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rd_data

        CALL reset_drv(td_dra, SIZE(rd_data))
        td_dra%dat = rd_data
    END SUBROUTINE assign_drv_vector

    !> \brief nullify pointer component of dynamic real vector
    !! \param[in, out] td_drv dynamic real vector to be nullify
    !<
    SUBROUTINE drv_nullify(td_drv)
        TYPE(dyn_rVector), INTENT(INOUT) :: td_drv

        NULLIFY(td_drv%dat)
    END SUBROUTINE drv_nullify

    SUBROUTINE reset_drv_array(tda_drv)
      TYPE(dyn_rVector), DIMENSION(:), INTENT(INOUT) :: tda_drv
      !local variables
      INTEGER :: ibi

      DO ibi=1,SIZE(tda_drv)
        CALL reset_drv( tda_drv(ibi), 0 )
      END DO
    END SUBROUTINE reset_drv_array

    SUBROUTINE reset_drv(td_dra, id_size, id_lb)
        TYPE(dyn_rVector), INTENT(INOUT) :: td_dra
        INTEGER, INTENT(IN) :: id_size
        INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound in the array if needed
        INTEGER :: il_lb, il_stat
        LOGICAL :: ll_allocate

        IF ( PRESENT(id_lb) ) THEN
            il_lb = id_lb
        ELSE
            il_lb = 1
        END IF

        IF( ASSOCIATED(td_dra%dat) )THEN
            IF( (id_size==SIZE(td_dra%dat)).AND.(il_lb == LBOUND(td_dra%dat, 1) ) )THEN
                ll_allocate = .FALSE.
                !PRINT*, '--------------------- 1'
            ELSE
                DEALLOCATE(td_dra%dat)
                NULLIFY(td_dra%dat)
                ll_allocate = .TRUE.
                !PRINT*, '--------------------- 2'
            END IF
        ELSE
            ll_allocate = .TRUE.
            !PRINT*, '--------------------- 3'
        END IF
        il_stat=0
        IF(ll_allocate.AND.(id_size>0))ALLOCATE( td_dra%dat(il_lb:il_lb+id_size-1))!, STAT = il_stat)
        IF(il_stat/=0)THEN
            PRINT*, 'In reset_drv, ALLOCATE fails, STAT = ',il_stat
            CALL ABORT
        END IF
    END SUBROUTINE reset_drv

    SUBROUTINE assign_div_vector(td_dia, id_data)
        TYPE(dyn_iVector), INTENT(INOUT) :: td_dia
        INTEGER, DIMENSION(:), INTENT(IN) :: id_data

        CALL reset_div(td_dia, SIZE(id_data))
        td_dia%dat = id_data
    END SUBROUTINE assign_div_vector

    SUBROUTINE reset_div(td_dia, id_size, id_lb)
        TYPE(dyn_iVector), INTENT(INOUT) :: td_dia
        INTEGER, INTENT(IN) :: id_size
        INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound in the array if needed
        INTEGER :: il_lb, il_stat
        LOGICAL :: ll_allocate

        IF ( PRESENT(id_lb) ) THEN
            il_lb = id_lb
        ELSE
            il_lb = 1
        END IF

        IF( ASSOCIATED(td_dia%dat) )THEN
            IF( (id_size==SIZE(td_dia%dat)).AND.(il_lb == LBOUND(td_dia%dat, 1) ) )THEN
                ll_allocate = .FALSE.
            ELSE
                DEALLOCATE(td_dia%dat)
                NULLIFY(td_dia%dat)
                ll_allocate = .TRUE.
            END IF
        ELSE
            ll_allocate = .TRUE.
        END IF
        il_stat=0
        IF(ll_allocate.AND.(id_size>0))ALLOCATE( td_dia%dat(il_lb:il_lb+id_size-1))!, STAT = il_stat)
        IF(il_stat/=0)THEN
            PRINT*, 'In reset_div, ALLOCATE fails, STAT = ',il_stat
            CALL ABORT
        END IF
    END SUBROUTINE reset_div

    SUBROUTINE assign_dsv_vector(td_dsa, ada_data)
        TYPE(dyn_sVector), INTENT(INOUT) :: td_dsa
        CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: ada_data
        INTEGER :: ibi, il_size

        il_size = SIZE(ada_data)
        CALL reset_dsv(td_dsa, il_size)
        DO ibi = 1, il_size
            td_dsa%dat(ibi) = ada_data(ibi)
        END DO
    END SUBROUTINE assign_dsv_vector

    FUNCTION dsv_maxLen(td_dsa) RESULT(il_max)
        TYPE(dyn_sVector), INTENT(IN) :: td_dsa
        INTEGER :: il_max, il_tmp, ibi
        il_max = LEN_TRIM( td_dsa%dat( LBOUND(td_dsa%dat,1) ) )
        DO ibi = LBOUND(td_dsa%dat,1)+1, UBOUND(td_dsa%dat,1)
            il_tmp = LEN_TRIM( td_dsa%dat(ibi) )
            IF (il_max < il_tmp) il_max = il_tmp
        END DO
    END FUNCTION dsv_maxLen

    SUBROUTINE reset_dsv(td_dsa, id_size, id_lb)
        TYPE(dyn_sVector), INTENT(INOUT) :: td_dsa
        INTEGER, INTENT(IN) :: id_size
        INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound in the array if needed
        INTEGER :: il_lb, il_stat
        LOGICAL :: ll_allocate

        IF ( PRESENT(id_lb) ) THEN
            il_lb = id_lb
        ELSE
            il_lb = 1
        END IF

        IF( ASSOCIATED(td_dsa%dat) )THEN
            IF( (id_size==SIZE(td_dsa%dat)).AND.(il_lb == LBOUND(td_dsa%dat, 1) ) )THEN
                ll_allocate = .FALSE.
            ELSE
                DEALLOCATE(td_dsa%dat)
                NULLIFY(td_dsa%dat)
                ll_allocate = .TRUE.
            END IF
        ELSE
            ll_allocate = .TRUE.
        END IF
        il_stat = 0
        IF(ll_allocate.AND.(id_size>0))ALLOCATE( td_dsa%dat(il_lb:il_lb+id_size-1))!, STAT = il_stat)
        IF(il_stat/=0)THEN
            PRINT*, 'In reset_dsv, ALLOCATE fails, STAT = ',il_stat
            CALL ABORT
        END IF
    END SUBROUTINE reset_dsv

    !< \brief assign 1D array data to a 2D dynamic vector, the second dimension of the dynamic array is resize to 1
    !! \param[in, out] td_dd 2D dynamic array to wich one wants to assign a 1D array
    !! \param[in] rd_data array data to be assigned
    !! \todo overload assignment operator
    !>
    SUBROUTINE assign_2DdrA_vector(td_dd, rda_data)
        TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
        REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data

        CALL reset_d2DrA(td_dd, 1)
        td_dd%dat(1) = rda_data
    END SUBROUTINE assign_2DdrA_vector

    !< \brief assign 2D array data to a 2D dynamic array
    !! \param[in, out] td_dd 2D dynamic array to wich one wants to assign a 2D array
    !! \param[in] rd_data array data to be assigned
    !! \todo overload assignment operator
    !! \details each element of td_dd receives a column of rda_data, recall that 2D dynamic array is a java style 2D array
    !! column major assumption
    !>
    SUBROUTINE assign_2DdrA_2DA(td_dd, rda_data)
        TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
        REAL(KIND=cp), DIMENSION(:, :), INTENT(IN) :: rda_data
        INTEGER :: ibi

        CALL reset_d2DrA(td_dd, size(rda_data, 2))
        DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
            td_dd%dat(ibi) = rda_data(:, ibi)
        END DO
    END SUBROUTINE assign_2DdrA_2DA

    !< \brief assign 1D array data to an element of a 2D dynamic array
    !! recall that, a 2D dynamic array is like java 2D array
    !! \param[in, out] td_dd 2D dynamic array to wich one wants to put a 1D vector
    !! \param[in] rd_data vector data to be assigned
    !! \param[in] id_idx index of the 1D vector in td_dd, recall that 2D dynamic array is a java style 2D array
    !! \details
    !! \remark call to this routine can causes segmentation fault
    !! if the first dimension of the td_dd is less than id_idx
    !>
    SUBROUTINE put_2DdrA_vector(td_dd, id_idx, rda_data)
        TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
        INTEGER, INTENT(IN)               :: id_idx
        REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data

        td_dd%dat(id_idx) = rda_data
    END SUBROUTINE put_2DdrA_vector

    !< \brief get a 1D array data from a 2D dynamic array
    !! recall that, a 2D dynamic array is like java 2D array
    !! \param[in] td_dd 2D dynamic array from wich one wants to get a 1D vector
    !! \param[out] rda_data vector data to be get from td_dd
    !! \param[in] id_idx index of the 1D vector in td_dd, recall that 2D dynamic array is a java style 2D array
    !! \details
    !! \remark call to this routine can cause segmentation fault
    !! if the first dimension of the td_dd is less than id_idx
    !! also make sure that rl_data has the same size as
    !>
    SUBROUTINE get_2DdrA_vector(td_dd, id_idx, rda_data)
        TYPE(dyn_2DrArray), INTENT(IN) :: td_dd
        INTEGER, INTENT(IN)            :: id_idx
        REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_data

        rda_data = td_dd%dat(id_idx)%dat
    END SUBROUTINE get_2DdrA_vector

    SUBROUTINE reset_d2DrA(td_dd, id_size, id_lb)
        TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
        INTEGER, INTENT(IN) :: id_size
        INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound of the first dimension if needed
        INTEGER :: il_lb, ibi
        IF ( PRESENT(id_lb) ) THEN
            il_lb = id_lb
        ELSE
            il_lb = 1
        END IF

        IF(ASSOCIATED(td_dd%dat)) THEN
            DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
                CALL reset_drv(td_dd%dat(ibi), 0)
            END DO
            DEALLOCATE(td_dd%dat)
            NULLIFY(td_dd%dat)
        END IF
        IF(id_size>0)THEN
            ALLOCATE( td_dd%dat(il_lb:il_lb+id_size-1) )
            DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
                NULLIFY(td_dd%dat(ibi)%dat)
            END DO
        END IF
    END SUBROUTINE reset_d2DrA

    !1D dynamic array minval
    FUNCTION drv_minval(td_dd) RESULT(rl_min)
        TYPE(dyn_rVector), INTENT(IN) :: td_dd
        REAL(cp) :: rl_min

        rl_min = MINVAL(td_dd%dat)
    END FUNCTION drv_minval

    !array of 1D dynamic array minval
    FUNCTION drv_array_minval(tda_dd) RESULT(rl_min)
        TYPE(dyn_rVector), DIMENSION(:), INTENT(IN) :: tda_dd
        REAL(cp) :: rl_min, rl_temp
        INTEGER :: ibi

        rl_min = MINVAL(tda_dd(1)%dat)
        DO ibi = 2, SIZE(tda_dd)
            rl_temp = MINVAL(tda_dd(ibi)%dat)
            IF(rl_min>rl_temp) rl_min = rl_temp
        END DO
    END FUNCTION drv_array_minval

    !2D dynamic array minval
    FUNCTION d2DrA_minval(td_dd) RESULT(rl_min)
        TYPE(dyn_2DrArray), INTENT(IN) :: td_dd
        REAL(cp) :: rl_min, rl_temp
        INTEGER :: ibi, il_lb, il_ub

        il_lb = LBOUND(td_dd%dat,1)
        il_ub = UBOUND(td_dd%dat,1)
        rl_min = MINVAL(td_dd%dat( il_lb )%dat)
        DO ibi = il_lb+1, il_ub
            rl_temp = MINVAL(td_dd%dat(ibi)%dat)
            IF(rl_min>rl_temp) rl_min = rl_temp
        END DO
    END FUNCTION d2DrA_minval

    !1D dynamic array maxval
    FUNCTION drv_maxval(td_dd) RESULT(rl_max)
        TYPE(dyn_rVector), INTENT(IN) :: td_dd
        REAL(cp) :: rl_max

        rl_max = MAXVAL(td_dd%dat)
    END FUNCTION drv_maxval

    !array of 1D dynamic array maxval
    FUNCTION drv_array_maxval(tda_dd) RESULT(rl_max)
        TYPE(dyn_rVector), DIMENSION(:), INTENT(IN) :: tda_dd
        REAL(cp) :: rl_max, rl_temp
        INTEGER :: ibi

        rl_max = MAXVAL(tda_dd(1)%dat)
        DO ibi = 2, SIZE(tda_dd)
            rl_temp = MAXVAL(tda_dd(ibi)%dat)
            IF(rl_max<rl_temp) rl_max = rl_temp
        END DO
    END FUNCTION drv_array_maxval

    !< 2D dynamic array maxval
    FUNCTION d2DrA_maxval(td_dd) RESULT(rl_max)
        TYPE(dyn_2DrArray), INTENT(IN) :: td_dd
        REAL(cp) :: rl_max, rl_temp
        INTEGER :: ibi, il_lb, il_ub

        il_lb = LBOUND(td_dd%dat,1)
        il_ub = UBOUND(td_dd%dat,1)
        rl_max = MAXVAL(td_dd%dat( il_lb )%dat)
        DO ibi = il_lb+1, il_ub
            rl_temp = MAXVAL(td_dd%dat(ibi)%dat)
            IF(rl_max<rl_temp) rl_max = rl_temp
        END DO
    END FUNCTION d2DrA_maxval

    !> \brief compute the number of regularly-spaced coordinates
    !! \param[in] ida_max maximum coordinates that can be generate in each dimension
    !! \param[in] shifts optional minimum coordinates that can be generate in each dimension, default is 1 for each dimension. They give the lowest coordinates.
    !! \param[in] steps step between successive coordinates in each dimension
    !! \param[out] ncoords optional coordinates count in each dimension
    !! \detail Asumming that dim_count is the number of dimensions of the problem, steps, ida_max [, shifts] [and ncoords]  are one-dimensional arrays of size dim_count each. This subroutine also compute the number of regularly-spaced coordinates in each dimension of the problem, see the optional argument ncoords
    !! \todo add a version for real coordinates and use the overloading, rename this as regular_int_coord_count
    !<
    FUNCTION regular_coord_count(ida_max, shifts, steps, ncoords) RESULT(il_count)
      INTEGER, DIMENSION(:), INTENT(IN) :: ida_max
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)  :: shifts, steps
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: ncoords
      INTEGER, DIMENSION(SIZE(ida_max)) :: ila_shift, ila_step, ila_count
      INTEGER :: il_count

      IF( PRESENT(shifts) ) THEN
        ila_shift = shifts
      ELSE
        ila_shift = 1
      END IF
      IF( PRESENT(steps) ) THEN
        ila_step = steps
      ELSE
        ila_step = 1
      END IF
      ila_count = 1 + (ida_max - ila_shift)/ila_step
      il_count = PRODUCT(ila_count)
      IF( PRESENT(ncoords) ) ncoords = ila_count
    END FUNCTION regular_coord_count

    !> \brief build regularly-spaced coordinates of type integer
    !! \param[out] ida_coord regularly-spaced coordinates generated by the subroutine
    !! \param[in] ncoords optional array giving the number of coordinates in each dimension
    !! \param[in] max_coords optional array giving the maximum coordinate in each dimension
    !! \param[in] shifts optional array giving the shift from zero in each dimension. Default is 1 to conform to one-base indexation in fortran.
    !! \param[in] steps optional array giving the step between successive coordinates in each dimension. Default is 1.
    !! \detail At least one of ncoords and max_coords must be provided. If ncoords is provided, max_coords is not necessary. Asumming that coord_count is the total number of coordinates to generate and dim_count is the number of dimensions of the problem, ida_coord is of shape (dim_count, coord_count). The second dimension can be larger than coord_count but not less. ida_step, ida_max and ida_min (if provided) are one-dimensional array of size dim_count.
    !! This routine is used to build regularly-spaced coordinates for twin observation. The coordinates are 1-based indexed for fortran default indexation. The integer coordinates can easilly be converted to real-value coordinates if uniformly-spaced mesh is used. Shift to 0-based index and multiply by the space step (each dimension separately if the space step if different for each), then add the real-coordinates of the origin if not 0.
    !<
    SUBROUTINE build_regular_int_coord(ida_coord, ncoords, max_coords, shifts, steps)
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: ida_coord
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)  :: max_coords, shifts, steps, ncoords
      INTEGER, DIMENSION(SIZE(ida_coord,1)) :: ila_shift, ila_step, ila_ncoord, ila_prodinf
      INTEGER :: il_ncoord, ibi, il_idx, il_ndim, il_val

      IF( PRESENT(shifts) ) THEN
        ila_shift = shifts
      ELSE
        ila_shift = 1
      END IF
      IF( PRESENT(steps) ) THEN
        ila_step = steps
      ELSE
        ila_step = 1
      END IF
      IF( PRESENT(ncoords) ) THEN
        ila_ncoord = ncoords
      ELSE IF ( PRESENT(max_coords) )THEN
        PRINT*, 'Calling regular_coord_count'
        il_ncoord = regular_coord_count(max_coords, ila_shift, ila_step, ila_ncoord)
      END IF
      PRINT*, 'checking PRESENT(ncoords).OR.PRESENT(max_coords)'
      IF( PRESENT(ncoords).OR.PRESENT(max_coords) ) THEN
        il_ndim = SIZE(ida_coord,1)
        il_ncoord = PRODUCT(ila_ncoord)
        !the index i of  variable ila_prodinf stores the product of the (i-1) first elts of il_ncoord, ila_obs actually give the shape of the non vectorized obs array
        ila_prodinf(1) = 1
        DO ibi = 2,il_ndim
          ila_prodinf(ibi) = ila_prodinf(ibi-1)*ila_ncoord(ibi-1)
        END DO
        PRINT*, 'starting loop'
        PRINT*, 'ila_ncoord   = ', ila_ncoord
        PRINT*, 'ila_step     = ', ila_step
        PRINT*, 'ila_shift    = ', ila_shift
        PRINT*, 'ila_prodinf  = ', ila_prodinf
        DO il_idx = 1,il_ncoord
          PRINT*, 'il_idx = ', il_idx
          il_val = il_idx-1 !compute zero-based coordinates
          DO ibi = il_ndim,2,-1
            PRINT*, '        ibi = ', ibi
            PRINT*, 'SHAPE(ida_coord) = ', SHAPE(ida_coord)

            !shift to one-based coordinated with appropriate shift
            ida_coord(ibi, il_idx) = il_val/ila_prodinf(ibi)*ila_step(ibi) + ila_shift(ibi)
            PRINT*, 'here  '
            il_val = mod( il_val, ila_prodinf(ibi) )
            PRINT*, 'ending nested loop = '
          END DO
          ida_coord(1, il_idx) = il_val*ila_step(ibi) + ila_shift(ibi)
            PRINT*, 'ending outside loop = '
        END DO
      ELSE
        ida_coord = -999
      END IF
    END SUBROUTINE build_regular_int_coord

    !> \brief Builds regularly spaced indices
    !! \param[in, out] ida_idx array of builded indices
    !! \param[in] id_idxMin minimal value for idx
    !! \param[in] id_idxMax maximal value for idx
    !! \param[in] id_obsDist distance between two observations
    !! \remark: this subroutine was previously named build_regular_idx
    !>
    SUBROUTINE build_regular_idx(ida_idx, id_idxMin, id_idxMax, id_obsDist)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ida_idx
        INTEGER, INTENT(IN) :: id_idxMin, id_idxMax, id_obsDist
        INTEGER :: ibi, ibj

        ibj = 1
        DO ibi = id_idxMin, id_idxMax, id_obsDist
            ida_idx(ibj) = ibi
            ibj = ibj + 1
        END DO
    END SUBROUTINE build_regular_idx

    !> \brief Save 1D trajectory in ascii file
    !! \param[in] tda_trj the trajectory to be saved
    !! \param[in] ada_fileName, name of the file in wich date should be saved
    !! \param[in] id_lb (optional) lower bound of the interresting part of data to be saved at each time step
    !! \param[in] id_ub (optional) upper bound of the interresting part of data to be saved at each time step
    !! \param[in] ld_interpolate (optional) says if interpolation is needed. This is the case if one need all data to be saved at the same location
    !! \param[in] ld_column (optional) says if each element of the trajectory should be considered as a column(usefull for gnuplot)
    !! \details Data are saved time step by time step : old the date from the first time step, then the second time step and so on.
    !! parameters id_lb, id_ub are useful in the case where one do not want to save everything, for exemple, the extension of the state variable for boundary condition
    !! Interpolation suppose that the variable is computed between grid point, the goal is to interpolate to grid point before saving. In this case, id_lb, id_ub are the bound of the computed variable.
    !<
    SUBROUTINE save_trj(tda_trj, ada_fileName, id_lb, id_ub, ld_interpolate, ld_column)
      INTEGER, PARAMETER :: ip_fid = 41
      TYPE(dyn_rVector), DIMENSION(:), INTENT(IN) :: tda_trj
      CHARACTER(LEN=*), INTENT(IN) :: ada_fileName
      INTEGER, OPTIONAL, INTENT(IN) :: id_lb, id_ub
      LOGICAL, OPTIONAL, INTENT(IN) :: ld_interpolate
      LOGICAL, OPTIONAL, INTENT(IN) :: ld_column
      !local variables
      REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: rla_tmp
      INTEGER :: ibk, ibj, il_lb, il_ub, il_ios
      LOGICAL :: ll_interpolate, ll_column

      !Checking optional argument
      IF( PRESENT(id_lb) ) THEN
         il_lb = id_lb
      ELSE
         il_lb = LBOUND( tda_trj(1)%dat, 1)
      END IF
      IF( PRESENT(id_ub) ) THEN
         il_ub = id_ub
      ELSE
         il_ub = UBOUND( tda_trj(1)%dat, 1)
      END IF
      IF( PRESENT(ld_interpolate) ) THEN
         ll_interpolate = ld_interpolate
      ELSE
         ll_interpolate = .FALSE.
      END IF
      IF( PRESENT(ld_column) ) THEN
         ll_column = ld_column
      ELSE
         ll_column = .FALSE.
      END IF

      OPEN( UNIT = ip_fid, FILE = ada_fileName, STATUS = FILE_REPLACE, FORM = FILE_FORMATTED, IOSTAT = il_ios)
      IF(il_ios/=0) THEN
        WRITE(*, *) 'In save_uTrj :: error creating the file ', ada_fileName
        STOP
      END IF

      IF (ll_column) THEN
        ALLOCATE( rla_tmp( SIZE(tda_trj) ) )
        !save the first row
        DO ibj =1, SIZE(rla_tmp)
          rla_tmp(ibj) = tda_trj(ibj)%dat(il_lb)
        END DO
        WRITE(UNIT = ip_fid, FMT=*) rla_tmp
        !save intermediate rows
        DO ibk = il_lb+1, il_ub
          DO ibj =1, SIZE(rla_tmp)
            IF( ll_interpolate ) THEN
              rla_tmp(ibj) = (tda_trj(ibj)%dat(ibk-1) + tda_trj(ibj)%dat(ibk) ) / 2.0_cp !interpolation
            ELSE
              rla_tmp(ibj) = tda_trj(ibj)%dat(ibk)
            END IF
          END DO
          WRITE(UNIT = ip_fid, FMT=*) rla_tmp
        END DO
        !save the last row
        IF( ll_interpolate ) THEN
          DO ibj =1, SIZE(rla_tmp)
            rla_tmp(ibj) = tda_trj(ibj)%dat(il_ub)
          END DO
          WRITE(UNIT = ip_fid, FMT=*) rla_tmp
        END IF
      ELSE
        IF( ll_interpolate ) THEN
          ALLOCATE( rla_tmp(il_lb:il_ub+1) )
        ELSE
          ALLOCATE( rla_tmp(il_lb:il_ub) )
        END IF

        DO ibk = 1, SIZE(tda_trj)
          IF (ll_interpolate) THEN ! interpolation
              rla_tmp(il_lb) = tda_trj(ibk)%dat(il_lb)
              rla_tmp(il_lb+1:il_ub) = ( tda_trj(ibk)%dat(il_lb:il_ub-1) + tda_trj(ibk)%dat(il_lb+1:il_ub) ) / 2.0_cp
              rla_tmp(il_ub+1) = tda_trj(ibk)%dat(il_ub)
          ELSE
              rla_tmp = tda_trj(ibk)%dat(il_lb:il_ub)
          END IF
          WRITE(UNIT = ip_fid, FMT=*) rla_tmp
        END DO
      END IF

      CLOSE(ip_fid)
      PRINT*, 'SIZE(rla_tmp) = ', SIZE(rla_tmp)
      DEALLOCATE(rla_tmp)
    END SUBROUTINE save_trj




  ! Compute the value of Whittaker function W_{(\alpha-1)/2, 1/2}
  ! For the cases alpha = 1 or alpha = 3 only
  ! /!\ Need to be completed by a general formula for all alpha
  ! see Gradshteyn, I.S., Ryshik, I.M. 1980, Tables of Integrals, Academic Press
  FUNCTION whittaker( r, rd_rvm )
    !!CALL as xx = whittaker (r, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: r, rd_rvm
    REAL(KIND=dp)               :: whittaker
    REAL(KIND=dp)               :: w1, w3
    ! (alpha = 1) non-isolated vortex (gaussian velocity profile)
    ! Warning : w1 does not reach its min and max values at r = pprvm !
    IF ( ABS(r) <= 1d-16 ) THEN
        w1 = 0.0_dp
    ELSE
        w1 = rd_rvm * rd_rvm * (1.0 - exp( -0.5_dp * ( r / rd_rvm ) ** 2 ) ) / r
    END IF

    ! (alpha = 3) isolated vortex (gaussian h profile)
    IF( ABS(r) <= 1d-16 ) THEN
        w3 = 0.0_dp
    ELSE
        w3 = rd_rvm * exp( 0.5_dp) * ( r / rd_rvm ) * exp( -0.5_dp * ( r / rd_rvm )** 2 )
    END IF

    ! Choose either w1 or w3 depending if the vortex is  isolated (w3) or non-isolated (w1)
    whittaker = w3

  END FUNCTION whittaker

  ! value of the y-coordinate of the initial velocity of the vortex at (x,y)
  ! Uses whittaker function
  FUNCTION u_init( x, y, rd_xvc, rd_yvc, rd_rvm)
    !!Call as xx = u_init(x, y, tg_swp%r_xvc, tg_swp%r_yvc, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: x, y, rd_xvc, rd_yvc, rd_rvm
    REAL(KIND=dp)               :: r
    REAL(KIND=dp)               :: u_init

    r = SQRT( ( x - rd_xvc ) ** 2 + ( y - rd_yvc ) ** 2 )

    IF( ABS(r) < 1d-16 ) THEN
        u_init = 0.0_dp
    ELSE
        u_init =  - ( y - rd_yvc ) * whittaker( r, rd_rvm ) / r
    END IF

    !* Simulation de l'initialisation par agitation dans un cylindre
    !* que l'on retire par la suite
    !IF( ( ABS(r) < 1E-16 ) .OR. ( ABS(r) > 5 * rd_rvm ) ) THEN
    !   u_init = 0.0
    !ELSE
    !   u_init =  - ( y - rd_yvc ) * whittaker( r ) / r
    !END IF

  END FUNCTION u_init

  ! value of the y-coordinate of the initial velocity of the vortex at (x,y)
  ! Uses whittaker function
  FUNCTION v_init( x, y, rd_xvc, rd_yvc, rd_rvm )
    !!Call as xx = v_init(x, y, tg_swp%r_xvc, tg_swp%r_yvc, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: x, y, rd_xvc, rd_yvc, rd_rvm
    REAL(KIND=dp)               :: r
    REAL(KIND=dp)               :: v_init

    r = SQRT( ( x - rd_xvc ) ** 2 + ( y - rd_yvc ) ** 2 )

    IF ( ABS(r) <= 1d-16 ) THEN
        v_init = 0.0_dp
    ELSE
        v_init =  ( x - rd_xvc ) * whittaker( r, rd_rvm ) / r
    END IF

    !* Simulation de l'initialisation par agitation dans un cylindre
    !* que l'on retire par la suite
    ! IF( ( ABS(r) <= 1E-16 ) .OR. ( ABS(r) > 5 * rd_rvm ) ) THEN
    !    v_init = 0.0
    ! ELSE
    !    v_init =  ( x - rd_xvc ) * whittaker( r ) / r
    ! END IF

  END FUNCTION v_init

!   ! intialialize the velocity at the node of the C-grid
!   ! u and v are matrices
!   ! Uses u_init and v_init functions
!   SUBROUTINE init_velocity( td_field, rd_xvc, rd_yvc, rd_rvm)
!     !!Call as init_velocity( tg_field, tg_swp%r_xvc, tg_swp%r_yvc, tg_swp%r_rvm)
!     IMPLICIT NONE
!     TYPE(modstate), INTENT( INOUT ) :: td_field
!     REAL(KIND=dp), INTENT( IN ) :: rd_xvc, rd_yvc, rd_rvm
!     REAL(KIND=dp)     :: xi
!     INTEGER    :: ibi,ibj
!     ! PRINT*, "Profil de vitesse initiale"
!     ! x-coordinate of velocity
!     !PRINT*, 'In init_velocity, rd_rvm = ', rd_rvm
!     DO ibi = 1, td_field%i
!
!         xi = ( ibi ) * tg_gp%dx
!
!         DO ibj = 0, td_field%j
!
!           td_field%u(ibi,ibj) = u_init( xi, ( ibj ) * tg_gp%dy + tg_gp%dy / 2.0_dp, rd_xvc, rd_yvc, rd_rvm )
!
!         END DO
!
!     END DO
!
!     ! y-coordinate of velocity
!     DO ibi = 0, td_field%i
!
!         xi = ( ibi ) * tg_gp%dx + tg_gp%dx / 2.0_dp
!
!         DO ibj = 1, td_field%j
!
!           td_field%v(ibi,ibj) = v_init( xi, ( ibj ) * tg_gp%dy, rd_xvc, rd_yvc, rd_rvm )
!
!         END DO
!
!     END DO
!
!   END SUBROUTINE init_velocity

END MODULE tools
