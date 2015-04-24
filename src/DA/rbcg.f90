!> @file rbcg.f90
!! Defines the Restricted B-preconditioned style conjugate gradient
!<

!> @brief Restricted B-preconditioned style conjugate gradient(RBCG) module
!!
!<
module rbcg
    use debug_tools
    use com_tools
implicit none

contains

    !> @brief Solve the Ensemble Kalman Filter system for one member using RBCG
    !! @param [in] td_ep exchange parameter data structure
    !! @param [in] obs_std vector of standard deviation of observation error
    !! @param [in] d right hand side
    !! @param [in,out] x solution of the system
    !! @param [in] matrixmult procedure that does the matric multiplication
    !! @param [in] nItrMax maximum number of iterations
    !! @param [in] tol stopping criteria based on the B-norm of the gradient
    !!
    !! @details
    !! This algorithm solve a system of the form (A+R)x = d
    !! where A and R (obs_std) are symmetric positive (semi-)definite matrices
    !! This implementation considers R to be diagonal
    !! @note
    !!   for now, there is no size or shape check on the array variables
    !<
    subroutine rbcg_solve(td_ep, obs_std, d, x, matrixmult, nItrMax, tol)
        interface
            subroutine matrixmult(ep, r,w,diverged)
                use general_constant
                use com_tools
                type(exchange_param), intent(inout) :: ep
                real(dp), dimension(:), intent(in) :: r
                real(dp), dimension(:), intent(inout) :: w
                logical, intent(out) :: diverged
            end subroutine matrixmult
        end interface
        type(exchange_param), intent(inout) :: td_ep
        real(dp), dimension(:), intent(in) :: obs_std
        real(dp), dimension(:), intent(in) :: d
        real(dp), dimension(:), intent(in out) :: x
        integer, intent(in) :: nItrMax
        real(dp), intent(in) :: tol
        ! local variables
        ! RBCG variables
        real(dp) :: cg_r(size(d))
        real(dp) :: cg_lambda(size(d))
        real(dp) :: cg_t(size(d))
        real(dp) :: cg_p(size(d))
        real(dp) :: cg_w(size(d))
        real(dp) :: cg_q(size(d))
        real(dp) :: cg_alpha
        real(dp) :: cg_beta
        ! others
        real(dp) :: r0_norm ! B-norm of the residue at the zeroth iteration
        real(dp) :: ri_norm ! B-norm of the residue at the current iteration
        real(dp) :: cg_wTr ! dot product w by r (i)
        real(dp) :: cg_qTt ! dot product q by t
        real(dp) :: cg_wTrp1 ! dot product w by r (i+1)
        real(dp) :: cg_omega, cg_epsilon
        integer :: cg_itr
        logical :: converged, hv_diverged
        character(*), parameter :: fmt0= "('*** RBCG   Iter ***',I3,' -- |r_0|=',ES14.5E3)"
        character(*), parameter :: fmti= "('*** RBCG   Iter ***',I3,' -- |r|/|r_0|=',ES14.5E3)"
        !
        call debug(tol, "RBCG - tol", tag=dALLWAYS)
        call debug(nItrMax, "RBCG - nItrMax", tag=dALLWAYS)
        cg_epsilon = epsilon(1.0_dp)
        cg_itr = 0
        cg_lambda = 0.0
        cg_r = d/obs_std**2
        cg_p = cg_r
        ! Computing w = HA(HA)^Tr
        call matrixmult(td_ep, cg_r, cg_w, hv_diverged)
        cg_t = cg_w

        cg_wTr = dot_product(cg_w, cg_r)
        ri_norm = sqrt(cg_wTr)
        r0_norm    = ri_norm
        cg_omega= ri_norm/r0_norm
        write(*,fmt=fmt0)cg_itr,r0_norm
        cg_itr = 1
        converged = .false.

        do while( (cg_itr<= nItrMax).and.(.not.converged) )
            cg_q = cg_t/obs_std**2 + cg_p
            ! The following two quantities are used as denominator
            ! So if one is zero, either, there is convergence
            !  or something is wrong.
            ! In addition, cg_wTr is the result of an inner product,
            ! so it must be greater than zero.
            ! Below we check is abs(cg_qTt) is two small
            ! and only check if cg_wTr is two small (it is already positive)
            cg_wTr = dot_product(cg_w, cg_r)
            cg_qTt = dot_product(cg_q, cg_t)
            if( (abs(cg_qTt)>cg_epsilon).and.(cg_wTr>cg_epsilon) )then
                cg_alpha = cg_wTr/cg_qTt
                cg_lambda = cg_lambda + cg_alpha*cg_p
                cg_r = cg_r - cg_alpha*cg_q
                !re-conjugate r
                call matrixmult(td_ep, cg_r, cg_w, hv_diverged)
                cg_wTrp1 = dot_product(cg_w, cg_r)
                cg_beta = cg_wTrp1/cg_wTr
                if( abs(cg_beta)<cg_epsilon )then
                    converged = .true.
                end if
                cg_p = cg_r + cg_beta*cg_p
                cg_t = cg_w + cg_beta*cg_t
                ri_norm = sqrt(cg_wTrp1)
                cg_omega= ri_norm/r0_norm
                call debug(cg_omega, "In loop, cg_omega", tag=dALLWAYS)

                !Check the convergence
                if (cg_omega > tol) then
                    write(*,fmt=fmti)cg_itr,cg_omega
                else
                    converged = .true.
                endif
                !
                !Moving to the next iteration
                !
                cg_itr=cg_itr+1
            else
                converged = .true. ! the denominator of alpha is zero or the denominator of beta is zero
            end if
        end do

        if( .not.converged ) then
            call debug( nItrMax, "Maximum number of iterations reached", tag=dALLWAYS )
        else
            call debug(cg_itr, "RBCG converged, total nb iterations", tag=dALLWAYS)
            call debug(r0_norm, "  r0_norm", tag=dALLWAYS)
            call debug(cg_wTrp1, "  cg_wTrp1", tag=dALLWAYS)
            call debug(cg_wTr, "  cg_wTr", tag=dALLWAYS)
            call debug(cg_qTt, "  cg_qTt", tag=dALLWAYS)
            call debug(cg_beta, "  cg_beta", tag=dALLWAYS)
            call debug(cg_omega, "RBCG converged, cg_omega", tag=dALLWAYS)
            call debug(sum(abs(cg_lambda)), "Solution sum(abs(cg_lambda))")
            if(abs(cg_qTt)<=cg_epsilon)call debug(cg_qTt, "q*t is too small", tag=dALLWAYS)
            if(abs(cg_wTr)<=cg_epsilon)call debug(cg_wTr, "w*r is too small", tag=dALLWAYS)
            if(abs(cg_beta)<=cg_epsilon)call debug(cg_beta, "beta is too small", tag=dALLWAYS)
        end if

        x = cg_lambda
    end subroutine rbcg_solve



    !> @brief Solve the linear system using Conjugate gradient
    !! @param [in] td_ep exchange parameter data structure
    !! @param [in] obs_std vector of standard deviation of observation error
    !! @param [in] d right hand side
    !! @param [in,out] x solution of the system
    !! @param [in] matrixmult procedure that does the matric multiplication
    !! @param [in] nItrMax maximum number of iterations
    !! @param [in] tol stopping criteria based on the B-norm of the gradient
    !!
    !! @details
    !! This algorithm solve a system of the form (A+R)x = d
    !! where A and R (obs_std) are symmetric positive (semi-)definite matrices
    !! This implementation considers R to be diagonal
    !! @note
    !!   for now, there is no size or shape check on the array variables
    !<
    subroutine cg_solve(td_ep, obs_std, d, x, matrixmult, nItrMax, tol)
        interface
            subroutine matrixmult(ep, r,w,diverged)
                use general_constant
                use com_tools
                type(exchange_param), intent(inout) :: ep
                real(dp), dimension(:), intent(in) :: r
                real(dp), dimension(:), intent(inout) :: w
                logical, intent(out) :: diverged
            end subroutine matrixmult
        end interface
        type(exchange_param), intent(inout) :: td_ep
        real(dp), dimension(:), intent(in) :: obs_std
        real(dp), dimension(:), intent(in) :: d
        real(dp), dimension(:), intent(in out) :: x
        integer, intent(in) :: nItrMax
        real(dp), intent(in) :: tol
        ! local variables
        ! RBCG variables
        real(dp) :: cg_r(size(d))
        real(dp) :: cg_lambda(size(d))
        real(dp) :: cg_p(size(d))
        real(dp) :: cg_q(size(d))
        real(dp) :: cg_alpha
        real(dp) :: cg_beta
        ! others
        real(dp) :: r0_norm ! norm of the residue at the zeroth iteration
        real(dp) :: ri_norm ! norm of the residue at the current iteration
        real(dp) :: cg_rTr ! dot product r by r (i)
        real(dp) :: cg_qTp ! dot product q by p
        real(dp) :: cg_rTrp1 ! dot product r by r (i+1)
        real(dp) :: cg_omega, cg_epsilon
        integer :: cg_itr
        logical :: converged, hv_diverged
        character(*), parameter :: fmt0= "('*** CG   Iter ***',I3,' -- |r_0|=',ES14.5E3)"
        character(*), parameter :: fmti= "('*** CG   Iter ***',I3,' -- |r|/|r_0|=',ES14.5E3)"
        !
        call debug(tol, "RBCG - tol", tag=dALLWAYS)
        call debug(nItrMax, "RBCG - nItrMax", tag=dALLWAYS)
        cg_epsilon = epsilon(1.0_dp)
        cg_itr = 0
        cg_lambda = 0.0
        cg_r = d
        cg_p = cg_r

        cg_rTr = dot_product(cg_r, cg_r)
        ri_norm = sqrt(cg_rTr)
        r0_norm    = ri_norm
        cg_omega= ri_norm/r0_norm
        write(*,fmt=fmt0)cg_itr,r0_norm
        cg_itr = 1
        converged = .false.

        do while( (cg_itr<= nItrMax).and.(.not.converged) )
            call matrixmult(td_ep, cg_p, cg_q, hv_diverged)
            cg_q = cg_q + obs_std*cg_p
            ! The following two quantities are used as denominator
            ! So if one is zero, either, there is convergence
            !  or something is wrong.
            ! In addition, cg_rTr is the result of an inner product,
            ! so it must be greater than zero.
            ! Below we check is abs(cg_qTr) is two small
            ! and only check if cg_wTr is two small (it is already positive)
            cg_rTr = dot_product(cg_r, cg_r)
            cg_qTp = dot_product(cg_q, cg_p)
            if( (abs(cg_qTp)>cg_epsilon).and.(cg_rTr>cg_epsilon) )then
                cg_alpha = cg_rTr/cg_qTp
                cg_lambda = cg_lambda + cg_alpha*cg_p
                cg_r = cg_r - cg_alpha*cg_q
                !re-conjugate r
                cg_rTrp1 = dot_product(cg_r, cg_r)
                cg_beta = cg_rTrp1/cg_rTr
                if( abs(cg_beta)<cg_epsilon )then
                    converged = .true.
                end if
                cg_p = cg_r + cg_beta*cg_p
                ri_norm = sqrt(cg_rTrp1)
                cg_omega= ri_norm/r0_norm

                !Check the convergence
                if (cg_omega > tol) then
                    write(*,fmt=fmti)cg_itr,cg_omega
                else
                    converged = .true.
                endif
                !
                !Moving to the next iteration
                !
                cg_itr=cg_itr+1
            else
                converged = .true. ! the denominator of alpha is zero or the denominator of beta is zero
            end if
        end do

        if( .not.converged ) then
            call debug( nItrMax, "Maximum number of iterations reached", tag=dALLWAYS )
        else
            call debug(cg_itr, "CG converged, total nb iterations", tag=dALLWAYS)
            call debug(r0_norm, "  r0_norm", tag=dALLWAYS)
            call debug(cg_rTrp1, "  cg_wTrp1", tag=dALLWAYS)
            call debug(cg_rTr, "  cg_wTr", tag=dALLWAYS)
            call debug(cg_qTp, "  cg_qTt", tag=dALLWAYS)
            call debug(cg_beta, "  cg_beta", tag=dALLWAYS)
            call debug(cg_omega, "CG converged, cg_omega", tag=dALLWAYS)
            call debug(sum(abs(cg_lambda)), "Solution sum(abs(cg_lambda))")
            if(abs(cg_qTp)<=cg_epsilon)call debug(cg_qTp, "q*p is too small", tag=dALLWAYS)
            if(abs(cg_rTr)<=cg_epsilon)call debug(cg_rTr, "r*r is too small", tag=dALLWAYS)
            if(abs(cg_beta)<=cg_epsilon)call debug(cg_beta, "beta is too small", tag=dALLWAYS)
        end if

        x = cg_lambda
    end subroutine cg_solve
end module rbcg