!> \file svd_EnKF.f90
!! Ensemble Kalman FIlter analysis using SVD decomposition
!! This module implements the Ensemble Kalman Filter as decribed
!! in Geir Evenson, The Ensemble Kalman Filter: theoretical formulation
!! and practical implementation. Ocean Dynamics (2003) 53: 343-367
!! DOI 10.1007/s10236-003-9
!! The paper is refered to as Evenson2003 in the module.
!<

!>\brief ensemble Kalman filter analysis
module svd_EnKF
  use debug_tools
  use general_tools
implicit none
!
!     interface EnKF_X4_analysis
!         module procedure EnKF_analysis_using_X4
!     end interface EnKF_X4_analysis
!
!     interface EnKF_X5_analysis
!         module procedure EnKF_analysis_using_X5
!     end interface EnKF_X5_analysis

contains

    !> @brief Computes the analysed ensemble in the EnKF, using SVD
    !! @param [in,out] X ensemble matrix
    !! @param [in] X_dev ensemble deviation from the ensemble mean
    !! @param [in] HX_dev observation operator applied to Xdev
    !! @param [in] Y_dev ensemble of observation perturbation (each column)
    !! @param [in] D ensemble innovations (each column is an innovation)
    !! @param [in,out] E ensemble model error
    !! @param [in] E_dev ensemble deviation from model error mean
    !! @param [in] incr boolean flag saying if increment is computed or not
    !!   The increment is the update to be applied to the background to get
    !!   the analysis. By default the increment is not computed. If computed,
    !!   the increment is returned in X_dex, resp. E_dev.
    !! @param [in, out] X4 X4 matrix
    !!
    !! @details
    !!   the number of rows of X defines the size of the state variable
    !!   the number of columns of X defines the size of the ensemble
    !!   \a E and \a E_dev are used in conjunction with state augmentation.
    !!   They have the same shape as X. If ensemble augmentation with model error
    !!   is required, both should be present, otherwise, both should be absent
    !<
    subroutine svd_EnKFAnalysis(X, X_dev, HX_dev, Y_dev, D, E, E_dev, incr, X4)
        real(dp), dimension(:,:), intent(inout) :: X
        real(dp), dimension(:,:), intent(inout) :: X_dev
        real(dp), dimension(:,:), intent(in) :: HX_dev, Y_dev, D
        real(dp), dimension(:,:), optional, intent(inout) :: E
        real(dp), dimension(:,:), optional, intent(inout) :: E_dev
        real(dp), dimension(:,:), optional, intent(inout) :: X4
        logical, optional, intent(in) :: incr
        !local variables
        real(dp), dimension( size(X,2), size(X,2) ) :: T4
        real(dp), dimension(:,:), allocatable ::  X_inc, E_inc
        !dimension variables
        integer :: ndim !! Dimension of model state
        integer :: nmem !! Number of ensemble members
        integer :: edim !! Dimension of model error
        !dgemm factor
        real(dp):: alpha!, beta
        logical :: ll_incr


        !initialization
        ndim = SIZE(X,1)
        nmem = SIZE(X,2)
        if( present(E) )then
            edim = size(E,1)
        else
            edim = -1
        end if
        if( present(incr) )then
            ll_incr = incr
        else
            ll_incr = .false.
        end if

        call debug(sum((X_dev)), 'svd_EnKFAnalysis; sum((X_dev)) = ')
        call svd_EnKF_X4(HX_dev, Y_dev, D, T4)
        call EnKF_analysis_using_X4(T4, X, X_dev, E, E_dev, incr)

!         !Updating X, X = X + X_dev*T4
!         alpha = 1.0_dp
!         call debug('updating X, X = X + X_dev*T4')
!         call dgemm('N', 'N', ndim, nmem, nmem, alpha, X_dev, ndim&
!             , T4, nmem, 1.0_dp, X, ndim)
!         !computing the increment if requested
!         if(ll_incr)then
!             allocate( X_inc(size(X,1),size(X,2) ) )
!             call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X_dev, ndim&
!                         , T4, nmem, 0.0_dp, X_inc, ndim)
!             X_dev = X_inc
!             deallocate(X_inc)
!         end if
!
!         !Updating the model error in the case of state augmentation
!         if( present(E) )then
!             call debug('updating E, E = E + E_dev*T4')
!             call dgemm('N', 'N', edim, nmem, nmem, alpha, E_dev, edim&
!                 , T4, nmem, 1.0_dp, E, edim)
!             !computing the increment if requested
!             if(ll_incr)then
!                 allocate( E_inc(size(E,1),size(E,2) ) )
!                 call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E_dev, edim&
!                 , T4, nmem, 0.0_dp, E_inc, edim)
!                 E_dev = E_inc
!                 deallocate(E_inc)
!             end if
!         end if

        if(present(X4)) X4 = T4
    end subroutine svd_EnKFAnalysis

    !> @brief Compute the matrix X4 eq. 65 [Evenson 2002]
    !! @param [in] HX_dev observation operator applied to Xdev
    !! @param [in] Y_dev ensemble of observation perturbation (each column)
    !! @param [in] D ensemble innovations (each column is an innovation)
    !! @param [in, out] T4 X4 matrix (named T4 for local purpose)
    !!
    !! @details
    !!   the number of columns of HX_dev defines the size of the ensemble
    !!   T4 must be a square matrix with the number of columns/rows given
    !!   by the number of the ensemble members
    !! @note for now there is no check on the shape of matrix parameters
    !<
    subroutine svd_EnKF_X4(HX_dev, Y_dev, D, T4)
        real(dp), dimension(:,:), intent(in) :: HX_dev, Y_dev, D
        real(dp), dimension(:,:), intent(in out) :: T4
        !local variables
        real(dp), dimension(:,:), allocatable :: T, U, VT, T1, T2, T3
        real(dp), dimension(:), allocatable :: sig, sig_inv
        !SVD variables
        integer, dimension(:), allocatable :: IWORK
        real(dp), dimension(:), allocatable :: WORK
        real(dp) :: sigsum, sigsum1
        integer :: nrsigma
        !dimension variables
        integer :: nmem !! Number of ensemble members
        integer :: edim !! Dimension of model error
        integer :: nobs !! Number of observations
        integer :: minMR !!min(nmem, nobs)
        integer :: lwork, info, i

        !initialization
        nmem = SIZE(HX_dev,2)
        nobs = SIZE(D,1)
        minMR = min(nmem, nobs)

        !debugging
        call debug(nobs, 'nobs = ')
        call debug(sum(abs(HX_dev)), 'svd_EnKF_X4; sum(abs(HX_dev))')
        call debug(sum(abs(Y_dev)) , 'svd_EnKF_X4; sum(abs(Y_dev)).')
        call debug(sum(abs(D))     , 'svd_EnKF_X4; sum(abs(D)).....')
        !end of debugging

        !computing T = HX_dev + Y_dev
        allocate( T(nobs, nmem) )
        T = HX_dev + Y_dev

        !computing the SVD of T
        allocate( U(nobs,minMR), VT(minMR, nmem), sig(minMR)&
            , sig_inv(minMR), iwork(8*minMR), work(1) )
        U = 0.0_dp
        VT = 0.0_dp
        sig = 0.0_dp
        lwork = -1 !! query the optimal size of WORK
        call dgesdd('S', nobs, nmem, T, nobs, sig, U, nobs, VT&
        , minMR, work, lwork, iwork, info )
        lwork = int( work(1) )
        deallocate( work )
        call debug('work ok')
        allocate( work(lwork) )
        call debug('Computing the SVD ...')
        call dgesdd('S', nobs, nmem, T, nobs, sig, U, nobs, VT&
            , minMR, work, lwork, iwork, info )

        if(info<0)then
            call stop_program(-info, 'call to dgesdd with illegal argument number ')
        end if
        if(info>0)then
            call stop_program('DBDSDC did not converge, updating process failed')
        end if

        call debug("Done, computing the SVD.")

        !convert singular values to eigenvalues
        sig = sig*sig
        ! Compute number of significant singular values
        sigsum = sum(sig)
        call debug(sigsum, "sum(sig) = ")
        if(isNaN(sigsum)) call debug(sig, 'sig = ')
        sigsum1 = 0.0_dp
        nrsigma = 0
        i = 1
        do while((sigsum1/sigsum < 0.9999_dp).and.(i<=size(sig,1)))
            nrsigma = nrsigma + 1
            sigsum1 = sigsum1 + sig(i)
            i = i+1
        end do

        sig_inv(1:nrsigma) = 1.0_dp/sig(1:nrsigma)
        if(nrsigma<size(sig,1))then
            sig(nrsigma+1:size(sig,1))  = 0.0_dp
            sig_inv(nrsigma+1:size(sig,1))  = 0.0_dp
        end if

        call debug(sig(1:nrsigma), 'significant singular values = ')
        call debug(sig_inv(1:nrsigma), 'significant sig_inv = ')
        call debug(nrsigma, 'nrsigma (number of keeped singular values) = ')

        !computing T1 = S^-1U^t
        allocate(T1(minMR, nobs))
        do i=1, minMR
            T1(i,:) = sig_inv(i)*U(:, i)
        end do

        !computing T2 = T1*D
        call debug('computing T2 = T1*D')
        allocate( T2( minMR, nmem ) )
        call dgemm('N', 'N', minMR, nmem, nobs, 1.0_dp, T1, minMR&
            , D, nobs, 0.0_dp, T2, minMR)

        !computing T3 = U*T2
        call debug('computing T3 = U*T2')
        allocate( T3( nobs, nmem ) )
        call dgemm('N', 'N', nobs, nmem, minMR, 1.0_dp, U, nobs&
            , T2, minMR, 0.0_dp, T3, nobs)

        !deallocate SVD variables
        deallocate( work, iwork, U, VT, sig, sig_inv )

        deallocate(T1, T2)

        !Computing T4 = (HX_dev)^T*T3
        call debug('computing T4 = (HT_dev)^T*T3')
        call dgemm('T', 'N', nmem, nmem, nobs, 1.0_dp, HX_dev, nobs&
            , T3, nobs, 0.0_dp, T4, nmem)
        deallocate(T3)
    end subroutine svd_EnKF_X4

    !> @brief Compute the EnKF analysis using eq. 66 [Evenson 2002]
    !! @param [in,out] T4 updating matrix, defined as X4 in eq. 65 [Evenson 2002]
    !! @param [in,out] X ensemble matrix
    !! @param [in,out] X_dev ensemble deviation from the ensemble mean
    !! @param [in,out] E ensemble model error
    !! @param [in,out] E_dev ensemble deviation from model error mean
    !! @param [in,out] incr (optional) increment flag, says if increment is needed
    !! if increment is needed, it overrides @a X_dev
    !!
    !<
    subroutine EnKF_analysis_using_X4(T4, X, X_dev, E, E_dev, incr)
        real(dp), dimension(:,:), intent(in) :: T4
        real(dp), dimension(:,:), intent(inout) :: X
        real(dp), dimension(:,:), intent(inout) :: X_dev
        real(dp), dimension(:,:), optional, intent(inout) :: E
        real(dp), dimension(:,:), optional, intent(inout) :: E_dev
        logical, optional, intent(in) :: incr
        !local variables
        real(dp), dimension(:,:), allocatable :: X_inc, E_inc
        logical :: ll_incr
        integer :: nmem, ndim, edim

!         if(.not.(
!             ( size(T4,1)==size(T4,2) ) & !T4 must be a square matrix
!             .and.( size(T4,1)==size(X ,2) ) & ! The second dimension of X and E
!             .and.( size(T4,1)==size(E ,2) ) & ! must be equal to the dim of T4
!             .and. all( shape(X)==shape(X_dev) ) & !X_dev has the same shape as X
!             .and. all( shape(E)==shape(E_dev) ) & !E_dev has the same shape as E
!             ) ) then
!             call debug("In EnKF_analysis_using_X4: Incompatible input shapes", tag=dALLWAYS)
!             call debug("Input matrices must the exact same number of columns", tag=dALLWAYS)
!             call debug("  and each matrix and the associated deviation matrix", tag=dALLWAYS)
!             call debug("  must be the exact same size", tag=dALLWAYS)
!             call debug("Got", tag=dALLWAYS)
!             call debug(shape(T4)   , "shape(T4)    = ", tag=dALLWAYS)
!             call debug(shape(X)    , "shape(X)     = ", tag=dALLWAYS)
!             call debug(shape(X_dev), "shape(X_dev) = ", tag=dALLWAYS)
!             call debug(shape(E)    , "shape(E)     = ", tag=dALLWAYS)
!             call debug(shape(E_dev), "shape(E_dev) = ", tag=dALLWAYS)
!         end if

        nmem = SIZE(T4,1)
        ndim = SIZE(X,1)
        edim = SIZE(E,1)

        if( present(incr) )then
            ll_incr = incr
        else
            ll_incr = .false.
        end if

        !Updating X, X = X + X_dev*T4
        call debug('Computing X as X = X + X_dev*T4')
        call debug(shape(T4), "shape(T4)=")
        !computing the increment if requested
        if(ll_incr)then
            call debug('Computing the ensemble increment X_inc + X_dev*T4')
            allocate( X_inc(size(X,1),size(X,2) ) )
            call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X_dev, ndim&
                        , T4, nmem, 0.0_dp, X_inc, ndim)
            X_dev = X_inc
            call debug('Computing X as X = X + X_inc')
            X = X + X_inc
            deallocate(X_inc)
        else
            call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X_dev, ndim&
                , T4, nmem, 1.0_dp, X, ndim)
        end if

        !Updating the model error in the case of state augmentation
        if( present(E) )then
            call debug('updating E, E = E + E_dev*T4')
            call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E_dev, edim&
                , T4, nmem, 1.0_dp, E, edim)
            !computing the increment if requested
            if(ll_incr)then
                allocate( E_inc(size(E,1),size(E,2) ) )
                call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E_dev, edim&
                            , T4, nmem, 0.0_dp, E_inc, edim)
                E_dev = E_inc
                E = E + E_inc
                deallocate(E_inc)
            else
                call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E_dev, edim&
                    , T4, nmem, 1.0_dp, E, edim)
            end if
        end if
    end subroutine EnKF_analysis_using_X4

    !> @brief Compute the EnKF analysis using eq. 72 [Evenson 2002]
    !! @param [in,out] T4 updating matrix, defined as X4 in eq. 65 [Evenson 2002]
    !! @param [in,out] X ensemble matrix
    !! @param [in,out] E ensemble model error
    !! @param [in,out] X_inc (optional) ensemble increment
    !! @param [in,out] E_inc (optional) ensemble model error increment
    !! @details Eq. 72 [Evenson 2002] defines a new version of the analysis
    !! equation that uses only the Ensemble states. If the equation is
    !! developped, one sees that it is the same as eq. 66 with the ensemble
    !! deviation replaced by the ensemble state. See [Evenson 2002] for
    !! the simplifications made
    !!
    !! @todo computes X5 explicitly and avoid using two ensembles.
    !!   With the right implementation of this version, it is not
    !!   possible to compute the increments
    !!
    !<
    subroutine EnKF_analysis_using_X5(T4, X, E, X_inc, E_inc)
        real(dp), dimension(:,:), intent(in) :: T4
        real(dp), dimension(:,:), intent(inout) :: X
        real(dp), dimension(:,:), optional, intent(inout) :: X_inc
        real(dp), dimension(:,:), optional, intent(inout) :: E
        real(dp), dimension(:,:), optional, intent(inout) :: E_inc
        !local variables
        real(dp), dimension(:,:), allocatable :: X_tmp, E_tmp
        integer :: nmem, ndim, edim

!         if(.not.(
!             ( size(T4,1)==size(T4,2) ) & !T4 must be a square matrix
!             .and.( size(T4,1)==size(X ,2) ) & ! The second dimension of X and E
!             .and.( size(T4,1)==size(E ,2) ) & ! must be equal to the dim of T4
!             .and. all( shape(X)==shape(X_dev) ) & !X_dev has the same shape as X
!             .and. all( shape(E)==shape(E_dev) ) & !E_dev has the same shape as E
!             ) ) then
!             call debug("In EnKF_analysis_using_X4: Incompatible input shapes", tag=dALLWAYS)
!             call debug("Input matrices must the exact same number of columns", tag=dALLWAYS)
!             call debug("  and each matrix and the associated deviation matrix", tag=dALLWAYS)
!             call debug("  must be the exact same size", tag=dALLWAYS)
!             call debug("Got", tag=dALLWAYS)
!             call debug(shape(T4)   , "shape(T4)    = ", tag=dALLWAYS)
!             call debug(shape(X)    , "shape(X)     = ", tag=dALLWAYS)
!             call debug(shape(X_dev), "shape(X_dev) = ", tag=dALLWAYS)
!             call debug(shape(E)    , "shape(E)     = ", tag=dALLWAYS)
!             call debug(shape(E_dev), "shape(E_dev) = ", tag=dALLWAYS)
!         end if

        nmem = SIZE(T4,1)
        ndim = SIZE(X,1)
        edim = SIZE(E,1)

        !Updating X, X = X*T5.
        !T5 = I + T4.
        !X = X + X*T4
        call debug('Computing X as X = X*T5')
        allocate( X_tmp(size(X,1),size(X,2) ) )
        call debug('Computing the ensemble increment X_inc = X*T4')
        call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X, ndim&
            , T4, nmem, 0.0_dp, X_tmp, ndim)
        call debug('Computing X as X = X + X_tmp')
        X = X + X_tmp
        !computing the increment if requested
        if( present(X_inc) )then
            X_inc = X_tmp
        end if
        deallocate(X_tmp)

        !Updating the model error in the case of state augmentation
        if( present(E) )then
            call debug('updating E, E = E + E_dev*T4')
            allocate( E_tmp(size(E,1),size(E,2) ) )
            call debug('Computing the model error increment E_inc = E*T4')
            call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E, edim&
                , T4, nmem, 0.0_dp, E_tmp, edim)
            E = E + E_tmp
            !computing the increment if requested
            if( present(E_inc) )then
                E_inc = E_tmp
            end if
            deallocate(E_tmp)
        end if
    end subroutine EnKF_analysis_using_X5

    !> @brief Computes the analysed ensemble in the EnKF, using SVD
    !! @param [in,out] X ensemble matrix
    !! @param [in] X_dev ensemble deviation from the ensemble mean
    !! @param [in] HX_dev observation operator applied to Xdev
    !! @param [in] Y_dev ensemble of observation perturbation (each column)
    !! @param [in] D ensemble innovations (each column is an innovation)
    !! @param [in,out] E ensemble model error
    !! @param [in] E_dev ensemble deviation from model error mean
    !! @param [in] incr boolean flag saying if increment is computed or not
    !!   The increment is the update to be applied to the background to get
    !!   the analysis. By default the increment is not computed. If computed,
    !!   the increment is returned in X_dex, resp. E_dev.
    !!
    !! @details
    !!   the number of rows of X defines the size of the state variable
    !!   the number of columns of X defines the size of the ensemble
    !!   \a E and \a E_dev are used in conjunction with state augmentation.
    !!   They have the same shape as X. If ensemble augmentation with model error
    !!   is required, both should be present, otherwise, both should be absent
    !!
    !!  This is the version before September 09th 2014
    !<
    subroutine svd_EnKFAnalysis_original(X, X_dev, HX_dev, Y_dev, D, E, E_dev, incr)
        real(dp), dimension(:,:), intent(inout) :: X
        real(dp), dimension(:,:), intent(inout) :: X_dev
        real(dp), dimension(:,:), intent(in) :: HX_dev, Y_dev, D
        real(dp), dimension(:,:), optional, intent(inout) :: E
        real(dp), dimension(:,:), optional, intent(inout) :: E_dev
        logical, optional, intent(in) :: incr
        !local variables
        real(dp), dimension(:,:), allocatable :: T, U, VT, T1, T2, T3, T4, T5
        real(dp), dimension(:), allocatable :: sig, sig_inv
        !debugging variables
        !real(dp), dimension(:,:), allocatable :: Tb, U2
        !real(dp), dimension(:), allocatable :: bt, xa, xt
        !integer :: j
        !end of debugging variables
        !SVD variables
        integer, dimension(:), allocatable :: IWORK
        real(dp), dimension(:), allocatable :: WORK
        real(dp) :: sigsum, sigsum1
        integer :: nrsigma
        !dimension variables
        integer :: ndim !! Dimension of model state
        integer :: nmem !! Number of ensemble members
        integer :: edim !! Dimension of model error
        integer :: nobs !! Number of observations
        integer :: minMR !!min(nmem, nobs)
        integer :: lwork, info, i
        !inflation factor: I do not remember exactely what it is, I reset it from 2 to 1
        real(dp), parameter :: iFact = 1.0_dp
        !dgemm factor
        real(dp):: alpha!, beta
        logical :: ll_incr


        !initialization
        ndim = SIZE(X,1)
        nmem = SIZE(X,2)
        nobs = SIZE(D,1)
        minMR = min(nmem, nobs)
        if( present(E) )then
            edim = size(E,1)
        else
            edim = -1
        end if
        if( present(incr) )then
            ll_incr = incr
        else
            ll_incr = .FALSE.
        end if

        !debugging
        call debug(nobs, 'nobs = ')
        call debug(sum(abs(D)), 'sum(abs(D)) = ')
        call debug(sum(abs(X)), 'sum(abs(X)) = ')
        call debug(sum(abs(X_dev)), 'sum(abs(X_dev)) = ')
        call debug(sum(abs(HX_dev)), 'sum(abs(HX_dev)) = ')
        call debug(sum(abs(Y_dev)), 'sum(abs(Y_dev)) = ')
        !end of debugging

        !computing T = HX_dev + Y_dev
        allocate( T(nobs, nmem) )
        T = sqrt(iFact)*HX_dev + Y_dev

    !     !debugging
    !     !test variables
    !       allocate( xt(nobs), xa(nobs), bt(nobs), Tb(nobs, nobs) )
    !       xt = 0.0_dp
    !       xa = 0.0_dp
    !       bt = 0.0_dp
    !       call random_number(xt)
    !       Tb = matmul( T,transpose(T) )
    !       bt = matmul( Tb, xt )
    !     !end of test variables
    !     !end of debugging

        !computing the SVD of T
        allocate( U(nobs,minMR), VT(minMR, nmem), sig(minMR)&
            , sig_inv(minMR), iwork(8*minMR), work(1) )
        U = 0.0_dp
        VT = 0.0_dp
        sig = 0.0_dp
        lwork = -1 !! query the optimal size of WORK
        call dgesdd('S', nobs, nmem, T, nobs, sig, U, nobs, VT&
            , minMR, work, lwork, iwork, info )
        lwork = int( work(1) )
        deallocate( work )
        call debug('work ok')
        allocate( work(lwork) )
        call debug('Computing the SVD ...')
        call dgesdd('S', nobs, nmem, T, nobs, sig, U, nobs, VT&
            , minMR, work, lwork, iwork, info )

        if(info<0)then
            call stop_program(-info, 'call to dgesdd with illegal argument number ')
        end if
        if(info>0)then
            call stop_program('DBDSDC did not converge, updating process failed')
        end if

        call debug("Done, computing the SVD.")

        !convert singular values to eigenvalues
        sig = sig*sig
        !call debug(sig, "sig = ")
        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! Compute number of significant singular values
        sigsum = sum(sig)
        call debug(sigsum, "sum(sig) = ")
        if(isNaN(sigsum)) call debug(sig, 'sig = ')
        sigsum1 = 0.0_dp
        nrsigma = 0
        i = 1
        do while((sigsum1/sigsum < 0.9999_dp).and.(i<=size(sig,1)))
            nrsigma = nrsigma + 1
            sigsum1 = sigsum1 + sig(i)
            i = i+1
        end do

        sig_inv(1:nrsigma) = 1.0_dp/sig(1:nrsigma)
        if(nrsigma<size(sig,1))then
            sig(nrsigma+1:size(sig,1))  = 0.0_dp
            sig_inv(nrsigma+1:size(sig,1))  = 0.0_dp
        end if

        call debug(sig(1:nrsigma), 'significant singular values = ')
        call debug(sig_inv(1:nrsigma), 'significant sig_inv = ')
        call debug(nrsigma, 'nrsigma (number of keeped singular values) = ')

        !computing T1 = S^-1U^t
        allocate(T1(minMR, nobs))
        do i=1, minMR
            T1(i,:) = sig_inv(i)*U(:, i)
        end do

    !     !debugging
    !       !computing U2 = US
    !       allocate( U2(nobs,minMR) )
    !       do j=1,size(sig,1)
    !         U2(:,j) = U(:,j)*sig(j)
    !       end do
    !       !call debug(sig, "sig = ")
    !       call debug( norm2( Tb - matmul(U2,transpose(U)) ), 'residual norm (Tb-USU^t)')
    !       deallocate( U2 )
    !       call debug("========================================")
    !
    !       !test variables
    !       xa = matmul(U, matmul(T1, bt) )
    !       call debug("========================================")
    !       call debug(norm2(xt-xa)/norm2(xa), 'norm2(xt-xa)/norm2(xa) = ')
    !       call debug("========================================")
    !       allocate( xt, xa, bt, Tb )
    !       !end of test variables
    !    !end of debugging

        !computing T2 = T1*D
        call debug('computing T2 = T1*D')
        allocate( T2( minMR, nmem ) )
        call dgemm('N', 'N', minMR, nmem, nobs, 1.0_dp, T1, minMR&
            , D, nobs, 0.0_dp, T2, minMR)

        !computing T3 = U*T2
        call debug('computing T3 = U*T2')
        allocate( T3( nobs, nmem ) )
        call dgemm('N', 'N', nobs, nmem, minMR, 1.0_dp, U, nobs&
            , T2, minMR, 0.0_dp, T3, nobs)

        !deallocate SVD variables
        deallocate( work, iwork, U, VT, sig, sig_inv )

        deallocate(T1, T2)

        !Computing T4 = (HX_dev)^T*T3
        call debug('computing T4 = (HT_dev)^T*T3')
        allocate( T4(nmem, nmem) )
        call dgemm('T', 'N', nmem, nmem, nobs, 1.0_dp, HX_dev, nobs&
            , T3, nobs, 0.0_dp, T4, nmem)
        deallocate(T3)

        !Updating X, X = X + X_dev*T4
        alpha = iFact
        call debug('updating X, X = X + X_dev*T4')
        call dgemm('N', 'N', ndim, nmem, nmem, alpha, X_dev, ndim&
            , T4, nmem, 1.0_dp, X, ndim)
        !computing the increment if requested
        if(ll_incr)then
        allocate( T5(size(X,1),size(X,2) ) )
        call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X_dev, ndim&
                    , T4, nmem, 0.0_dp, T5, ndim)
        X_dev = T5
        deallocate(T5)
        end if

        !Updating the model error in the case of state augmentation
        if( present(E) )then
        call debug('updating E, E = E + E_dev*T4')
        call dgemm('N', 'N', edim, nmem, nmem, alpha, E_dev, edim&
            , T4, nmem, 1.0_dp, E, edim)
        !computing the increment if requested
        if(ll_incr)then
            allocate( T5(size(E,1),size(E,2) ) )
            call dgemm('N', 'N', edim, nmem, nmem, 1.0_dp, E_dev, edim&
                        , T4, nmem, 0.0_dp, T5, edim)
            E_dev = T5
            deallocate(T5)
        end if
        end if

        deallocate(T4)
    end subroutine svd_EnKFAnalysis_original

    !> \brief Pure SVD version, must be pure eigenvalues version...
    !! \param[in,out] X ensemble matrix
    !! \param[in] X_dev ensemble deviation from the ensemble mean
    !! \param[in] HX_dev observation operator applied to Xdev
    !! \param[in] Y_dev ensemble of observation perturbation (each column)
    !! \param[in] D ensemble innovations (each column is an innovation)
    !! \details
    !!   the number of rows of X defines the size of the state variable
    !!   the number of columns of X defines the size of the ensemble
    !!   Computes the analysed ensemble in the EnKF, using SVD
    !<
    subroutine svd_EnKFAnalysis_first(X, X_dev, HX_dev, Y_dev, D)
        real(dp), dimension(:,:), intent(inout) :: X
        real(dp), dimension(:,:), intent(in) :: X_dev, HX_dev, Y_dev, D
        !local variables
        real(dp), dimension(:,:), allocatable :: T, U, U2, VT, T1, T2, T3, Tb
        real(dp), dimension(:), allocatable :: sig, sig_inv, xt, bt, xa
        !SVD variables
        integer, dimension(:), allocatable :: IWORK
        real(dp), dimension(:), allocatable :: WORK
        real(dp) :: sigsum, sigsum1, alpha
        integer :: nrsigma
        !dimension variables
        integer :: ndim !! Dimension of model state
        integer :: nmem !! Number of ensemble members
        integer :: nobs !! Number of observations
        integer :: lwork, info, i

        !initialization
        ndim = SIZE(X,1)
        nmem = SIZE(X,2)
        nobs = SIZE(D,1)

        !debugging
        call debug(nobs, 'nobs = ')
        call debug(sum(abs(D)), 'sum(abs(D)) = ')
        call debug(sum(abs(X)), 'sum(abs(X)) = ')
        call debug(sum(abs(X_dev)), 'sum(abs(X_dev)) = ')
        call debug(sum(abs(HX_dev)), 'sum(abs(HX_dev)) = ')
        call debug(sum(abs(Y_dev)), 'sum(abs(Y_dev)) = ')
        !end of debugging
        !computing T = (HX_dev)*(HX_dev)^T
        allocate( T(nobs,nobs) )
        alpha = 1.0_dp/real(nmem-1,dp)
        T = 0.0_dp
        call dgemm('N', 'T', nobs, nobs, nmem, 1.0_dp, HX_dev, nobs&
                    , HX_dev, nobs, 0.0_dp, T, nobs)

        !computing T = T + Y_dev * (Y_dev)^T
        call dgemm('N', 'T', nobs, nobs, nmem, 1.0_dp, Y_dev, nobs&
                    , Y_dev, nobs, 1.0_dp, T, nobs)

        !test variables
        allocate( xt(nobs), xa(nobs), bt(nobs), Tb(nobs,nobs) )
        xt = 0.0_dp
        xa = 0.0_dp
        bt = 0.0_dp
        call random_number(xt)
        bt = matmul(T, xt)
        Tb = T
        !end of test variables

        !computing the SVD of T
        allocate( U(nobs,nobs), VT(nobs,nobs), sig(nobs)&
            , sig_inv(nobs), iwork(8*min(nobs,nobs)), work(1) )
        U = 0.0_dp
        VT = 0.0_dp
        sig = 0.0_dp
        lwork = -1 !! query the optimal size of WORK
        call dgesdd('A', nobs, nobs, T, nobs, sig, U, nobs, VT, nobs&
            , work, lwork, iwork, info )
        lwork = int( work(1) )
        deallocate( work )
        allocate( work(lwork) )
        call dgesdd('A', nobs, nobs, T, nobs, sig, U, nobs, VT, nobs&
            , work, lwork, iwork, info )

        if(info<0)then
            call stop_program(-info, 'call to dgesdd with illegal argument number ')
        end if
        if(info>0)then
            call stop_program('DBDSDC did not converge, updating process failed')
        end if

        call debug("after SVD")

        !call debug(sig, "sig = ")
        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
        ! Compute number of significant singular values
        sigsum = sum(sig)
        call debug(sigsum, "sum(sig) = ")
        if(isNaN(sigsum)) call debug(sig, 'sig = ')
        sigsum1 = 0.0_dp
        nrsigma = 0
        i = 1
        do while((sigsum1/sigsum < 0.999_dp).and.(i<=size(sig,1)))
            nrsigma = nrsigma + 1
            sigsum1 = sigsum1 + sig(i)
            i = i+1
        end do

        sig_inv(1:nrsigma) = 1.0_dp/sig(1:nrsigma)
        if(nrsigma<size(sig,1))then
            sig(nrsigma+1:size(sig,1))  = 0.0_dp
            sig_inv(nrsigma+1:size(sig,1))  = 0.0_dp
        end if

        call debug(sig(1:nrsigma), 'significant singular values = ')
        call debug(sig_inv(1:nrsigma), 'significant sig_inv = ')
        call debug(nrsigma, 'nrsigma (number of keeped singular values) = ')

        !computing U2 = US
        allocate( U2(nobs,nobs) )
        do i=1,size(sig,1)
            U2(:,i) = U(:,i)*sig(i)
        end do
        !call debug(sig, "sig = ")
        call debug( norm2( Tb - matmul(U2,VT) ), 'residual norm (Tb-USV^t)' )
        deallocate( U2 )
        call debug( "========================================" )
        !computing U = US^-1
        do i=1,size(sig,1)
            U(:,i) = U(:,i)*sig_inv(i)
        end do

        !test variables
        xa = matmul( transpose(VT), matmul(transpose(U), bt) )
        call debug("========================================")
        call debug(norm2(xt-xa), 'norm2(xt-xa) = ')
        call debug("========================================")
        !end of test variables

        !computing T1 = UT*D
        call debug('computing T1 = UT*D')
        allocate( T1( nobs, nmem ) )
        call dgemm('T', 'N', nobs, nmem, nobs, 1.0_dp, U, nobs&
            , D, nobs, 0.0_dp, T1, nobs)

        !computing T2 = V*T1
        call debug('computing T2 = V*T1')
        allocate( T2( nobs, nmem ) )
        call dgemm('T', 'N', nobs, nmem, nobs, 1.0_dp, VT, nobs&
            , T1, nobs, 0.0_dp, T2, nobs)

        !deallocate SVD variables
        deallocate( work, iwork, U, VT, sig, sig_inv )

        deallocate(T1)

        !computing T3 = (HX_dev)^T*T2
        call debug('computing T3 = (HT_dev)^T*T2')
        allocate( T3(nmem, nmem) )
        call dgemm('T', 'N', nmem, nmem, nobs, 1.0_dp, HX_dev, nobs&
            , T2, nobs, 0.0_dp, T3, nmem)
        deallocate(T2)

        !updating X, X = X + X_dev*T3
        call debug('updating X, X = X + X_dev*T3')
        call dgemm('N', 'N', ndim, nmem, nmem, 1.0_dp, X_dev, ndim&
            , T3, nmem, 1.0_dp, X, ndim)
        deallocate(T3)
    end subroutine svd_EnKFAnalysis_first

    ! +==========================================================
    ! RBCG version
    !! Ensemble Kalman Filter analysis using conjugate gradient algorithm
    !! This module implements the Ensemble Kalman Filter as decribed
    !! in Geir Evenson, The Ensemble Kalman Filter: theoretical formulation
    !! and practical implementation. Ocean Dynamics (2003) 53: 343-367
    !! DOI 10.1007/s10236-003-9
    !! The paper is refered to as Evenson2003 in the module.
    !! The main difference is that the RBCG algorithm
    !! is used to solved the system in the place of the SVD
    !<

    !> @brief Compute the matrix X4 eq. 65 [Evenson 2002] using rbcg
    !! @param [in] HX_dev observation operator applied to Xdev
    !! @param [in] Y_dev ensemble of observation perturbation (each column)
    !! @param [in] D ensemble innovations (each column is an innovation)
    !! @param [in, out] T4 X4 matrix (named T4 for local purpose)
    !! @param [in] obs_std standard deviation of observations
    !! @param [in] nItrMax optional maximum number of iteration
    !!    of the conjugate gradient, default 100
    !! @param [in] tol optional tolerance minimum of the CG
    !!    ratio of gradient, default 10^-16
    !!
    !! @details
    !!   the number of columns of HX_dev defines the size of the ensemble
    !!   T4 must be a square matrix with the number of columns/rows given
    !!   by the number of the ensemble members
    !!   @a state_coord, @a obs_coord and @a cls are used for localization
    !!   in the covariance. There is no special requirement on the units
    !!   of the coordinates, but they must be consistent with the @a cls.
    !!   for example, the coordinates can be the indices in the model grid
    !!   and the @a cls be given in terms of number of grid points.
    !!   the coordinates can also be the geographical coordinates
    !!   (lon/lat/alt etc.).
    !!   This subroutine does not know anything about the number of physical
    !!   dimensions of the system, this allows the user to apply the
    !!   localization only in the chosen direction(s) for exemple, one can
    !!   choose to localize only along the x direction, or along x and y.
    !!   The number of directions for the localization is defined by the
    !!   number of column (second dimension) of the coordinate vectors.
    !!   The firs dimension correspond to the number of element of state
    !!   or observation vector.
    !!   @a cls is a 1-D vector of size the number of directions of
    !!   the localization
    !!
    !! @note for now there is no check on the shape of matrix parameters
    !<
    subroutine rbcg_EnKF_X4(HX_dev, Y_dev, D, T4, obs_std, nItrMax, tol)
        real(dp), dimension(:,:), intent(in) :: HX_dev, Y_dev, D
        real(dp), dimension(:,:), intent(in out) :: T4
        real(dp), dimension(:), intent(in) :: obs_std
        integer, optional, intent(in) :: nItrMax
        real(dp), optional, intent(in) :: tol
        !local variables
        real(dp), dimension(size(HX_dev,1)) :: x
        !dimension variables
        integer :: nmem !! Number of ensemble members
        integer :: edim !! Dimension of model error
        integer :: nobs !! Number of observations
        integer :: i,j
        !
        integer :: il_nItrMax
        real(dp):: rl_tol

        !initialization
        nmem = SIZE(HX_dev,2)
        nobs = SIZE(D,1)
        if( present(nItrMax) )then
            il_nItrMax = nItrMax
        else
            il_nItrMax = 100
        end if
        if( present(tol) )then
            rl_tol = tol
        else
            rl_tol = 1.0d-16
        end if

        !debugging
        call debug(nobs, 'nobs = ')
        call debug(sum(abs(HX_dev)), 'svd_EnKF_X4; sum(abs(HX_dev))')
        call debug(sum(abs(Y_dev)) , 'svd_EnKF_X4; sum(abs(Y_dev)).')
        call debug(sum(abs(D))     , 'svd_EnKF_X4; sum(abs(D)).....')
        !end of debugging
        T4 = 999.0
        do j = 1, nmem
            !Solve the linear system for one innovation member
            call debug(j, "RBCG, Solving the linear system for member", tag=dALLWAYS)
            call rbcg_solve(nobs, nmem, HX_dev, obs_std, D(:,j), x, il_nItrMax, rl_tol)

            ! multiply the result by the transpose of HX_dev
            ! this does not apply the localization
            do i = 1, nmem
                T4(i,j) = dot_product( HX_dev(:,i), x )
            end do
        end do
    end subroutine rbcg_EnKF_X4

    !> @brief Solve the Ensemble Kalman Filter system for one member using RBCG
    !! @param [in] nRow number of row of HA, this is the size of the system
    !! @param [in] nCol number of columns of HA
    !! @param [in] HA observation operator applied to ensemble deviation from the mean
    !! @param [in] obs_std vector of standard deviation of observation error
    !! @param [in] d innovation of one member of the ensemble
    !! @param [in,out] x solution of the system
    !! @param [in] nItrMax maximum number of iterations
    !! @param [in] tol stopping criteria based on the B-norm of the gradient
    !!
    !! @details This algorithm solves only for one member, the user must call it
    !! in a loop to solve for all members
    !! It is important to recall that is only the linear systems in the
    !! EnKF analysis equation: (HAA^TH^T + R)^-1D, see eq 54 Evenson2003
    !<
    subroutine rbcg_solve(nRow, nCol, HA, obs_std, d, x, nItrMax, tol)
        integer, intent(in) :: nRow, nCol
        real(dp), dimension(nRow, nCol), intent(in) :: HA
        real(dp), dimension(nRow), intent(in) :: obs_std
        real(dp), dimension(nRow), intent(in) :: d
        real(dp), dimension(nRow), intent(in out) :: x
        integer, optional, intent(in) :: nItrMax
        real(dp), optional, intent(in) :: tol
        !local variables
        !RBCG variables
        real(dp) :: cg_r(nRow)
        real(dp) :: cg_rtmp(nRow)
        real(dp) :: cg_lambda(nRow)
        real(dp) :: cg_t(nRow)
        real(dp) :: cg_beta
        real(dp) :: cg_p(nRow)
        real(dp) :: cg_alpha
        real(dp) :: cg_w(nRow)
        real(dp) :: cg_q(nRow)
        !others
        real(dp) :: r0_norm ! B-norm of the residue at the zeroth iteration
        real(dp) :: ri_norm ! B-norm of the residue at the current iteration
        real(dp) :: cg_wTr !dot product w by r (i)
        real(dp) :: cg_qTt !dot product q by t
        real(dp) :: cg_wTrp1 !dot product w by r (i+1)
        real(dp) :: cg_gNorm !b-norm of the gradient
        real(dp) :: cg_omega, cg_epsilon
        integer :: cg_itr
        integer :: j
        logical :: converged
        character(*), parameter :: fmt0= "('*** RBCG   Iter ***',I3,' -- |r_0|=', ,ES14.5E3)"
        character(*), parameter :: fmti= "('*** RBCG   Iter ***',I3,' -- |r|/|r_0|=', ,ES14.5E3)"

        call debug(tol, "RBCG - tol=", tag=dALLWAYS)
        call debug(nItrMax, "RBCG - nItrMax=", tag=dALLWAYS)
        cg_epsilon = epsilon(1.0_dp)
        cg_itr = 0
        cg_lambda = 0.0
        cg_r = d/obs_std**2
        cg_p = cg_r
        !computing HA(HA)^Tr
        call matrixmult(nRow, nCol, HA, cg_r, cg_w)
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
            !The following two quantities are used as denominator
            !So if one is zero, either, there is convergence or something is wrong
            cg_wTr = dot_product(cg_w, cg_r)
            cg_qTt = dot_product(cg_q, cg_t)
            if( (cg_qTt>cg_epsilon).and.(cg_wTr>cg_epsilon) )then
                cg_alpha = cg_wTr/cg_qTt
                cg_lambda = cg_lambda + cg_alpha*cg_p
                cg_r = cg_r - cg_alpha*cg_q
                !re-conjugate r
                call matrixmult(nRow, nCol, HA, cg_r, cg_w)
                cg_wTrp1 = dot_product(cg_w, cg_r)
                cg_beta = cg_wTrp1/cg_wTr
                if( cg_beta<cg_epsilon )then
                    converged = .true.
                end if
                cg_p = cg_r + cg_beta*cg_p
                cg_t = cg_w + cg_beta*cg_t
                ri_norm = sqrt(cg_wTrp1)
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
            call debug(cg_itr, "RBCG converged, total nb iterations", tag=dALLWAYS)
            call debug(sum(abs(cg_lambda)), "Solution sum(abs(cg_lambda))")
            if(cg_qTt<=cg_epsilon)call debug(cg_qTt, "q*t is too small", tag=dALLWAYS)
            if(cg_wTr<=cg_epsilon)call debug(cg_wTr, "w*r is too small", tag=dALLWAYS)
            if(cg_beta<=cg_epsilon)call debug(cg_beta, "beta is too small", tag=dALLWAYS)
        end if

        x = cg_lambda
    end subroutine rbcg_solve

    !> @brief apply the the matrix AA^T to x and returns the result in y
    !! @param [in] nRow number of row of A, this is the size of the system
    !! @param [in] nCol number of columns of A
    !! @param [in] A factor of the matrix of the system
    !! @param [in] x input vector
    !! @param [in,out] y output vector
    !<
    subroutine matrixmult(nRow, nCol, A, x, y)
        integer, intent(in) :: nRow, nCol
        real(dp), dimension(nRow, nCol), intent(in) :: A
        real(dp), dimension(nRow), intent(in) :: x
        real(dp), dimension(nRow), intent(in out) :: y
        !local variables
        integer :: j

        y = 0.0
        do j = 1, nCol
            y = y + dot_product( A(:,j), x)*A(:,j)
        end do
    end subroutine matrixmult

end module svd_EnKF
