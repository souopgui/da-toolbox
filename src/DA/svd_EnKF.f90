!> \file svd_EnKF.f90
!! Ensemble Kalman FIlter analysis using SVD decomposition
!<

!>\brief ensemble Kalman filter analysis
module svd_EnKF
  use debug_tools
  use general_tools
implicit none

contains

  !> \brief Computes the analysed ensemble in the EnKF, using SVD
  !! \param[in,out] X ensemble matrix
  !! \param[in] X_dev ensemble deviation from the ensemble mean
  !! \param[in] HX_dev observation operator applied to Xdev
  !! \param[in] Y_dev ensemble of observation perturbation (each column)
  !! \param[in] D ensemble innovations (each column is an innovation)
  !! @param [in,out] E ensemble model error
  !! @param [in] E_dev ensemble deviation from model error mean
  !!
  !! \details
  !!   the number of rows of X defines the size of the state variable
  !!   the number of columns of X defines the size of the ensemble
  !!   \a E and \a E_dev are used in conjunction with state augmentation. They have the same shape as X. If ensemble augmentation with model error is required, both should be present, otherwise, both should be absent
  !<
  subroutine svd_EnKFAnalysis(X, X_dev, HX_dev, Y_dev, D, E, E_dev)
    real(dp), dimension(:,:), intent(inout) :: X
    real(dp), dimension(:,:), intent(in) :: X_dev, HX_dev, Y_dev, D
    real(dp), dimension(:,:), optional, intent(inout) :: E
    real(dp), dimension(:,:), optional, intent(in) :: E_dev
    !local variables
    real(dp), dimension(:,:), allocatable :: T, U, VT, T1, T2, T3, T4
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
    integer :: nobs !! Number of observations
    integer :: minMR !!min(nmem, nobs)
    integer :: lwork, info, i
    !inflation factor
    real(dp), parameter :: iFact = 2.0_dp
    !dgemm factor
    real(dp):: alpha!, beta


    !initialization
    ndim = SIZE(X,1)
    nmem = SIZE(X,2)
    nobs = SIZE(D,1)
    minMR = min(nmem, nobs)

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
    do while((sigsum1/sigsum < 0.99_dp).and.(i<=size(sig,1)))
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

    !computing T4 = (HX_dev)^T*T3
    call debug('computing T4 = (HT_dev)^T*T3')
    allocate( T4(nmem, nmem) )
    call dgemm('T', 'N', nmem, nmem, nobs, 1.0_dp, HX_dev, nobs&
      , T3, nobs, 0.0_dp, T4, nmem)
    deallocate(T3)

    !updating X, X = X + X_dev*T4
    alpha = iFact
    call debug('updating X, X = X + X_dev*T4')
    call dgemm('N', 'N', ndim, nmem, nmem, alpha, X_dev, ndim&
      , T4, nmem, 1.0_dp, X, ndim)
    !updating model error in case of state augmentation
    if( present(E) )then
      call debug('updating E, E = E + E_dev*T4')
      call dgemm('N', 'N', ndim, nmem, nmem, alpha, E_dev, ndim&
        , T4, nmem, 1.0_dp, E, ndim)
    end if
    deallocate(T4)
  end subroutine svd_EnKFAnalysis

  !> \brief Pure SVD version, must be pure eigenvalues version... Computes the analysed ensemble in the EnKF, using SVD
  !! \param[in,out] X ensemble matrix
  !! \param[in] X_dev ensemble deviation from the ensemble mean
  !! \param[in] HX_dev observation operator applied to Xdev
  !! \param[in] Y_dev ensemble of observation perturbation (each column)
  !! \param[in] D ensemble innovations (each column is an innovation)
  !! \details
  !!   the number of rows of X defines the size of the state variable
  !!   the number of columns of X defines the size of the ensemble
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

end module svd_EnKF
