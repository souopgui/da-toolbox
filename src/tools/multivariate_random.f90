!> \file multivariate_random.f90
!! Defines tools to sample multivariate random variables
!<

!> Module defining procedure to sample multivariate random variables
!!
!<
module multivariate_random
  use general_constant
  use evd_covariance
  use random
implicit none

  !> Interface for multivariate normal variables
  interface sample_normal
    module procedure sample_evd_normal_zero_mean
    module procedure sample_evd_normal_constant_mean
    module procedure sample_evd_normal_variable_mean
  end interface sample_normal
contains

  !> Samples a multivariate random variable \a x of mean 0 and covariance matrix \a C
  !! @param [out] x is a real vector that contains in returns the random (vector) variable
  !! @param [in] C is a \a evd_cov_matrix that gives the covariance matrix
  !!
  !! \details The subroutine assumes that all variable are of the proper shape and size
  !!
  !<
  subroutine sample_evd_normal_zero_mean(x, C)
    real(dp), dimension(:), intent(out) :: x
    type(evd_cov_matrix), intent(in) :: C
    !local variables
    real(dp), dimension(:), allocatable :: rla_rand

    allocate( rla_rand(get_evd_cov_nsv(C)) )
    !samples a vector of independent standard normal random variables
    CALL rand_normal( rla_rand )
    call apply_evd_square_cov(C, rla_rand, x)
    deallocate(rla_rand)
  end subroutine sample_evd_normal_zero_mean

  !> Samples a multivariate random variable \a x of constant scalar mean \a mean and covariance matrix \a C
  !! @param [out] x is a real vector that contains in returns the random (vector) variable
  !! @param [in] mean is a real number that gives the mean of the distribution
  !! @param [in] C is a \a evd_cov_matrix that gives the covariance matrix
  !!
  !! \details The subroutine assumes that all variable are of the proper shape and size
  !!
  !<
  subroutine sample_evd_normal_constant_mean(x, mean, C)
    real(dp), dimension(:), intent(out) :: x
    real(dp), intent(in) :: mean
    type(evd_cov_matrix), intent(in) :: C
    !local variables
    real(dp), dimension(:), allocatable :: rla_rand

    call sample_evd_normal_zero_mean(x, C)
    x = x + mean
    deallocate(rla_rand)
  end subroutine sample_evd_normal_constant_mean

  !> Samples a multivariate random variable \a x of mean \a mean and covariance matrix \a C
  !! @param [out] x is a real vector that contains in returns the random (vector) variable
  !! @param [in] mean is a real number that gives the mean of the distribution
  !! @param [in] C is a \a evd_cov_matrix that gives the covariance matrix
  !!
  !! \details The subroutine assumes that all variable are of the proper shape and size
  !!
  !<
  subroutine sample_evd_normal_variable_mean(x, mean, C)
    real(dp), dimension(:), intent(out) :: x
    real(dp), dimension(:), intent(in)  :: mean
    type(evd_cov_matrix), intent(in) :: C
    !local variables
    real(dp), dimension(:), allocatable :: rla_rand

    call sample_evd_normal_zero_mean(x, C)
    x = x + mean
    deallocate(rla_rand)
  end subroutine sample_evd_normal_variable_mean
end module multivariate_random