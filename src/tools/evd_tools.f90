!> \file evd_tools.f90
!! Eigenvalues decomposition tools
!<

module evd_tools
  use general_constant
  !use general_tools
  use debug_tools

implicit none

contains

  !> Computes the k larges eigenvalues and associated eigenvectors
  !! @param [in, out] A square symmetric matrix whose eigenvalues and eigenvectors are desired
  !! @param [in] k Number of eigenvalues and eigenvectors to be returned
  !! @param [out] U matrix whose columns are the eigenvectors
  !! @param [out] L vector of eigenvalues
  !!
  !! \details This procedure makes a call to lapack dsyev and retains only the \a k largest
  !! The matrix A is destroyed on return. There is no check on parameters size
  !! For more details on the lapack symmetric eigenvalues decomposition. see http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html
  subroutine symmetric_ev(A, k, U, L)
    real(dp), dimension(:,:), intent(in out) :: A
    integer, intent(in) :: k
    real(dp), dimension(:,:), intent(out) :: U
    real(dp), dimension(:), intent(out) :: L
    !local variables
    real(dp), dimension(:), allocatable :: W, work
    integer :: n, info, lwork

    n = size(A,1)
    allocate( W(n), work(1) )
    lwork = -1
    call dsyev ('V', 'U', n, A, n, W, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsyev ('V', 'U', n, A, n, W, work, lwork, info)
    if(info<0)then
      call stop_program(-info, 'call to dsyev with illegal argument number')
    end if
    if(info>0)then
      call stop_program(info, 'dsyev did not converge, number of non converged elts: ')
    end if

    U(:, 1:k) = A(:, n-k+1:n)
    L(1:k) = W(n-k+1:n)
  end subroutine symmetric_ev

  !> Computes the \a k significant eigenvalues
  !! @param [in, out] L vector of eigenvalues
  !! @param [in] order is an integer giving the order of element in \a L. The following name constant values from \a general_constant are available: \a ASC or \a ASCENDING for ascending order; \a DESC or \a DESCENDING for descending order
  !<
  subroutine ev_signifivant(L, order)
    real(dp), dimension(:), intent(in out) :: L
    integer, intent(in) :: order
    !local variables
    real(dp) :: lsum, sigsum
    integer :: nsig, i, incr, first, last

    lsum = sum(L)
    call debug(lsum, "ev_signifivant: sum(L) = ")
    if(isNaN(lsum)) call debug(L, 'L = ')
    sigsum = 0.0_dp
    nsig = 0
    select case(order)
      case (ASC)
        incr = -1
        first = size(L,1)
        last = 1
      case (DESC)
        incr = 1
        first = 1
        last = size(L,1)
    end select
    i = first
    do while((sigsum/lsum < 0.999_dp).and.(i/=last+incr))
      nsig = nsig + 1
      sigsum = sigsum + L(i)
      i = i+incr
    end do
    !zeroing non significant values
    if(nsig<size(L,1))then
      L(first+nsig*incr:last:incr)  = 0.0_dp
    end if
  end subroutine ev_signifivant

!   !> Computes the increment, the indexes of the smallest and the largest element in a sorted vector
!   !!
!   !! @param [in] order Order of elements in the vector
!   !! @param [in] incr unit increment to move from the smallest to the largest element (+1 for ASC ordered and -1 for descending ordered vector )
!   !! @param [in] order Order of elements in the vector
!   !! @param [in] order Order of elements in the vector
!   !!\details
!   !<
!   subroutine incr_first_last(order, incr, first, last)
!   end subroutine incr_first_last
!
!   !> Computes the number of significant values in a vector
!   interface vec_signifivant_count
!     module procedure vsc_with_nlimit
!     module procedure vsc_with_sum_ratio
!     module procedure vsc_with_bound_ratio
!   end interface vec_signifivant_count
!
!   !> Computes the number of significant values in a sorted vector
!   !!
!   !<
!   subroutine vec_signifivant_count(L, order)
!   end subroutine vec_signifivant_count
end module evd_tools