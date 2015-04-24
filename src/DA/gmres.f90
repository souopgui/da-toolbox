!> @file gmres.f90
!! Defines the GMRES solver for a specialized system
!<

!> @brief GMRES module
!!
!<
module gmres
    use general_constant
    use debug_tools
    use com_tools
implicit none


contains

    !> @brief Solve a given linear system using the GMRES approach
    !! @param [in, out] td_ep exchange parameter data structure
    !! @param [in] obs_std vector of standard deviation of observation error
    !! @param [in] b right hand side of the system
    !! @param [out] x solution of the system
    !! @param [in] matrixmult procedure that applies the matrix of the linear system
    !! @param [in] nItrMax maximum number of iterations of the Arnoldi process
    !! @param [in] tol prescibed tolerance, nomr of the residual
    !! @param [in] nRestart maximum number of restart
    !! @details
    !!  This procedure applies the restarted GMRES method to solve a given
    !!  linear system. The fixed number of iteration for the Arnoldi
    !!  process is nItrMax. And the maximum number of restart is nRestart.
    !!  The precess is restarted until the tolerance is reached on the
    !!  residual or the maximum number of restart is reached.
    !!  obs_std is added for convenience, to have the exact same signature
    !!  as RBCG. Otherwise, it will simply be build into matrixmult
    !<
    subroutine gmres_solve(td_ep, obs_std, b, x, matrixmult, nItrMax, tol, nRestart)
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
        real(dp), dimension(:), intent(in) :: b
        real(dp), dimension(:), intent(in out) :: x
        integer , intent(in) :: nItrMax
        real(dp), intent(in) :: tol
        integer , optional, intent(in) :: nRestart
        ! local variables
        real(dp), dimension(size(b),nItrMax) :: V
        real(dp), dimension(nItrMax,nItrMax) :: H
        real(dp), dimension(size(b)) :: vj, Avj, Ax, v_hat, r, Vy
        real(dp), dimension(nItrMax) :: y, c, s
        real(dp), dimension(nItrMax+1) :: g ! RSH of the triangular system
        ! others
        real(dp) :: r0_norm ! norm of the residue at the zeroth iteration
        real(dp) :: rj_norm ! norm of the residue at the current iteration
        real(dp) :: hij, hipj, hjpj !elements of the last column of H
        real(dp) :: denum, beta, omega
        integer  :: i, j, k, il_nRestart
        logical  :: converged, hv_diverged
        character(*), parameter :: fmt0= "('*** GMRES   Iter ***',I3,' -- |r_0|=',ES14.5E3)"
        character(*), parameter :: fmti= "('*** GMRES * restart', I3,' * Iter ***',I3,' -- |r|/|r_0|=',ES14.5E3)"
        !
        if( present(nRestart) )then
            il_nRestart = nRestart
        else
            il_nRestart = 0
        end if
        !
        call debug(tol, "GMRES - tol", tag=dALLWAYS)
        call debug(nItrMax, "GMRES - nItrMax", tag=dALLWAYS)
        call debug(il_nRestart, "GMRES - nRestart", tag=dALLWAYS)
        ! initialization of local variables
        c = 0.0_dp
        s = 0.0_dp
        !
        x = 0
        !restart loop
        k = 0

        converged = .false.
        !position of the restarting part
        do while( (k<= il_nRestart).and.(.not.converged) )
            if(k == 0)then
                r = b
                r0_norm = sqrt( dot_product(r,r) )
            else
                call matrixmult( td_ep, x, Ax, hv_diverged )
                Ax = Ax + obs_std*x
                r = b - Ax
            end if
            !
            rj_norm = sqrt( dot_product(r,r) )
            vj = r/rj_norm
            !
            j = 0
            g = 0.0_dp
            g(1) = rj_norm
            do while( (j< nItrMax).and.(.not.converged) )
                j = j+1
                V(:,j) = vj
                call matrixmult( td_ep, vj, Avj, hv_diverged )
                Avj = Avj + obs_std*vj
                v_hat = Avj
                do i = 1, j
                    H(i,j) = dot_product( Avj, V(:,i) )
                    v_hat = v_hat - H(i,j)*V(:,i)
                end do
                hjpj = sqrt( dot_product(v_hat,v_hat) )
                !
                if(hjpj>epsilon(hjpj))then
                    vj = v_hat/hjpj
                else
                    ! A^j r is in the Krilov subspace
                    call debug(hjpj, 'hjpj too small')
                    hjpj = 0.0
                end if
                !
                ! Applying the previous plane rotations to the last column of H
                do i = 1, j-1
                    hij      = H(i,j)
                    hipj     = H(i+1,j)
                    H(i,j)   = c(i)*hij - s(i)*hipj
                    H(i+1,j) = s(i)*hij + c(i)*hipj
                end do
                ! Defining the plane rotation to eliminate the subdiagonal of H
                denum = sqrt( H(j,j)**2 + hjpj**2 )
                c(j) = H(j,j)/denum
                s(j) = -hjpj/denum
                ! print*,'j, c, s = ', j, c(j), s(j)
                ! Updating the diagonal element of the last column of H
                H(j,j) = c(j)*H(j,j) - s(j)*hjpj
                ! Applying the plane rotation to the rhs of the triangular system
                ! this defines the entry j+1 and modifies the entry j
                ! Before this the entry j+1 is zero
                g(j+1) = s(j)*g(j)
                g(j) = c(j)*g(j)
                ! Computing the residual
                rj_norm = abs( g(j+1) )
                !Check the convergence
                omega= rj_norm/r0_norm
                if (omega > tol) then
                    write(*,fmt=fmti)k,j,omega
                else
                    !
                    ! Stopping
                    !
                    converged = .true.
                endif
            end do ! Arnoldi loop
            !
            ! Solving the triangular system
            !
            call backward_subst(j, H, g, y)
            !
            ! Computing the solution as linear combination of V
            !
            do i = 1, j
                x = x + V(:,i) * y(i)
            end do
            call matrixmult( td_ep, x, Ax, hv_diverged )
            Ax = Ax + obs_std*x
            rj_norm = sqrt( dot_product(r,r) )
            omega= rj_norm/r0_norm
            converged = ( omega <= tol )
            call debug(omega, 'omega')
            !
            !
            k = k+1
        end do !restart loop
        !
    end subroutine gmres_solve


    !> @brief Solve a upper triangular linear system
    !! @param [in] n size of the system
    !! @param [in] U upper triangular matrix
    !! @param [in] b right hand side vector
    !! @param [out] x solution of theupper triangular system
    !!
    !! @details U, d, x can be larger the size of the system
    !!  It is done this way because the procedure was written to solve
    !!  the triangular system in GMRES, and the size can vary according to
    !!  the requested tolerance
    !!  This procedure declare all those variables as assumed-shape arrays
    !<
    subroutine backward_subst( n, U, b, x )
        implicit none
        integer :: n
        real(dp), dimension(:,:), intent(in)  :: U
        real(dp), dimension(:), intent(in)  :: b
        real(dp), dimension(:), intent(out) :: x
        !local variables
        real :: tmp
        integer :: i,j

        x = 0.0_dp
        do i = n, 1, -1
            tmp = b(i)
            do j = i+1, n
                tmp = tmp - U(i,j)*x(j)
            end do
            x(i) = tmp/U(i,i)
        end do
    end subroutine backward_subst

end module gmres