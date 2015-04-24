!> @file lanczos_solver.f90
!!
!<

!> @brief Fortran module implementing RBLanczos solvers
!!
!<
program lanczos_solver
implicit none
    integer, parameter :: N = 100 !state dimension
    integer, parameter :: M = 100 !obs dimension
    integer, parameter :: MaxNIter = 500 !number of iterations
    real, parameter :: cg_stop = 1.0e-16

    real, dimension(N,N) :: B !background covariance matrix
    real, dimension(M) :: R_cov !observation error covariance matrix
    real, dimension(M) :: R_std !observation error standard deviation
    real, dimension(M,N) :: H !observation operator
    real, dimension(N,N) :: F !model operator
    integer, parameter :: mnproc = 1

    print*, "MaxNIter = ", MaxNIter
    call test_identity()
contains

    subroutine cg_solve(n_state, n_obs, obs_anm, obs_err,  u, iterMax)

!         n,m,nt,mt,l,ls,lz,nr,nq,nest,nobmax,ntyp,
!      &                      ntc,nrvmax,mod_dir,par_dir,io_dir,trj_dir,
!      &                      n_obs,beta,obs_anm,obs_dtg,obs_tme,obs_var,
!      &                      obs_err,obs_xi,obs_yj,obs_lvl,niter,adj_dt,
!      &                      obs_dt,ic_date,ic_time,iec,sigz,rep,weak_const,
!      &                      io_prec,ntime)
!     implicit none
! c
! c declare passed variables (primary variables)
! c
!     integer n,m,nt,mt,l,ls,lz,nr,nq
!     integer nest,nobmax,ntyp,ntc,nrvmax
!     character(256) :: par_dir,io_dir
!     character(256) :: trj_dir,mod_dir
!     character(2)   :: io_prec
!     integer n_obs,ntime
!     real    obs_anm(n_obs)
!     real    beta(   n_obs)
!     integer obs_dtg(n_obs)
!     integer obs_tme(n_obs)
!     integer obs_var(n_obs)
!     real    obs_err(n_obs)
!     real    obs_xi( n_obs)
!     real    obs_yj( n_obs)
!     real    obs_lvl(n_obs)
!     integer niter,ic_date,ic_time
!     real    obs_dt,adj_dt
!     integer iec(8)
!     logical sigz,rep,weak_const


    integer, intent(in) :: n_state, n_obs, iterMax
    real, intent(in out) :: obs_anm(n_obs)
    real, intent(in) :: obs_err(n_obs)
    real, intent(in out) :: u(n_state)

! c
! c declare conjugate gradient variables
! c
    integer cg_itr
    real    cg_rzOld
    real    cg_beta
    real    cg_alpha
    real    cg_omega
    real    cgAp( n_obs)
    real    cgr(  n_obs)
    real    cgp(  n_obs)
    real    cgp1( n_obs)
    real    cg_bsum1,cg_bsum
    real    cg_rz1,cg_rz
    real    cg_rAz1,cg_rAz
    real    cg_rho1,cg_rho

    real    beta(   n_obs)
    real    true_r(  n_obs)
    real    true_r0_norm, true_omega, omega2
! c
! c declare local variables
! c
    real    loc_sum,J1,Jmin
    real    n_obs1,nobs1
    character(256) :: sfx,fname
    integer i,j,k,len,nobs,cond,fdate,ftime,var



    call matrixmult2(beta/obs_err, u)
    call true_obsgap(u*0, true_r)
    true_r0_norm = global_norm( n_obs,true_r )
! c
! c scale innovations by observation error standard deviation
! c
    do i=1,n_obs
       obs_anm(i) = obs_anm(i)/obs_err(i)
       beta(   i) = 0.0
    enddo
    !call costfun(1,0,n_obs,obs_var,obs_err,obs_anm,beta)
! c
! c initialize variables
! c

       beta = 0.0
       do i=1,n_obs
          beta(i) = 0.0
          cgr( i) = obs_anm(i)
       enddo
       cg_bsum1 = sum(cgr*cgr)
       !call xcglobsum(cg_bsum1,cg_bsum)

       cg_bsum = sqrt(cg_bsum1)
       if(mnproc.eq. 1) write(*,*)'CG_SUM',cg_bsum
       cg_itr = 1
       cond   = 1
       cg_omega = 0.0
! c
! c perform cg iterations
! c
    do while(cond.eq. 1 .and.cg_itr.le.itermax)
       cg_rz1 = sum(cgr*cgr)
       !call xcglobsum(cg_rz1,cg_rz)
       cg_rz = cg_rz1
       if(cg_itr.eq. 1) then
          cgp       = cgr
          cg_rzOld  = cg_rz
       else
          cg_beta   = cg_rz/cg_rzOld
          cg_rzOld  = cg_rz
          do i=1,n_obs
             cgp(i) = cgr(i) + cg_beta*cgp(i)
          enddo
       endif
! c
! c perform matrix multiplication
! c
       cgAp=0.0; cgp1=0.0
       cgp1=cgp/obs_err !! pre-scale by observation error standard deviation

       call matrixmult(cgp1,cgAp)

!        call matrixmult(1,rep,n,m,nt,mt,l,ls,lz,nr,nq,nest,nobmax,
!      &                     ntyp,ntc,nrvmax,mod_dir,par_dir,io_dir,
!      &                     trj_dir,n_obs,cgp1,cgAp,obs_dtg,obs_tme,
!      &                     obs_var,obs_err,obs_xi,obs_yj,obs_lvl,
!      &                     niter,niter,adj_dt,obs_dt,ic_date,ic_time,
!      &                     iec,weak_const,io_prec,ntime)
       cgAp=cgAp/obs_err !! post-scale by observation error standard deviation
!        if(sigz) then
!           lz = 0
!        endif
! c
! c perform I acting on y (and add to Ry)
! c
       do i=1,n_obs
          cgAp(i) = cgAp(i) + cgp(i)
       enddo
! c
! c perform minimization steps
! c
       cg_rAz1  = sum(cgp*cgAp)
       !call xcglobsum(cg_rAz1,cg_rAz)
       cg_rAz = cg_rAz1
       cg_alpha = cg_rz/cg_rAz
       do i=1,n_obs
          beta(i) = beta(i) + cg_alpha*cgp(i)
          cgr(i)  = cgr(i)  - cg_alpha*cgAp(i)
       enddo
       cg_rho1  = sum(cgr*cgr)
       !call xcglobsum(cg_rho1,cg_rho)
       cg_rho = cg_rho1

       cg_omega = sqrt(cg_rho)/cg_bsum
! c
! c save values for restart
! c
!        write(sfx,'(''CG_rstart_proc'',i4.4,''.dat'')')mnproc
!        fname = trim(trj_dir) // '/' // trim(sfx)
!        len   = len_trim(fname)
!        open( 10,file=fname(1:len),form='unformatted')
!        write(10)n_obs,cg_itr
!        write(10)cg_rzOld,cg_bsum
!        write(10)(beta(i),i=1,n_obs)
!        write(10)(cgr( i),i=1,n_obs)
!        write(10)(cgp( i),i=1,n_obs)
!        close(10)
! c
! c check for convergence
! c
        call matrixmult2(beta/obs_err, u)
        call true_obsgap(u, true_r)
        true_omega = sqrt( dot_product(true_r,true_r) )/true_r0_norm
        omega2 = sqrt(dot_product(cgr*obs_err,cgr*obs_err))/true_r0_norm
       if(cg_omega .gt. cg_stop) then
          if(mnproc.eq. 1) then
            write(*,*)'*** CG Iteration   ***',cg_itr,'rrtol value',cg_omega, 't2', omega2, 'tt',true_omega
          endif
          !call costfun(2,cg_itr,n_obs,obs_var,obs_err,obs_anm,beta)
          cg_itr = cg_itr + 1
       else
          cond = 0
       endif
! c
! c end cg itertion loop
! c
    enddo
    if(mnproc.eq. 1) then
       write(*,*)'*** CG Convergence ***',cg_itr,'rrtol value',cg_omega
    endif
    !call costfun(2,cg_itr,n_obs,obs_var,obs_err,obs_anm,beta)
! c
! c recover beta by scaling by observation error standard deviation
! c
    do i=1,n_obs
       beta(i)    = beta(i)/obs_err(i)
       obs_anm(i) = obs_anm(i)*obs_err(i) !! undo scaling for proper computation of J1 and Jmin
    enddo
! c
! c perform final sweep...
! c
    call matrixmult2(beta, u)

!     call matrixmult(2,rep,n,m,nt,mt,l,ls,lz,nr,nq,nest,nobmax,
!      &                  ntyp,ntc,nrvmax,mod_dir,par_dir,io_dir,trj_dir,
!      &                  n_obs,beta,obs_anm,obs_dtg,obs_tme,obs_var,
!      &                  obs_err,obs_xi,obs_yj,obs_lvl,niter,niter,
!      &                  adj_dt,obs_dt,ic_date,ic_time,iec,weak_const,
!      &                  io_prec,ntime)

    !end subroutine cg_solve
    end subroutine cg_solve

    !> @brief Lanczos algorithm with re-orthogonalization for 4D-VAR
    !!
    !<
    subroutine lanczos_solve(n_state, n_obs, obs_anm, obs_err,  u, iterMax)
    implicit none
        integer, intent(in) :: n_state, n_obs, iterMax
        real, intent(in) :: obs_anm(n_obs)
        real, intent(in) :: obs_err(n_obs) ! standard deviation of obs error
        real, intent(in out) :: u(n_state)
        !local variables
        real :: rbl_r(n_obs)
        real :: rbl_lambda(n_obs)
        real :: rbl_alpha
        real :: rbl_beta0
        real :: rbl_beta
        real :: rbl_betap1 !rbl_beta at rbl_itr+1 iteration
        real :: rbl_zeta
        real :: rbl_vm1(n_obs) !v at rbl_itr-1 iteration
        real :: rbl_v(n_obs) !
        real :: rbl_z(n_obs)
        real :: rbl_w(n_obs)
        real :: rbl_tmp1(n_obs)
        real :: rbl_tmp2(n_obs)
        !matrix variables
        real :: rbl_Vmat(n_obs,iterMax)
        !tridiagonal solve variables
        real :: T_diag(iterMax)
        real :: T_off(iterMax)
        real :: rbl_d(iterMax)
        real :: rbl_s(iterMax)
        real :: rbl_tmp(n_obs)
        !diagnostic variables
        real :: ri_norm
        real :: r0_norm
        real :: rbl_omega
        !real loc_sum
        real    true_r(  n_obs)
        real    true_r0_norm, true_omega, omega2

        integer :: rbl_itr, rbl_nIter, rbl_oi
        real    :: rbl_zv
        logical :: converged
        character(*), parameter :: fmt0=&
              "('*** RBLanczos-O with re-orthogonalization r0_norm',ES14.5E3)"
        character(*), parameter :: fmti=&
              "('*** RBLanczos-O Iteration ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"
        character(*), parameter :: fmtl=&
              "('*** RBLanczos-O Converged ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"


        call true_obsgap(u*0, true_r)
        true_r0_norm = global_norm( n_obs,true_r )

        !rbl_itr = 0
        rbl_r = obs_anm/obs_err
        !diagnostic
        r0_norm = global_norm( n_obs,rbl_r )
        if(mnproc.eq. 1) then
            write(*,fmt=fmt0) r0_norm
        end if
        !end of diagnostic
        rbl_beta0 = global_norm(n_obs, rbl_r)
        !rbl_itr=1
        rbl_v = rbl_r/rbl_beta0
        rbl_beta   = 0.0
        rbl_vm1    = 0.0
        rbl_Vmat(:,1) = rbl_v

        T_diag = 0.0
        T_off = 0.0

        converged = .false.
        rbl_itr = 1
        do while( (rbl_itr<=iterMax).and.(.not.converged) )
            if(rbl_itr>1) rbl_beta = rbl_betap1

            rbl_tmp1 = rbl_v/obs_err
            call matrixmult(rbl_tmp1, rbl_tmp2)
            rbl_w = (rbl_tmp2/obs_err + rbl_v) - rbl_beta*rbl_vm1

            rbl_vm1  = rbl_v !updating vm1
            rbl_alpha = global_dot_product(n_obs, rbl_w, rbl_v )
            rbl_w = rbl_w - rbl_alpha*rbl_v

            !Re-orthogonalize rbl_w using columns of V and Z
            do rbl_oi = 1, rbl_itr-1
                rbl_zv = global_dot_product(n_obs, rbl_Vmat(:,rbl_oi), rbl_w )
                rbl_w = rbl_w - rbl_zv*rbl_Vmat(:,rbl_oi)
            end do
            !End of re-orthogonalization

            rbl_betap1 = global_norm(n_obs, rbl_w)
            if( rbl_betap1>epsilon(1.0) )then
                rbl_v    = rbl_w/rbl_betap1
                if(rbl_itr<iterMax)then
                    rbl_Vmat(:,rbl_itr+1) = rbl_v
                end if
            else
                converged = .true.
                !print*, "rbl_beta = ",  rbl_beta
                !print*, "rbl_itr = ",  rbl_itr
            end if
            T_diag(rbl_itr) = rbl_alpha
            if(rbl_itr>1)then
                T_off(rbl_itr-1) = rbl_beta
            end if
            !diagnostic
                rbl_d=0.0
                rbl_d(1) = rbl_beta0
                call tridiag_solve( rbl_itr, T_off, T_diag, T_off, rbl_d, rbl_s )
                rbl_r = - rbl_betap1*rbl_s(rbl_itr) *rbl_v !ei^t s
                !multiply by obs_err to get similar result to CG
                ri_norm = global_norm(n_obs,rbl_r)
                rbl_omega = ri_norm/r0_norm

                rbl_lambda =  matmul( rbl_Vmat(:, 1:rbl_itr), rbl_s(1:rbl_itr) )

                !
                ! check for convergence
                !
                rbl_lambda = rbl_lambda/obs_err
                call matrixmult2(rbl_lambda, u)
                call true_obsgap(u, true_r)
                true_omega = global_norm(n_obs,true_r)/true_r0_norm
                omega2 = global_norm(n_obs,rbl_r*obs_err)/true_r0_norm
                if(rbl_omega .gt. cg_stop) then
                    if(mnproc.eq. 1) then
                        write(*,fmt=fmti)rbl_itr,rbl_omega,omega2,true_omega
                    endif
                    !call costfun(2,cg_itr,n_obs,obs_var,obs_err,obs_anm,beta)
                else
                    converged = .true.
                endif
            !end of diagnostic
            rbl_itr = rbl_itr+ 1
        end do !rblanczos iteration

        rbl_nIter = min(rbl_itr-1, iterMax)

        if(mnproc.eq. 1) then
            write(*,fmt=fmtl)rbl_nIter,rbl_omega,omega2,true_omega
        endif
        call matrixmult2( rbl_lambda, u )
    end subroutine lanczos_solve

    !> @brief RBLanczos algorithm for 4DVAR
    subroutine rblw_solve(n_state, n_obs, obs_anm, obs_err,  u, iterMax)
    implicit none
        integer, intent(in) :: n_state, n_obs, iterMax
        real, intent(in) :: obs_anm(n_obs)
        real, intent(in) :: obs_err(n_obs) ! standard deviation of obs error
        real, intent(in out) :: u(n_state)
        !local variables
        real :: rbl_r(n_obs)
        real :: rbl_lambda(n_obs)
        real :: rbl_t(n_obs)
        real :: rbl_beta0
        real :: rbl_beta
        real :: rbl_betap1 !rbl_beta at rbl_itr+1 iteration
        real :: rbl_zeta
        real :: rbl_vm1(n_obs) !v at i-1 iteration
        real :: rbl_v(n_obs) !
        real :: rbl_z(n_obs)
        real :: rbl_gamma
        real :: rbl_pm1(n_obs) !p at iteration i-1
        real :: rbl_ptmp(n_obs)
        real :: rbl_p(n_obs)
        real :: rbl_alpha
        real :: rbl_eta
        real :: rbl_w(n_obs)
        real :: rbl_q(n_obs)
        !diagnostic variables
        real :: T_diag(iterMax)
        real :: T_off(iterMax)
        real :: d(iterMax)
        real :: s(iterMax)
        real :: tmp(n_obs)
        real :: r0_norm
        real :: ri_norm
        real :: rbl_omega
        real    true_r(  n_obs)
        real    true_r0_norm, true_omega, omega2


        integer :: rbl_itr, rbl_nIter
        logical :: converged
        character(*), parameter :: fmt0=&
              "('*** RBLanczos-W without ro-thogonalization r0_norm',ES14.5E3)"
        character(*), parameter :: fmti=&
              "('*** RBLanczos Iteration ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"
        character(*), parameter :: fmtl=&
              "('*** RBLanczos Converged ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"

        call true_obsgap(u*0, true_r)
        true_r0_norm = sqrt( dot_product(true_r,true_r) )

        !rbl_itr = 0
        rbl_r = obs_anm/obs_err**2
        !diagnostic
            r0_norm = sqrt( dot_product(rbl_r, rbl_r) )
            if(mnproc.eq. 1) then
                write(*,fmt=fmt0)r0_norm
            end if
        !end of diagnostic
        rbl_lambda = 0.0
        call matrixmult(rbl_r, rbl_t)
        rbl_beta0 = sqrt( dot_product(rbl_t, rbl_r) )
        !rbl_itr=1
        rbl_zeta = rbl_beta0
        rbl_v = rbl_r/rbl_beta0
        rbl_z = rbl_t/rbl_beta0
        rbl_beta = 0.0
        rbl_vm1 = 0.0
        rbl_gamma = 0.0
        rbl_pm1 = 0.0
        !initialization of diagnostic variables
        T_diag = 0.0
        T_off = 0.0
        d = 0.0
        d(1) = rbl_beta0
        !End of initialization of diagnostic variables
        converged = .false.
        rbl_itr = 1
        do while( (rbl_itr<= iterMax).and.(.not.converged) )
            !print*, "rbl_itr = ", rbl_itr, "beta = ", rbl_beta
            if(rbl_itr>1) rbl_beta = rbl_betap1

            rbl_q = (rbl_v + rbl_z/obs_err**2) - rbl_beta*rbl_vm1
            rbl_vm1  = rbl_v !updating vm1
            rbl_alpha = dot_product( rbl_q, rbl_z )
            if ( rbl_itr > 1 ) then
                rbl_gamma = rbl_beta/rbl_eta
                rbl_zeta  = -rbl_gamma * rbl_zeta
            end if
            rbl_eta  = rbl_alpha - rbl_gamma * rbl_beta
            rbl_ptmp = rbl_p !saving previous p
            rbl_p    = (rbl_v - rbl_beta*rbl_p)/rbl_eta
            rbl_pm1  = rbl_ptmp !updating pm1
            rbl_lambda = rbl_lambda + rbl_zeta * rbl_p
            rbl_w      = rbl_q - rbl_alpha * rbl_v
            call matrixmult(rbl_w, rbl_t)
            !rbl_itr+1 varibles in ith loop
            rbl_betap1 = sqrt( dot_product(rbl_t, rbl_w) )
            if( rbl_betap1>epsilon(1.0) )then
                rbl_v    = rbl_w/rbl_betap1
                rbl_z    = rbl_t/rbl_betap1
            else
                converged = .true.
            end if

            !diagnostic
                T_diag(rbl_itr) = rbl_alpha
                if(rbl_itr>1)then
                    T_off(rbl_itr-1) = rbl_beta
                end if
                d=0.0
                d(1) = rbl_beta0
                call tridiag_solve( rbl_itr, T_off, T_diag, T_off, d, s )
                rbl_r = - rbl_betap1*s(rbl_itr) *rbl_v !ei^t s
                !ultiply by obs_err to get the same result as in CG
                ri_norm = global_norm(n_obs,rbl_r)
                rbl_omega = ri_norm/r0_norm
                !
                ! check for convergence
                !
                call matrixmult2(rbl_lambda, u)
                call true_obsgap(u, true_r)
                true_omega = global_norm(n_obs,true_r)/true_r0_norm
                omega2 = global_norm(n_obs,rbl_r*obs_err**2)/true_r0_norm
                if(rbl_omega .gt. cg_stop) then
                    if(mnproc.eq. 1) then
                        write(*,fmt=fmti)rbl_itr,rbl_omega,omega2,true_omega
                    endif
                    !call costfun(2,cg_itr,n_obs,obs_var,obs_err,obs_anm,beta)
                else
                    converged = .true.
                endif
            !end of diagnostic

            rbl_itr=rbl_itr+1
        end do
        rbl_nIter = min(rbl_itr-1, iterMax)
        if(mnproc.eq. 1) then
            write(*,fmt=fmtl)rbl_nIter,rbl_omega,omega2,true_omega
        endif

        print*, "rbl_nIter = ",  rbl_nIter
        call matrixmult2(rbl_lambda, u)
    end subroutine rblw_solve

    !> @brief RBLanczos algorithm with re-orthogonalization for 4D-VAR
    !!
    !<
    subroutine rblo_solve(n_state, n_obs, obs_anm, obs_err,  u, iterMax)
    implicit none
        integer, intent(in) :: n_state, n_obs, iterMax
        real, intent(in) :: obs_anm(n_obs)
        real, intent(in) :: obs_err(n_obs) ! standard deviation of obs error
        real, intent(in out) :: u(n_state)
        !local variables
        real :: rbl_r(n_obs)
        real :: rbl_lambda(n_obs)
        real :: rbl_t(n_obs)
        real :: rbl_beta0
        real :: rbl_beta
        real :: rbl_betap1 !rbl_beta at rbl_itr+1 iteration
        real :: rbl_zeta
        real :: rbl_vm1(n_obs) !v at rbl_itr-1 iteration
        real :: rbl_v(n_obs) !
        real :: rbl_z(n_obs)
        real :: rbl_gamma
        real :: rbl_pm1(n_obs) !p at iteration rbl_itr-1
        real :: rbl_ptmp(n_obs)
        real :: rbl_p(n_obs)
        real :: rbl_alpha
        real :: rbl_eta
        real :: rbl_w(n_obs)
        real :: rbl_q(n_obs)
        !matrix variables
        real :: V(n_obs,iterMax)
        real :: Z(n_obs,iterMax)
        !tridiagonal solve variables
        real :: T_diag(iterMax)
        real :: T_off(iterMax)
        real :: d(iterMax)
        real :: s(iterMax)
        !diagnostic variables
        real :: ri_norm
        real :: r0_norm
        real :: rbl_omega
        real :: J0, J, Jb, J_obs
        !real loc_sum
        real    true_r(  n_obs)
        real    true_r0_norm, true_omega, omega2

        integer :: rbl_itr, rbl_nIter
        logical :: converged
        character(*), parameter :: fmt0=&
              "('*** RBLanczos-O with re-orthogonalization r0_norm',ES14.5E3)"
        character(*), parameter :: fmti=&
              "('*** RBLanczos-O Iteration ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"
        character(*), parameter :: fmtl=&
              "('*** RBLanczos-O Converged ***',I3,' * rrtol value'"//&
              ",ES14.5E3,' * t2',ES14.5E3,' * tt',ES14.5E3)"


        call true_obsgap(u*0, true_r)
        true_r0_norm = global_norm(n_obs,true_r)

        !rbl_itr = 0
        rbl_r = obs_anm/obs_err**2
        !diagnostic
        !multiply by obs_err to get similar result to CG
        r0_norm = global_norm(n_obs, rbl_r)
        if(mnproc.eq. 1) then
            write(*,fmt=fmt0) r0_norm
        end if
        !end of diagnostic
        call matrixmult(rbl_r, rbl_t)
        rbl_beta0 = global_dot_product(n_obs, rbl_t, rbl_r)
        rbl_beta0 = sqrt(rbl_beta0)
        !rbl_itr=1
        rbl_v = rbl_r/rbl_beta0
        rbl_z = rbl_t/rbl_beta0
        rbl_beta   = 0.0
        rbl_vm1    = 0.0
        Z(:,1) = rbl_z
        V(:,1) = rbl_v

        T_diag = 0.0
        T_off = 0.0

        converged = .false.
        rbl_itr = 1
        do while( (rbl_itr<=iterMax).and.(.not.converged) )
            if(rbl_itr>1) rbl_beta = rbl_betap1
            !print*, "rbl_itr = ", rbl_itr
            !print*, "iterMax = ", iterMax
            rbl_q = (rbl_v + rbl_z/obs_err**2) - rbl_beta*rbl_vm1
            rbl_vm1  = rbl_v !updating vm1
            rbl_alpha = global_dot_product(n_obs, rbl_q, rbl_z )
            !print*, "rbl_q(1:5 ) = ",  rbl_q(1:5 )
            !print*, "rbl_z(1:5 ) = ",  rbl_z(1:5 )
            !print*, "rbl_v(1:5 ) = ",  rbl_v(1:5 )
            !print*, "rbl_beta0 = ",  rbl_beta0
            rbl_w = rbl_q - rbl_alpha*rbl_v
            !Re-orthogonalize rbl_w using columns of V and Z
            call rbl_reorthogonalization( V(:,1:rbl_itr), Z(:,1:rbl_itr), rbl_w )
            !End of re-orthogonalization
            call matrixmult(rbl_w, rbl_t)
            rbl_betap1 = global_dot_product(n_obs, rbl_t, rbl_w)
            rbl_betap1 = sqrt(rbl_betap1)
            !print*, "rbl_itr = ", rbl_itr, "beta = ", rbl_beta
            if( rbl_betap1>epsilon(1.0) )then
                rbl_v    = rbl_w/rbl_betap1
                rbl_z    = rbl_t/rbl_betap1
                if(rbl_itr<iterMax)then
                    Z(:,rbl_itr+1) = rbl_z
                    V(:,rbl_itr+1) = rbl_v
                end if
            else
                converged = .true.
                !print*, "rbl_beta = ",  rbl_beta
                !print*, "rbl_itr = ",  rbl_itr
            end if
            T_diag(rbl_itr) = rbl_alpha
            if(rbl_itr>1)then
                T_off(rbl_itr-1) = rbl_beta
            end if
            !diagnostic
                d=0.0
                d(1) = rbl_beta0
                call tridiag_solve( rbl_itr, T_off, T_diag, T_off, d, s )
                rbl_r = - rbl_betap1*s(rbl_itr) *rbl_v !ei^t s
                !multiply by obs_err to get similar result to CG
                ri_norm = global_norm(n_obs,rbl_r)
                rbl_omega = ri_norm/r0_norm

                rbl_lambda =  matmul( V(:, 1:rbl_itr), s(1:rbl_itr) )

                !
                ! check for convergence
                !

                call matrixmult2(rbl_lambda, u)
                call true_obsgap(u, true_r)
                true_omega = global_norm(n_obs,true_r)/true_r0_norm
                omega2 = global_norm(n_obs,rbl_r*obs_err**2)/true_r0_norm
                if(rbl_omega .gt. cg_stop) then
                    if(mnproc.eq. 1) then
                        write(*,fmt=fmti)rbl_itr,rbl_omega,omega2,true_omega
                    endif
                    !call costfun(2,cg_itr,n_obs,obs_var,obs_err,obs_anm,beta)
                else
                    converged = .true.
                endif
            !end of diagnostic
            rbl_itr = rbl_itr+ 1
        end do !rblanczos iteration

        rbl_nIter = min(rbl_itr-1, iterMax)

        if(mnproc.eq. 1) then
            write(*,fmt=fmtl)rbl_nIter,rbl_omega,omega2,true_omega
        endif

        call matrixmult2( rbl_lambda, u )
    end subroutine rblo_solve

    !> @brief re-orthogonalization in RBLanczos
    !! @param[in] V V matrix in the RBLanczos algorithm
    !! @param[in] Z Z matrix in the RBLanczos algorithm
    !! @param[in] w Lanczos vector of the current iteration in RBLanczos
    !!
    !! @details This algorith is based on eq. 34 and 35 in Gurol et al.
    !<
    subroutine rbl_reorthogonalization(V, Z, w)
    implicit none
        real, dimension(:,:), intent(in) :: V, Z
        real, dimension(:), intent(in out) :: w
        !local variables
        integer :: j

        do j = 1, size(V,2)
            w = w - dot_product( Z(:,j), w )*V(:,j)
        end do
    end subroutine rbl_reorthogonalization

    !> @brief Solve a tridiagonal systems
    !! @param[in] n size of vectors
    !! @param[in] a vector of subdiagonal elements
    !! @param[in] b_in vector of diagonal elements
    !! @param[in] c vector of updiagonal elements
    !! @param[in] d_in right hand side vector
    !! @param[out] x solution of tridiagonal system
    !!
    !! @details If @a a and @a c have the same size as @a b, the usefull element
    !!   must be packed in the first part of the vectors and the last element
    !!   be unused
    !<
    subroutine tridiag_solve(n, a, b_in, c, d_in, x)
    implicit none
        integer, intent(in) :: n
        real, dimension(:), intent(in) :: a, b_in, c, d_in
        real, dimension(:), intent(in out) :: x
        !local variables
        real, dimension(n) :: b, d
        real :: t
        integer :: k

        !copy
        b = b_in(1:n)
        d = d_in(1:n)

        !forward elimination process
        do k = 2, n
            t = a(k-1)/b(k-1)
            b(k) = b(k) - t*c(k-1)
            d(k) = d(k) - t*d(k-1)
        end do

        !Backward substitution phase

        x(n) = d(n)/b(n)
        do k = n-1, 1, -1
            x(k) = ( d(k) - c(k)*x(k+1) )/b(k)
        end do
    end subroutine tridiag_solve

    subroutine truth(x)
        real, dimension(M), intent(in out) :: x

        x = 1.0
    end subroutine truth

    subroutine true_obsgap(x, r)
        real, dimension(N), intent(in) :: x
        real, dimension(M), intent(in out) :: r
        !local variables
        real, dimension(N) :: xt
        real, dimension(M) :: b, y

        call truth(xt)
        call HF(xt, b)
        call HF(x, y)
        r = b - y
    end subroutine true_obsgap

    !y = (HF)B(HF)^t x =  HFBF^tH^tx
    subroutine matrixmult( x, y )
    implicit none
        real, dimension(M), intent(in) :: x
        real, dimension(M), intent(out) :: y
        !local variables
        real, dimension(N) :: tmp1, tmp2, tmp3, tmp4

        tmp1 = matmul( transpose(H), x ) !H^tx
        tmp2 = matmul( transpose(F), tmp1 ) !F^tH^tx
        tmp3 = matmul( B, tmp2 ) !BF^tH^tx
        tmp4 = matmul( F, tmp3 ) !FBF^tH^tx
        y = matmul( H, tmp4 )    !HFBF^tH^tx

    end subroutine matrixmult

    subroutine matrixmult2( x, y )
    implicit none
        real, dimension(M), intent(in) :: x
        real, dimension(N), intent(out) :: y
        !local variables
        real, dimension(N) :: tmp1, tmp2

        !print*, "matrixmult2, shape(x)", shape(x)
        !print*, "matrixmult2, shape(y)", shape(y)
        !print*, "matrixmult2, N", N
        tmp1 = matmul( transpose(H), x ) !H^tx
        !print*, "calling matmul 2"
        tmp2 = matmul( transpose(F), tmp1 ) !F^tH^tx
        !print*, "calling matmul 3"
        !print*, "matrixmult2, shape(B)", shape(B)
        y = matmul( B, tmp2 ) !BF^tH^tx

        !print*, "matrixmult2 done"
    end subroutine matrixmult2

    subroutine test_identity()
    implicit none
        integer :: i, j, l

        !truth state and observations
        real, dimension(N) :: x0t, xtt
        real, dimension(M) :: yot
        !forecasted state and observations
        real, dimension(N) :: x0, xt
        real, dimension(M) :: yo, d
        !analysed state and observations
        real, dimension(N) :: x0_rblw, x0_rblo, x0_lanczos, x0_cg, xta, u
        real, dimension(M) :: yoa, y_rblw, y_rblo, y_lanczos, y_cg, res

        F = 0.0
        H = 0.0
        !R = 0.0
        B = 0.0
        do i = 1, N
            F(i,i) = 1.0
            B(i,i) = 10.0
        end do
        do i = 1, M
            H(i,i) = 1.0
        end do
        call random_number(R_cov)
        R_cov = R_cov**2/sqrt(10.0) !
        !R_cov = 1.0d-1
        R_std = sqrt(R_cov)

        call random_number(F)
        F = matmul(F,transpose(F))
        do i = 1, N
            F(i,i) = 2*sum(abs(F(i,:) ) )
        end do

        print*, global_dot_product(N, x0, xt)

        call truth(x0t)
        call HF(x0t, yot)

        x0 = 0
        call HF(x0, yo)
        res = yot - yoa
        print*, "res = ", res
        print*, "residual norm = ", sqrt(dot_product(res, res))
        !read(*,*)

        print*, "--------------------------------------"
        print*, "Without re-orthogonalization"
        print*, "--------------------------------------"
        call truth(x0t)
        call HF(x0t, yot)

        x0 = 0
        call HF(x0, yo)
        d = yot
        !call rbl_solve_saad(N, M, d, R,  x0a, MaxNIter)
        call rblw_solve(N, M, d, R_std,  x0_rblw, MaxNIter)

        call HF(x0_rblw, y_rblw)
        res = y_rblw - yot
        print*, "residual norm = ", sqrt(dot_product(res, res))


        print*, "--------------------------------------"
        print*, "With re-orthogonalization"
        print*, "--------------------------------------"
        call truth(x0t)
        call HF(x0t, yot)

        x0 = 0
        call HF(x0, yo)
        d = yot
        call rblo_solve(N, M, d, R_std,  x0_rblo, MaxNIter)
        call HF(x0_rblo, y_rblo)
        res = y_rblo - yot
        print*, "residual norm = ", sqrt(dot_product(res, res))


        print*, "--------------------------------------"
        print*, "Lanczos with splitted preconditioner"
        print*, "--------------------------------------"
        call truth(x0t)
        call HF(x0t, yot)

        x0 = 0
        call HF(x0, yo)
        d = yot
        call lanczos_solve(N, M, d, R_std,  x0_lanczos, MaxNIter)
        call HF(x0_lanczos, y_lanczos)
        res = y_lanczos - yot
        print*, "residual norm = ", sqrt(dot_product(res, res))


        print*, "--------------------------------------"
        print*, "Conjugate gradient"
        print*, "--------------------------------------"
        x0t = 1.0
        call HF(x0t, yot)


        x0 = 0
        call HF(x0, yo)
        d = yot
        call cg_solve(N, M, d, R_std,  x0_cg, MaxNIter)
        call HF(x0_cg, y_cg)
        res = y_cg - yot
        print*, "residual norm = ", sqrt(dot_product(res, res))


        print*,''
        print*,'Comparison RBLO Lanczos'
        res = y_lanczos-y_rblo
        print*, "Difference norm = ", sqrt(dot_product(res, res))
        print*,''
        print*,'Comparison RBLW RBLO'
        res = y_rblw-y_rblo
        print*, "Difference norm = ", sqrt(dot_product(res, res))
        print*,''
        print*,'Comparison CG RBLO'
        res = y_cg-y_rblo
        print*, "Difference norm = ", sqrt(dot_product(res, res))
        print*,''
        print*,'Comparison CG RBLW'
        res = y_cg-y_rblw
        print*, "Difference norm = ", sqrt(dot_product(res, res))
    end subroutine test_identity

    subroutine model_integration(x0, xt)
        real, dimension(N), intent(in) :: x0
        real, dimension(N), intent(out) :: xt

        xt =  matmul(F, x0)
    end subroutine model_integration

    subroutine sample_obs(x, yo)
        real, dimension(N), intent(in) :: x
        real, dimension(M), intent(out) :: yo

        yo =  matmul(H, x)
    end subroutine sample_obs

    !integrate the model and sample observations
    subroutine HF(x0, yo)
        real, dimension(N), intent(in) :: x0
        real, dimension(M), intent(out) :: yo

        yo =  matmul( H, matmul(F, x0) )
    end subroutine HF


        real function global_dot_product(n, a, b)
        implicit none
            integer n
            real a(n)
            real b(n)
            !local variables
            real loc_dot, glob_dot

            global_dot_product = dot_product(a,b)
        return
        end function global_dot_product


        real function global_norm(n, a)
            implicit none
            integer n
            real a(n)
            !local variables

            global_norm = sqrt(dot_product(a,a))
            return
        end function global_norm

end program lanczos_solver
