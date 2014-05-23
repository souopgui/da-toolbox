!> module for noise (error) addition
!! used for Data Assimilation purpose
!!
!<
PROGRAM rand_test
  USE randmod
  USE random
  USE general_tools
  USE debug_tools
IMPLICIT NONE
  INTEGER, PARAMETER :: IMAX=8
  INTEGER, DIMENSION(IMAX) :: ila_vec
  REAL(KIND=cp), DIMENSION(IMAX) :: rla_vec
  REAL(KIND=cp), DIMENSION(1000) :: rla_rand

	PRINT*, "normal random generator "
	!(mu=0.0_cp, sigma=0.1_cp, radius=0.5_cp, nright=100, stat=.TRUE.)
	CALL test_grand(stat=.TRUE.)
	PRINT*, "extern random generator "
	CALL test_extern_nrand()
!   CALL rand_normal_vector(rla_rand, stat=.TRUE.)
!   CALL write_vector_for_plot(rla_rand, "rt_rand.dat")

  CALL finilize_normal_rand()


CONTAINS

  !> \brief Test the Gaussian random numbers generator
  !! \param[in] mu mean of the distribution
  !! \param[in] sigma Standard deviation of the distribution
  !! \param[in] radius Distance from the mean to the maximum/minimum number that can be generated, default is set to 4*rd_sigma
  !! \param[in] nright Number of possible values greater than the mean, default is set to 127
  !! \param[in] stat says if it is necessary to save statistics
  !! \details the total number of possible values that can be generated is 2*id_nright + 1
  !<
	SUBROUTINE test_grand(mu, sigma, radius, nright, stat)
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: mu, sigma
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: radius
    INTEGER, OPTIONAL, INTENT(IN) :: nright
    LOGICAL, OPTIONAL, INTENT(IN) :: stat
    !local var
		INTEGER :: ibi, clock1, clock2
		REAL(KIND=cp) :: rl_rand

		CALL init_normal_rand(mu=mu, sigma=sigma, radius=radius, nright=nright, stat=stat)
		CALL SYSTEM_CLOCK(COUNT=clock1)
		DO ibi = 1, 10000000
			CALL normal_rand(rl_rand)
		END DO
		CALL SYSTEM_CLOCK(COUNT=clock2)
		CALL debug(clock2-clock1, 'Total time = ', tag=dALLWAYS)
		CALL save_grand_stat()
		CALL finilize_normal_rand()
	END SUBROUTINE test_grand

	!> \brief Test the extern normal random numbers generator
	SUBROUTINE test_extern_nrand()
    !local var
		INTEGER :: ibi, clock1, clock2
		REAL(KIND=cp) :: rl_rand
		TYPE(rng_environment) :: tl_rnge

		CALL init_normal_rand(stat=.TRUE.)
		CALL copy_normal_rnge(tl_rnge)
		CALL finilize_normal_rand()

		CALL SYSTEM_CLOCK(COUNT=clock1)
		DO ibi = 1, 10000000
			rl_rand = random_normal()
			!CALL update_rnge_stat(tl_rnge, rl_rand)
		END DO
		CALL SYSTEM_CLOCK(COUNT=clock2)
		CALL debug(clock2-clock1, 'Total time = ', tag=dALLWAYS)
		CALL save_extern_nstat(tl_rnge)
		CALL finilize_rng_environment(tl_rnge)
	END SUBROUTINE test_extern_nrand

  !> \brief generate a vector of Gaussian random numbers. Each element of the vector is an independant Gaussian with mean mu and std deviation sigma
  !! \param[out] rda_rand vector that contains the generated numbers
  !! \param[in] mu mean of the distribution
  !! \param[in] sigma Standard deviation of the distribution
  !! \param[in] radius Distance from the mean to the maximum/minimum number that can be generated, default is set to 4*rd_sigma
  !! \param[in] nright Number of possible values greater than the mean, default is set to 127
  !! \param[in] stat says if it is necessary to save statistics
  !! \details the total number of possible values that can be generated is 2*id_nright + 1
  !<
	SUBROUTINE rand_normal_vector(rda_rand, mu, sigma, radius, nright, stat)
		REAL(KIND=cp), DIMENSION(:) :: rda_rand
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: mu, sigma
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: radius
    INTEGER, OPTIONAL, INTENT(IN) :: nright
    LOGICAL, OPTIONAL, INTENT(IN) :: stat
    !local var
    INTEGER :: ibi

    CALL init_normal_rand(mu=mu, sigma=sigma, radius=radius, nright=nright, stat=stat)
    DO ibi = 1, SIZE(rda_rand)
			CALL init_random_seed()
			CALL normal_rand( rda_rand(ibi) )
    END DO
    CALL finilize_normal_rand()
	END SUBROUTINE rand_normal_vector

	SUBROUTINE save_grand_stat()
		INTEGER :: ibi, il_size
		REAL(cp), DIMENSION( : ), ALLOCATABLE :: rla_x, rla_f
		!getting and saving stats
		CALL get_normal_stat(ssize=il_size)
		ALLOCATE( rla_x(il_size), rla_f(il_size) )
		CALL get_normal_stat( x = rla_x, y = rla_f )
		CALL write_vector_for_plot( "rt_grand_stat.dat", rla_f, rla_x )
		CALL get_normal_stat( x = rla_x, pdf = rla_f )
		CALL write_vector_for_plot( "rt_grand_pdf.dat", rla_f, rla_x )
	END SUBROUTINE save_grand_stat

	SUBROUTINE save_extern_nstat(td_rnge)
		TYPE(rng_environment), INTENT(IN) :: td_rnge
		INTEGER :: ibi, il_size
		REAL(cp), DIMENSION( : ), ALLOCATABLE :: rla_x, rla_f
		!getting and saving stats
		CALL get_rng_stat(td_rnge, ssize=il_size)
		ALLOCATE( rla_x( il_size), rla_f(il_size) )
		CALL get_rng_stat(td_rnge,  x = rla_x, y = rla_f )
		CALL write_vector_for_plot( "rt_ext_stat.dat", rla_f, rla_x )
		CALL get_rng_stat( td_rnge, x = rla_x, pdf = rla_f )
		CALL write_vector_for_plot( "rt_ext_pdf.dat", rla_f, rla_x )
	END SUBROUTINE save_extern_nstat

END PROGRAM rand_test