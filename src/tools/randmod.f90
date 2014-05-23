!> \file randmod.f90
!! Innocent defined module for normal random numbers
!<

!> Module for random numbers generation
!! used for Data Assimilation purpose
!! It is mandatory to call init_normal_rand to initialize the module before using it
!! One may call finilize_normal_rand to finilize the module
!! The method used here consists in a quasi inversion of the CDF. The space is subdivided into small intervals. Given a uniform random value between 0 and 1, the CDF of the normal distribution is inverted at the value. When the inversion falls into an interval, the middle of the interval is considered as the associated normal random value.
!<
MODULE randmod
  USE general_constant
  USE general_tools, only: search_pos_inc, init_random_seed
  USE debug_tools
IMPLICIT NONE
  !> default number of possible values at the right of the mean
  INTEGER, PARAMETER :: DEFAULT_NRIGHT = 4096!1024

  !> user defined type to keep and manage information on a random numbers generator
  TYPE rng_environment
    LOGICAL :: l_init = .FALSE. !< Says if the radom number generator environment is initialized or not
    LOGICAL :: l_stat = .FALSE. !< Save statistics (frequency) for analysis? Statistics can be used to check the Gaussianity of a gaussian generator
    !for statistics
    INTEGER :: i_nCall  = 0     !< number of call to the generator, used for statistics
    INTEGER :: i_lb, i_ub       !< bounds of the statistics array
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: ra_rands!<Vector of possible values
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: ra_cdf  !< CDF at possibles values

    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: ra_pdf  !< PDF of the generator
    INTEGER, DIMENSION(:), ALLOCATABLE :: ia_stats !< frequency of random numbers generated in a set of predefined intervals
  END TYPE rng_environment

  !> user type to store parameter of random Gaussian number generator
  TYPE, PRIVATE :: gaussian_param
    INTEGER :: i_nright = 0     !< Number of Number of possible values (greater than the mean) that can be generated
    REAL(KIND=cp) :: r_radius   !< Distance from the mean to the maximum/minimum number that can be generated
    REAL(KIND=cp) :: r_mean     !< mean of the distribution
    REAL(KIND=cp) :: r_sigma    !< Standard deviation of the distribution
    REAL(KIND=cp) :: r_step     !< step between possible successive values
    TYPE(rng_environment) :: t_env  !< environment of the random number generator
  END TYPE gaussian_param

  INTERFACE normal_random
    MODULE PROCEDURE normal_random_scalar
    MODULE PROCEDURE normal_random_vector
  END INTERFACE normal_random

  TYPE(gaussian_param), PRIVATE, SAVE :: tm_gaussp !< parameter for the Gaussian generator

CONTAINS

!--- info
  !> \brief Initializes the environment for a random number generator
  !! \param[in,out] td_rnge environment to be initialized
  !! \param[in] id_lb lower bound index of the environments arrays
  !! \param[in] id_ub upper bound index of the environments arrays
  !! \param[in] stat status of the statistics
  !<
  SUBROUTINE init_rng_environment(td_rnge, id_lb, id_ub, stat)
    TYPE(rng_environment), INTENT(INOUT) :: td_rnge
    INTEGER, INTENT(IN) :: id_lb, id_ub
    LOGICAL, INTENT(IN), OPTIONAL :: stat

    IF ( PRESENT( stat) ) THEN
      td_rnge%l_stat = stat
    ELSE
      td_rnge%l_stat = .FALSE.
    END IF
    td_rnge%i_lb = id_lb
    td_rnge%i_ub = id_ub
    IF( ALLOCATED(td_rnge%ra_rands) )DEALLOCATE(td_rnge%ra_rands)
    IF( ALLOCATED(td_rnge%ra_cdf  ) )DEALLOCATE(td_rnge%ra_cdf  )
    ALLOCATE( td_rnge%ra_rands(td_rnge%i_lb:td_rnge%i_ub),&
              td_rnge%ra_cdf  (td_rnge%i_lb:td_rnge%i_ub) &
    )
    td_rnge%ra_rands = 0.0_cp
    td_rnge%ra_cdf   = 0.0_cp
    !stats
    IF(td_rnge%l_stat)THEN
      IF( ALLOCATED(td_rnge%ra_pdf) )DEALLOCATE(td_rnge%ra_pdf)
      IF( ALLOCATED(td_rnge%ia_stats) )DEALLOCATE(td_rnge%ia_stats)
      ALLOCATE(td_rnge%ra_pdf(td_rnge%i_lb:td_rnge%i_ub),&
              td_rnge%ia_stats(td_rnge%i_lb:td_rnge%i_ub) &
      )
      td_rnge%i_nCall  = 0
      td_rnge%ia_stats = 0
      td_rnge%ra_pdf   = 0.0_cp
    END IF
    td_rnge%l_init = .TRUE.
  END SUBROUTINE init_rng_environment

  !> \brief finalizes the environment for a random number generator
  !! \param[in,out] td_rnge environment to be initialized
  !<
  SUBROUTINE finilize_rng_environment(td_rnge)
    TYPE(rng_environment), INTENT(INOUT) :: td_rnge

    !statistics
    IF(td_rnge%l_init)THEN
      IF(td_rnge%l_stat)THEN
        DEALLOCATE(td_rnge%ra_pdf, td_rnge%ia_stats )
      END IF
      DEALLOCATE( td_rnge%ra_rands, td_rnge%ra_cdf )
      td_rnge%i_lb =  0
      td_rnge%i_ub = -1! an upper bound less than the lower bound
      td_rnge%i_nCall  = -1
      td_rnge%l_init = .FALSE.
    END IF
  END SUBROUTINE finilize_rng_environment

  !> \brief generates a random number from the distribution associated with the random number generator environment
  !! \param[in,out] td_rnge environment of the random number generator to be used
  !<
  FUNCTION generate_rand(td_rnge) RESULT(rl_rand)
    TYPE(rng_environment), INTENT(INOUT) :: td_rnge
    REAL(KIND=cp) :: rl_cdf, rl_rand
    INTEGER :: il_pos, il_ntry

    IF(.NOT.td_rnge%l_init)THEN
      CALL debug("In generate_rand: random number generator is not initialized", tag=dALLWAYS)
      CALL ABORT()
    END IF

    CALL RANDOM_NUMBER(rl_cdf)
    il_pos = search_pos_inc(rl_cdf, td_rnge%ra_cdf, il_ntry) + td_rnge%i_lb - 1
    !WRITE(*,FMT="(A,ES13.5E2,A,I5,A,I5,A)")"Find ", rl_cdf, " At ", il_pos, " with ", il_ntry, " tries"
    rl_rand = td_rnge%ra_rands(il_pos)
    IF(td_rnge%l_stat)THEN
      td_rnge%i_nCall = td_rnge%i_nCall + 1
      td_rnge%ia_stats(il_pos) = td_rnge%ia_stats(il_pos) + 1!for stats
    END IF
  END FUNCTION generate_rand

!   !> \brief updates the statistics on a random number generator using the cfd
!   !! \param[in,out] td_rnge structure that keeps stats
!   !! \param[in] index of the position to be updated
!   !! this function is used to update statistics for external random number generator
!   !<
!   SUBROUTINE update_rnge_stat(td_rnge, pos)
!     TYPE(rng_environment), INTENT(INOUT) :: td_rnge
!     INTEGER, INTENT(IN) :: pos
!     !local variables
!     INTEGER :: il_ntry
!
!     il_pos = search_pos_inc(rd_rand, td_rnge%ra_rands, il_ntry) + td_rnge%i_lb - 1
!     IF(td_rnge%l_stat)THEN
!       td_rnge%i_nCall = td_rnge%i_nCall + 1
!       td_rnge%ia_stats(il_pos) = td_rnge%ia_stats(il_pos) + 1!for stats
!     END IF
!   END SUBROUTINE update_rnge_stat

  !> returns the statistics of the random number generator
  !< \param[in] td_rnge, environment of the random number generator to be used
  !! \param[out] ssize, size of the statistic vector
  !! \param[out] x vector of abcissa of the statistic points
  !! \param[out] y vector of estimated pdf based on calls to ramdom_normal since the last call to init_normal_rand
  !! \param[out] pdf vector of true pdf
  !! \param[out] stat call stat, return 0 if OK, -1 if statistics are not saved (call to init_normal_rand with stat = .FALSE., default)
  !! \detail this routine return usefull information to check the gaussianity of the seauence of number generated since the last call to init_normal_rand
  !! \this routine should be called only if init_normal_rand has been called with stat = .TRUE.
  !<
  SUBROUTINE get_rng_stat(td_rnge, ssize, x, y, pdf, stat)
    TYPE(rng_environment), INTENT(IN) :: td_rnge
    INTEGER, INTENT(OUT), OPTIONAL :: ssize, stat
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT), OPTIONAL :: x, y, pdf
    !local var
    INTEGER :: il_tmp

    IF( PRESENT(ssize) ) ssize = SIZE(td_rnge%ia_stats)
    IF(td_rnge%l_stat)THEN
      IF( PRESENT(x) ) x = td_rnge%ra_rands
      IF( PRESENT(y) )THEN
        !the count of values in a small interval is assumed to be the aproximation (rectangle methods) of the probability of a random values to be in this interval,
        il_tmp = td_rnge%i_lb
        y = ( REAL(td_rnge%ia_stats, cp)/REAL(td_rnge%i_nCall,cp) )/( td_rnge%ra_rands(il_tmp+1)-td_rnge%ra_rands(il_tmp) )
      END IF
      IF( PRESENT(pdf) ) pdf = td_rnge%ra_pdf
      IF( PRESENT(stat) ) stat = 0
    ELSE
      IF( PRESENT(stat) ) stat = -1
    END IF
  END SUBROUTINE get_rng_stat

!------------------
  !> \brief Generates a gaussian random number of preset standard deviation and mean
  !! initially named normal_rand_scalar_func and defined in the interface normal_rand
  !! Since the other porcedure of the interface normal_rand are subroutines, it did not respect the standard, although working with intel compilers
  !<
  FUNCTION normal_rand() RESULT(rl_rand)
    !local var
    REAL(KIND=cp) :: rl_rand
    !REAL(KIND=cp) :: rl_sign, rl_cdf
    !INTEGER :: il_ntry, il_pos

    rl_rand = generate_rand(tm_gaussp%t_env)
    !il_pos = search_pos_inc(rl_cdf, tm_gaussp%t_env%ra_cdf, il_ntry) + tm_gaussp%t_env%i_lb - 1
    !rl_rand = tm_gaussp%t_env%ra_rands(il_pos)
  END FUNCTION normal_rand

  !> \brief Generates a gaussian random number of preset standard deviation and mean
  !!
  !<
  SUBROUTINE normal_random_scalar(rd_rand)
    REAL(KIND=cp), INTENT(OUT) :: rd_rand
    !REAL(KIND=cp) :: rl_sign, rl_cdf
    !INTEGER :: il_ntry, il_pos

    rd_rand = generate_rand(tm_gaussp%t_env)
    !il_pos = search_pos_inc(rl_cdf, tm_gaussp%t_env%ra_cdf, il_ntry) + tm_gaussp%t_env%i_lb - 1
    !rd_rand = tm_gaussp%t_env%ra_rands(il_pos)
  END SUBROUTINE normal_random_scalar

  !> \brief Generates a vector of independant gaussian random numbers of preset standard deviation and mean
  !!
  !<
  SUBROUTINE normal_random_vector(rda_rand)
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_rand
    !REAL(KIND=cp) :: rl_sign, rl_cdf
    INTEGER :: ibi!, il_ntry, il_pos

    DO ibi = 1, SIZE(rda_rand)
      CALL init_random_seed()
      rda_rand(ibi) = generate_rand(tm_gaussp%t_env)
      !il_pos = search_pos_inc(rl_cdf, tm_gaussp%t_env%ra_cdf, il_ntry) + tm_gaussp%t_env%i_lb - 1
      !rd_rand = tm_gaussp%t_env%ra_rands(il_pos)
    END DO
  END SUBROUTINE normal_random_vector

  !> \brief Initializes the Gaussian generator
  !! \param[in] mu mean of the distribution
  !! \param[in] sigma Standard deviation of the distribution
  !! \param[in] radius Distance from the mean to the maximum/minimum number that can be generated, default is set to 4*rd_sigma
  !! \param[in] nright Number of possible values greater than the mean, default is set to 127
  !! \param[in] stat says if it is necessary to save statistics
  !! \details the total number of possible values that can be generated is 2*id_nright + 1
  !<
  SUBROUTINE init_normal_rand(mu, sigma, radius, nright, stat)
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: mu, sigma
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: radius
    INTEGER, OPTIONAL, INTENT(IN) :: nright
    LOGICAL, OPTIONAL, INTENT(IN) :: stat
    REAL(KIND=cp) :: rl_hstep
    INTEGER :: ibi, il_lb, il_ub

    IF ( PRESENT(mu) ) THEN
      tm_gaussp%r_mean  = mu
    ELSE
      tm_gaussp%r_mean = 0.0_cp
    END IF
    IF ( PRESENT(sigma) ) THEN
      tm_gaussp%r_sigma = sigma
    ELSE
      tm_gaussp%r_sigma = 1.0_cp
    END IF
    IF( PRESENT(radius) ) THEN
      tm_gaussp%r_radius = radius
    ELSE
      tm_gaussp%r_radius = 4.0_cp*tm_gaussp%r_sigma
    END IF
    IF( PRESENT(nright) ) THEN
      tm_gaussp%i_nright = nright
    ELSE
      tm_gaussp%i_nright = DEFAULT_NRIGHT
    END IF
    il_lb = -tm_gaussp%i_nright
    il_ub =  tm_gaussp%i_nright-1

    CALL init_rng_environment(tm_gaussp%t_env, il_lb, il_ub, stat)

    tm_gaussp%r_step = tm_gaussp%r_radius/( REAL(tm_gaussp%i_nright, cp)-0.5_cp)
    rl_hstep= tm_gaussp%r_step/2.0_cp
    DO ibi = il_lb, il_ub
      tm_gaussp%t_env%ra_rands(ibi) = (ibi)*tm_gaussp%r_step + tm_gaussp%r_mean
      tm_gaussp%t_env%ra_cdf(ibi)   = normal_cdf( tm_gaussp%t_env%ra_rands(ibi)-rl_hstep, tm_gaussp%r_mean, tm_gaussp%r_sigma)
      !stats
      IF(tm_gaussp%t_env%l_stat)THEN
        tm_gaussp%t_env%ia_stats(ibi) = 0
        tm_gaussp%t_env%ra_pdf(ibi) = &
          normal_pdf( tm_gaussp%t_env%ra_rands(ibi), tm_gaussp%r_mean, tm_gaussp%r_sigma )
      END IF
      !CALL debug(ibi, 'init_normal_rand loop..', tag=dALLWAYS)
    END DO
    CALL init_random_seed()
  END SUBROUTINE init_normal_rand

  SUBROUTINE copy_normal_rnge(td_rnge)
    TYPE(rng_environment), INTENT(INOUT) :: td_rnge

    CALL finilize_rng_environment(td_rnge)
    td_rnge%l_init = tm_gaussp%t_env%l_init
    IF(tm_gaussp%t_env%l_init)THEN
      CALL init_rng_environment(td_rnge, tm_gaussp%t_env%i_lb, tm_gaussp%t_env%i_ub, tm_gaussp%t_env%l_stat)
      td_rnge%ra_rands = tm_gaussp%t_env%ra_rands
      td_rnge%ra_cdf   = tm_gaussp%t_env%ra_cdf
      IF (td_rnge%l_stat) THEN
        td_rnge%ra_pdf   = tm_gaussp%t_env%ra_pdf
        td_rnge%ia_stats = tm_gaussp%t_env%ia_stats
        td_rnge%i_nCall  = tm_gaussp%t_env%i_nCall
      END IF
    END IF
  END SUBROUTINE copy_normal_rnge

  !> returns statistics of the random module
  !! \param[out] ssize size of the statistic vector
  !! \param[out] x vector of abcissa of the statistic points
  !! \param[out] y vector of estimated pdf based on calls to ramdom_normal since the last call to init_normal_rand
  !! \param[out] pdf vector of true pdf
  !! \param[out] stat call stat, return 0 if OK, -1 if statistics are not saved (call to init_normal_rand with stat = .FALSE., default)
  !! \details this routine return usefull information to check the gaussianity of the seauence of number generated since the last call to \a init_normal_rand.
  !! This routine should be called only if \a init_normal_rand has been called with stat = .TRUE.
  !<
  SUBROUTINE get_normal_stat(ssize, x, y, pdf, stat)
    INTEGER, INTENT(OUT), OPTIONAL :: ssize, stat
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT), OPTIONAL :: x, y, pdf
    !INTEGER :: ibi

    CALL get_rng_stat(tm_gaussp%t_env, ssize, x, y, pdf, stat)
  END SUBROUTINE get_normal_stat

  SUBROUTINE finilize_normal_rand()
    !TYPE(dyn_rVector), DIMENSION(3) :: tla_save
    IF(tm_gaussp%i_nright>0)THEN
      tm_gaussp%i_nright = 0
      CALL finilize_rng_environment(tm_gaussp%t_env)
    END IF
  END SUBROUTINE finilize_normal_rand

  !> \brief Computes the normal probability density function at point x
  !! \param[in] rd_x point where the PDF is calculated
  !! \param[in] mu (optional, default 0) mean of the distribution
  !! \param[in] sigma (optional, default 1) standard deviation of the distribution
  !<
  FUNCTION normal_pdf(rd_x, mu, sigma) RESULT(rl_pdf)
    REAL(KIND=cp), INTENT(IN) :: rd_x
    REAL(KIND=cp), INTENT(IN), OPTIONAL :: mu, sigma
    REAL(KIND=cp) :: rl_pdf, rl_mu, rl_sigma

    IF ( PRESENT(mu) ) THEN
      rl_mu = mu
    ELSE
      rl_mu = 0.0
    END IF
    IF ( PRESENT(sigma) ) THEN
      rl_sigma = sigma
    ELSE
      rl_sigma = 1.0
    END IF

    rl_pdf = ( 1/( rl_sigma*SQRT(2*rp_pi) ) )*exp(-(rd_x-rl_mu)**2/(2.0*rl_sigma**2) )
  END FUNCTION normal_pdf

  !> \brief Computes the normal Cumulative distribution function at point x
  !! \param[in] rd_x point where the CDF is calculated
  !! \param[in] mu (optional, default 0) mean of the distribution
  !! \param[in] sigma (optional, default 1) standard deviation of the distribution
  !<
  FUNCTION normal_cdf(rd_x, mu, sigma) RESULT(rl_cdf)
    REAL(KIND=cp), INTENT(IN) :: rd_x
    REAL(KIND=cp), INTENT(IN), OPTIONAL :: mu, sigma
    REAL(KIND=cp) :: rl_cdf, rl_mu, rl_sigma

    IF ( PRESENT(mu) ) THEN
      rl_mu = mu
    ELSE
      rl_mu = 0.0
    END IF
    IF ( PRESENT(sigma) ) THEN
      rl_sigma = sigma
    ELSE
      rl_sigma = 1.0
    END IF

    rl_cdf = 0.5*( 1 + erf( (rd_x-rl_mu)/(rl_sigma*SQRT(2.0)) ) )
  END FUNCTION normal_cdf

END MODULE randmod