!> @file compiler_dependant.F90
!!
!<
!the file extension is capital to force preprocessing
!> @brief defines compiler specifics processing
!!
!<
module compiler_dependant
#if defined(__GFORTRAN__)
!gfortran specifics
   use iso_fortran_env, only: int64
#elif defined(__INTEL_COMPILER)
!ifort specifics
   use iso_fortran_env, only: int64
#elif defined(__IBMC__)
!ibm specifics
   use iso_fortran_env, only: int64
#elif defined(__SUNPRO_F90)
!sunpro specifics
   use iso_fortran_env, only: int64
#else
!this option is especially used for pgf90 to provide a function to test for NaN values
use ieee_arithmetic, only : isnan => ieee_is_nan
#endif
implicit none

#if defined(__GFORTRAN__)
!gfortran specifics
#elif defined(__INTEL_COMPILER)
!ifort specifics
#elif defined(__IBMC__)
!ibm specifics
#elif defined(__SUNPRO_F90)
!sunpro specifics
#else
!this option is especially used with pgf90 12.3-0 (old compilers) for NRL NCOM
   integer, parameter :: int64 = 8
#endif
contains

#if defined(__GFORTRAN__)
!gfortran specifics
#elif defined(__INTEL_COMPILER)
!ifort specifics
#elif defined(__IBMC__)
!ibm specifics
#elif defined(__SUNPRO_F90)
!sunpro specifics
#else
!this option is especially used for pgf90 to provide a getpid() function
   !> @brief Returns the process ID of the current process
   !! @todo write the actual code, for now returns a fixed value
   !<
   function getpid()result(pid)
      integer pid
      pid = 53 !just a prime number, no special meaning
   end function getpid
#endif

end module compiler_dependant