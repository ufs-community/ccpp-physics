!> This module provides the Fortran KIND parameters for REAL variables in the SAMF scheme.
module mo_conv_kind
  use, intrinsic :: iso_c_binding, only: c_float, c_double
  implicit none
  public

  ! Define standard single and double precision kinds
  integer, parameter :: dp = c_double, sp = c_float

  ! Floating point working precision
  ! This is controlled by the -DSAMFDEEP_USE_SP macro
#ifdef SAMFDEEP_USE_SP
  integer, parameter :: conv_wp = sp
#else
  integer, parameter :: conv_wp = dp
#endif

end module mo_conv_kind
