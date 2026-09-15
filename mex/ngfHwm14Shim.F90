! NGFHWM14SHIM  A C-callable entry point into the UNMODIFIED HWM14 Fortran (NRL).
!
! Built by NonGravForces.buildMex together with hwm14.f90 as fetched by
! setupNonGravForces into mex/hwm14, and linked into
! +NonGravForces/+env/private/hwm14Gateway, visible to NonGravForces.env.wind and
! to nothing else.
!
! THIS FILE CONTAINS NO MODEL. It is this toolbox's own code: one bind(C) wrapper
! that calls the documented entry HWM14(IYD, SEC, ALT, GLAT, GLON, STL, F107A,
! F107, AP, W) exactly as its header describes it. No NRL file is modified.
!
! THE MODEL IS SINGLE PRECISION. Its arguments are real(4) and so is its
! arithmetic; the conversions below are the model's own resolution, which NRL's
! readme quotes as +-1 mm/s between compilers. STL, F107A, F107 and AP(1) are
! declared "not used" by the model and are passed as zero.
!
!   ngf_hwm14_calc(iyd, sec, alt, glat, glon, ap2, w)
!     iyd   int      year and day of year as yyddd
!     sec   double   UT seconds of that day
!     alt   double   geodetic altitude, km
!     glat  double   geodetic latitude, deg
!     glon  double   geodetic longitude, deg
!     ap2   double   current 3-hour ap; negative: quiet-time wind only (no DWM07)
!     w     double(2) out: meridional (+north) and zonal (+east) wind, m/s
!
! Created 2026-09-15 by Claude Opus 5 (Anthropic).
module ngf_hwm14_shim
  use, intrinsic :: iso_c_binding, only : c_int, c_double
  implicit none
  interface
    subroutine hwm14(iyd, sec, alt, glat, glon, stl, f107a, f107, ap, w)
      integer(4), intent(in) :: iyd
      real(4), intent(in)    :: sec, alt, glat, glon, stl, f107a, f107
      real(4), intent(in)    :: ap(2)
      real(4), intent(out)   :: w(2)
    end subroutine hwm14
  end interface
contains

  subroutine ngf_hwm14_calc(iyd, sec, alt, glat, glon, ap2, w) bind(C, name='ngf_hwm14_calc')
    integer(c_int), value, intent(in) :: iyd
    real(c_double), value, intent(in) :: sec, alt, glat, glon, ap2
    real(c_double), intent(out)       :: w(2)
    real(4)                           :: ap(2), w4(2)
    ap(1) = 0.0
    ap(2) = real(ap2, 4)
    call hwm14(int(iyd, 4), real(sec, 4), real(alt, 4), real(glat, 4), real(glon, 4), &
               0.0, 0.0, 0.0, ap, w4)
    w = real(w4, c_double)
  end subroutine ngf_hwm14_calc

end module ngf_hwm14_shim
