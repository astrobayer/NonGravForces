! NGFDTM2020SHIM  C-callable entry points into CNES's UNMODIFIED DTM2020 (operational: F10.7, Kp).
!
! Built by NonGravForces.buildMex together with the CNES subroutine file fetched by
! setupNonGravForces into mex/dtm2020 (dtm2020_F107_Kp-subr_MCM.f90 of SWAMI's MCM
! repository, pinned by commit and SHA-256) and linked into
! +NonGravForces/+env/private/dtm2020Gateway, visible to NonGravForces.env.density
! and to nothing else.
!
! THIS FILE CONTAINS NO MODEL. It is this toolbox's own code: two bind(C) wrappers
! that call LECDTM and DTM3 exactly as their headers and SWAMI's documentation
! (docs/dtm2020.rst) describe. No CNES file is modified. The DTM software is
! property of CNES; its licence (mex/dtm2020/LICENSE_MCM_DTM.txt) grants academic,
! non-commercial use only and forbids modification and redistribution.
!
! PRECISION. CNES's benchmark program says "compile [double precision]", and the
! target compiles the CNES file with -fdefault-real-8. The locals below are DEFAULT
! real on purpose, so that they have exactly the kind DTM3 was compiled with whatever
! the flag; the conversions from and to c_double are then exact in the double build.
!
! WHAT THIS FILE HAS TO WORK AROUND, without editing the CNES code:
!   * LECDTM reads its unit with no IOSTAT, so a missing file would end the MATLAB
!     process that loaded the MEX. The file is opened here with IOSTAT and a failure is
!     reported as a status instead.
!
!   ngf_dtm2020_init(path, n, status)
!     path    char(n)  full name of DTM_2020_F107_Kp.dat
!     status  int      0 success, 3 the file cannot be opened
!   ngf_dtm2020_calc(day, f, fbar, kp3, kp24, alt, hl, lat, lon, tz, tinf, ro, d, wmm)
!     as DTM3: day of year, F10.7 at t - 24 h, mean F10.7 of the last 81 days, Kp
!     delayed by 3 h, mean Kp of the last 24 h, altitude km (> 120), local time rad,
!     latitude and longitude rad; tz and tinf K, ro and d(6) (H He O N2 O2 N) g/cm^3,
!     wmm the mean molecular mass in g.
!
! Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
module ngf_dtm2020_shim
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char
  implicit none
  external :: lecdtm, dtm3
contains

  subroutine ngf_dtm2020_init(path, n, status) bind(C, name='ngf_dtm2020_init')
    integer(c_int), value, intent(in)  :: n
    character(kind=c_char), intent(in) :: path(*)
    integer(c_int), intent(out)        :: status
    character(len=n)                   :: full
    integer                            :: i, u, ios
    do i = 1, n
      full(i:i) = path(i)
    end do
    open(newunit=u, file=full, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      status = 3
      return
    end if
    call lecdtm(u)
    close(u)
    status = 0
  end subroutine ngf_dtm2020_init

  subroutine ngf_dtm2020_calc(day, f107, f107bar, kp3, kp24, alt, hl, lat, lon, tz, tinf, ro, d, wmm) &
      bind(C, name='ngf_dtm2020_calc')
    real(c_double), value, intent(in) :: day, f107, f107bar, kp3, kp24, alt, hl, lat, lon
    real(c_double), intent(out)       :: tz, tinf, ro, d(6), wmm
    real :: xf(2), xfbar(2), xkp(4), xd(6), xtz, xtinf, xro, xwmm
    xf = [real(f107), 0.]
    xfbar = [real(f107bar), 0.]
    xkp = [real(kp3), 0., real(kp24), 0.]
    call dtm3(real(day), xf, xfbar, xkp, real(alt), real(hl), real(lat), real(lon), &
              xtz, xtinf, xro, xd, xwmm)
    tz = real(xtz, c_double)
    tinf = real(xtinf, c_double)
    ro = real(xro, c_double)
    d = real(xd, c_double)
    wmm = real(xwmm, c_double)
  end subroutine ngf_dtm2020_calc

end module ngf_dtm2020_shim
