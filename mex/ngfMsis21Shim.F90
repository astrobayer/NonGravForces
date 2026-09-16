! NGFMSIS21SHIM  C-callable entry points into the UNMODIFIED NRLMSIS 2.1 Fortran.
!
! Built by NonGravForces.buildMex together with the NRL sources fetched by
! setupNonGravForces into mex/nrlmsis21 and linked into
! +NonGravForces/+env/private/msis21Gateway, visible to NonGravForces.env.density
! and to nothing else.
!
! THIS FILE CONTAINS NO MODEL. It is this toolbox's own code: two bind(C)
! wrappers that call MSISINIT and MSISCALC exactly as documented in their headers
! (msis_init.F90, msis_calc.F90). No NRL file is modified. This software
! incorporates the MSIS(R) empirical atmospheric model software designed and
! provided by NRL; its use is governed by the Open Source Academic Research License
! Agreement in nrlmsis2.1_license.txt (research, academic and non-profit use only).
!
! The model sources are compiled with -DDBLE, so real(kind=rp) is double
! precision and the conversions below are exact.
!
! TWO THINGS IN THE NRL CODE THIS FILE HAS TO WORK AROUND, both without editing it:
!   * MSISINIT keeps the parameter path and the file name in two character(len=128)
!     variables and opens trim(path)//trim(file). A toolbox installed below a long
!     folder (the session scratch directory measured 146 characters) would be cut
!     silently. The full file name is therefore split across the two arguments,
!     never at a blank (trim would eat it), which admits 256 characters.
!   * LOADPARMSET executes STOP when the file is absent, which ends the MATLAB
!     process that loaded the MEX. The existence of exactly the name it will open
!     is checked here first and reported as a status instead.
!
!   ngf_msis21_init(dir, n, aphistory, anomalous, status)
!     dir        char(n)  directory holding msis21.parm, WITH its trailing separator
!     aphistory  int      nonzero: legacy switch 9 = -1 (the 7-element ap history is read);
!                         zero: switch 9 = 1 (daily Ap only)
!     anomalous  int      nonzero: the total mass density includes anomalous oxygen (the
!                         model's default); zero: lmass_include(9) = .false.
!     status     int      0 success, 2 the full file name exceeds 256 characters or cannot
!                         be split, 3 msis21.parm is not there
!   ngf_msis21_calc(day, utsec, z, lat, lon, f107a, f107, ap, tn, dn, tex)
!     as MSISCALC: day of year, UT seconds, geodetic altitude km, geodetic latitude and
!     longitude deg, 81-day centred and previous-day F10.7, ap(7); tn K, dn(10) SI in
!     MSIS 2.1 order (rho kg/m3, N2 O2 O He H Ar N anomalous-O NO in m^-3), tex K.
!
! Created 2026-09-15 by Claude Opus 5 (Anthropic).
module ngf_msis21_shim
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_char
  use msis_constants, only : rp
  use msis_init, only : msisinit
  use msis_calc, only : msiscalc
  implicit none
contains

  subroutine ngf_msis21_init(dir, n, aphistory, anomalous, status) bind(C, name='ngf_msis21_init')
    integer(c_int), value, intent(in)  :: n
    character(kind=c_char), intent(in) :: dir(*)
    integer(c_int), value, intent(in)  :: aphistory, anomalous
    integer(c_int), intent(out)        :: status
    character(len=*), parameter        :: fname = 'msis21.parm'
    character(len=256)                 :: full
    real(4)                            :: sw(25)
    logical                            :: mass(10), here
    integer                            :: i, m, k
    m = n + len(fname)
    if (n < 1 .or. m > 256) then
      status = 2
      return
    end if
    full = ' '
    do i = 1, n
      full(i:i) = dir(i)
    end do
    full(n+1:m) = fname
    inquire(file=full(1:m), exist=here)
    if (.not. here) then
      status = 3
      return
    end if
    ! The first part keeps as many characters as fit, and never ends in a blank.
    k = min(m, 128)
    do while (k > 0)
      if (full(k:k) /= ' ') exit
      k = k - 1
    end do
    if (k < 1 .or. m - k > 128) then
      status = 2
      return
    end if
    sw = 1.0
    if (aphistory /= 0) sw(9) = -1.0
    mass = .true.
    if (anomalous == 0) mass(9) = .false.
    call msisinit(parmpath=full(1:k), parmfile=full(k+1:m), switch_legacy=sw, lmass_include=mass)
    status = 0
  end subroutine ngf_msis21_init

  subroutine ngf_msis21_calc(day, utsec, z, lat, lon, f107a, f107, ap, tn, dn, tex) &
      bind(C, name='ngf_msis21_calc')
    real(c_double), value, intent(in) :: day, utsec, z, lat, lon, f107a, f107
    real(c_double), intent(in)        :: ap(7)
    real(c_double), intent(out)       :: tn, dn(10), tex
    real(kind=rp)                     :: xap(7), xdn(10), xtn, xtex
    xap = real(ap, rp)
    call msiscalc(real(day, rp), real(utsec, rp), real(z, rp), real(lat, rp), real(lon, rp), &
                  real(f107a, rp), real(f107, rp), xap, xtn, xdn, tex=xtex)
    tn = real(xtn, c_double)
    dn = real(xdn, c_double)
    tex = real(xtex, c_double)
  end subroutine ngf_msis21_calc

end module ngf_msis21_shim
