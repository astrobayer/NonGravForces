function c = const()
%CONST Constants owned by NonGravForces (everything geodetic lives in referenceSystem).
%
%   c = NonGravForces.const
%
%   One authority per fact (rule R1). This module holds ONLY the constants
%   that no dependency publishes: radiative and gas-kinetic quantities and the
%   mission-specific handbook values. Earth GM, rotation rate, speed of light
%   and ellipsoid parameters are read from referenceSystem.iersConstants /
%   ellipsoidParameters at the call site and are NOT repeated here.
%
% OUTPUT ARGUMENTS
%   Name              Size  Type    Units     Description
%   ----              ----  ----    -----     -----------
%   c.S0              1x1   double  W/m^2     Total solar irradiance at 1 au,
%                                             1361 (IAU 2015 B3; SORCE/TIM,
%                                             TSIS-1). Legacy 1367 is an
%                                             Options field, not a constant.
%   c.au              1x1   double  m         Astronomical unit, IAU 2012
%                                             (149597870700). referenceSystem
%                                             publishes none (verified 2026-08).
%   c.Rsun            1x1   double  m         Nominal solar radius, IAU 2015 B3.
%   c.Rmoon           1x1   double  m         Mean lunar radius, IAU WGCCRE 2015
%                                             (1737.4 km); lunar-eclipse geometry.
%   c.sigmaSB         1x1   double  W/m^2/K^4 Stefan-Boltzmann, CODATA 2018 exact.
%   c.kB              1x1   double  J/K       Boltzmann, exact.
%   c.NA              1x1   double  1/mol     Avogadro, exact.
%   c.Rgas            1x1   double  J/mol/K   Molar gas constant, exact.
%   c.amu             1x1   double  kg        Atomic mass unit, CODATA 2018.
%   c.eV              1x1   double  J         Electron volt, exact (SI 2019).
%   c.torr            1x1   double  Pa        Torr, exact (101325/760); SESAM's Langmuir
%                                             parameters are published in torr^-1.
%   c.speciesName     1x8   string  -         He O N2 O2 Ar H N O_anomalous
%   c.speciesMass     1x8   double  kg/mol    Molar masses of the MSIS species
%   c.lambertDiffuse  1x1   double  -         2/3, Lambertian recoil factor
%   c.graceFO         1x1   struct  -         Handbook (JPL D-56935 Tab.4):
%                                             .massC .massD [kg]
%   c.grace           1x1   struct  -         Launch mass (Bettadpur 2012, GRACE 327-720
%                                             Sec. I.5.3.3): .massA .massB [kg]
%
% ERRORS
%   (none) The function takes no input and cannot fail.
%
% ACCURACY
%   All values are defining or nominal constants; none carries a model error.
%   S0: the observed TSI varies by 0.1 % over a solar cycle; the toolbox
%   therefore exposes Options.srp.tsi to override S0 per epoch if measured
%   TSI is supplied.
%
% REFERENCES
%   IAU 2015 Resolution B3; CODATA 2018; Wen et al. 2019 JPL D-56935 Table 4;
%   Bettadpur 2012 GRACE 327-720 Rev 4.6 Sec. I.5.3.3.
%
% EXAMPLE
%   c = NonGravForces.const; P0 = c.S0 / 299792458
%
% See also referenceSystem.iersConstants, referenceSystem.ellipsoidParameters.
%
% KEYWORDS constants solar irradiance species
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
c.S0      = 1361.0;
c.au      = 149597870700.0;
c.Rsun    = 6.957e8;
c.Rmoon   = 1737.4e3;
c.sigmaSB = 5.670374419e-8;
c.kB      = 1.380649e-23;
c.NA      = 6.02214076e23;
c.Rgas    = 8.31446261815324;
c.amu     = 1.66053906660e-27;
c.eV      = 1.602176634e-19;         % exact since the 2019 SI redefinition of the coulomb
c.torr    = 101325 / 760;            % exact by definition of the standard atmosphere
c.speciesName = ["He","O","N2","O2","Ar","H","N","O_anomalous"];
c.speciesMass = [4.002602, 15.9994, 28.0134, 31.9988, 39.948, 1.00794, 14.0067, 15.9994] * 1e-3;
c.lambertDiffuse = 2/3;
c.graceFO = struct('massC', 601.214, 'massD', 601.209);
c.grace   = struct('massA', 487.2, 'massB', 487.2);       % launch mass; classic MAS1B holds the mission value
end
