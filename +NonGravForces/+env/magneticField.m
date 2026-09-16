function [B, info] = magneticField(coeff, rItrs, options)
%MAGNETICFIELD IGRF main geomagnetic field at Earth-fixed positions, in tesla, Earth-fixed or rotated to any frame.
%
%   B = NonGravForces.env.magneticField(coeff, rItrs)
%   B = NonGravForces.env.magneticField(coeff, rItrs, 'Rotation', M)          % e.g. ITRS -> body
%   [B, info] = NonGravForces.env.magneticField(coeff, rItrs, 'MaxDegree', 1)  % the dipole alone
%
%   Synthesises the internal field of a spherical-harmonic model -- IGRF-14 as
%   io.igrfCoefficients delivers it -- at K positions:
%
%       V = a sum_n (a/r)^(n+1) sum_m [g_n^m cos(m phi) + h_n^m sin(m phi)] P_n^m(cos theta),  B = -grad V
%
%   with a = 6371.2 km, theta the geocentric colatitude, phi the east longitude
%   and P_n^m the Schmidt semi-normalised associated Legendre functions (no
%   Condon-Shortley phase). THE TIME IS THE COEFFICIENTS': pass the struct
%   io.igrfCoefficients(file, 'Epoch', t) returns, whose g and h already follow
%   the model's own time rule. Over one day the secular variation moves the
%   field by less than 0.05 nT, so one coefficient set per day is exact to far
%   below the model's own error; a longer arc takes one set per day.
%
%   THE POLES ARE ORDINARY POINTS. The east component contains m P_n^m / sin(theta);
%   the functions are evaluated as sin^m(theta) T_n^m(cos theta) with T a
%   polynomial from its three-term recursion, so that quotient is
%   sin^(m-1)(theta) T_n^m and nothing is divided by zero.
%
%   UNITS ARE SI HERE. The coefficients are published in nanotesla and stay so in
%   io.igrfCoefficients; B is returned in tesla, the conversion made once, here.
%
% INPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   coeff   1x1   struct  nT     io.igrfCoefficients result evaluated at an epoch (.g, .h
%                                (nMax+1)x(nMax+1), indexed [n+1, m+1]; .nMax; .year)
%   rItrs   3xK   double  m      Positions, Earth-fixed (ITRS), geocentric, |r| > 0
%
% NAME-VALUE ARGUMENTS
%   Name       Size     Type    Units  Description
%   ----       ----     ----    -----  -----------
%   MaxDegree  1x1      double  -      Truncation degree, 1..coeff.nMax; default coeff.nMax
%   Rotation   3x3|3x3xK double -      Maps Earth-fixed vectors to the frame B is returned in
%                                      (default []: B stays Earth-fixed)
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   B      3xK   double  T      Field vector, Earth-fixed or rotated by 'Rotation'
%   info   1x1   struct  -      .itrs_T (3xK Earth-fixed), .northEastDown_T (3xK, geocentric local
%                               frame: X north, Y east, Z down), .intensity_T (1xK),
%                               .maxDegree, .year (coeff.year), .radius_m (1xK), .provenance
%
% ACCURACY
%   The synthesis is exact to rounding for the coefficients given (degree 13:
%   about 1e-13 relative, checked against an independent standard-library
%   synthesis). The MODEL is the main field only: no crustal field beyond degree
%   13, no ionospheric or magnetospheric field, which at LEO leave a few tens to
%   hundreds of nT (a few 1e-3 of the intensity; the GRACE-C MAG1B day in
%   data/samples measures the difference).
%
% ERRORS
%   NonGravForces:magneticField:inputCoefficients  coeff has no g/h at an epoch (pass 'Epoch' to io.igrfCoefficients)
%   NonGravForces:magneticField:inputShape         rItrs not 3xK, or Rotation not 3x3 or 3x3xK
%   NonGravForces:magneticField:inputRange         |r| = 0, or MaxDegree outside 1..nMax
%   NonGravForces:magneticField:numericNaN         NaN/Inf in rItrs, the coefficients or Rotation
%
% REFERENCES
%   Alken P. et al. (2021) International Geomagnetic Reference Field: the thirteenth generation,
%   Earth Planets Space 73:49, eqs 1-3; Winch D.E. et al. (2005) Geophys. J. Int. 160:487
%   (Schmidt semi-normalisation); Langel R.A. (1987) in Geomagnetism vol. 1, ch. 4 (the synthesis).
%
% EXAMPLE
%   c = NonGravForces.io.igrfCoefficients("data/igrf/igrf14coeffs.txt", 'Epoch', 2020.0);
%   B = NonGravForces.env.magneticField(c, [6861.2e3; 0; 0]);     % over the equator at 490 km
%   norm(B) * 1e9                                                 % nT
%
% See also NonGravForces.io.igrfCoefficients, NonGravForces.io.readMAG1B, NonGravForces.kernel.bodyTorques.
%
% KEYWORDS IGRF geomagnetic field spherical harmonics Schmidt magnetic torque
% Created 2026-08-27 by Claude Fable 5 (Anthropic) as a stub; implemented 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    coeff (1,1) struct
    rItrs double
    options.MaxDegree double = []
    options.Rotation double = []
end
unit = "magneticField";
if ~all(isfield(coeff, ["g" "h" "nMax" "year"])) || isempty(coeff.g) || ~isequal(size(coeff.g), size(coeff.h)) || size(coeff.g, 1) < 2
    NonGravForces.fail(unit, "inputCoefficients", "coeff must be io.igrfCoefficients(file, 'Epoch', t): the coefficients AT an epoch -> the tables alone have no time");
end
if ~(ismatrix(rItrs) && size(rItrs, 1) == 3)
    NonGravForces.fail(unit, "inputShape", "rItrs must be 3xK -> got %s", mat2str(size(rItrs)));
end
K = size(rItrs, 2);
if any(~isfinite(rItrs(:))) || any(~isfinite(coeff.g(:))) || any(~isfinite(coeff.h(:)))
    NonGravForces.fail(unit, "numericNaN", "rItrs or the coefficients contain NaN/Inf -> screen the positions");
end
N = coeff.nMax;
if ~isempty(options.MaxDegree)
    N = options.MaxDegree;
end
if ~(isscalar(N) && N == round(N) && N >= 1 && N <= coeff.nMax && N <= size(coeff.g, 1) - 1)
    NonGravForces.fail(unit, "inputRange", "MaxDegree must be an integer in 1..%d -> got %s", coeff.nMax, mat2str(N));
end
rad = vecnorm(rItrs);
if any(rad == 0)
    NonGravForces.fail(unit, "inputRange", "a position is the geocentre (column %d) -> the field is not defined there", find(rad == 0, 1));
end
a = 6371.2e3;                                               % IGRF reference radius [m]
x = rItrs(3, :) ./ rad;                                     % cos(theta)
s = hypot(rItrs(1, :), rItrs(2, :)) ./ rad;                 % sin(theta) >= 0, exact 0 on the axis
phi = atan2(rItrs(2, :), rItrs(1, :));
g = coeff.g * 1e-9; h = coeff.h * 1e-9;                     % nT -> T, once
% T(n+1, m+1, :) with P_n^m = s^m T_n^m (unnormalised Ferrers, no Condon-Shortley), and dT/dx.
T = zeros(N + 1, N + 1, K); dT = zeros(N + 1, N + 1, K);
x3 = reshape(x, 1, 1, K);
for m = 0:N
    Tmm = prod(1:2:(2 * m - 1));                            % (2m-1)!!, 1 for m = 0
    T(m + 1, m + 1, :) = Tmm;
    if m + 1 <= N
        T(m + 2, m + 1, :) = (2 * m + 1) * x3 * Tmm;
        dT(m + 2, m + 1, :) = (2 * m + 1) * Tmm;
    end
    for n = m + 2:N
        T(n + 1, m + 1, :) = ((2 * n - 1) * x3 .* T(n, m + 1, :) - (n + m - 1) * T(n - 1, m + 1, :)) / (n - m);
        dT(n + 1, m + 1, :) = ((2 * n - 1) * (T(n, m + 1, :) + x3 .* dT(n, m + 1, :)) - (n + m - 1) * dT(n - 1, m + 1, :)) / (n - m);
    end
end
Br = zeros(1, K); Bt = zeros(1, K); Bp = zeros(1, K);
rho = a ./ rad;
for n = 1:N
    rn = rho .^ (n + 2);
    for m = 0:n
        if m == 0
            S = 1;
        else
            S = sqrt(2 * factorial(n - m) / factorial(n + m));
        end
        Tn = reshape(T(n + 1, m + 1, :), 1, K) * S;
        dTn = reshape(dT(n + 1, m + 1, :), 1, K) * S;
        c = cos(m * phi); sn = sin(m * phi);
        G = g(n + 1, m + 1) * c + h(n + 1, m + 1) * sn;     % 1xK
        H = -g(n + 1, m + 1) * sn + h(n + 1, m + 1) * c;    % d/dphi of G, over m
        if m == 0
            P = Tn;                                          % s^0
            dPdt = -s .* dTn;                                % dP/dtheta = -s dT/dx
            PoverS = zeros(1, K);                            % multiplied by m = 0 below
        else
            smm1 = s .^ (m - 1);
            P = smm1 .* s .* Tn;
            dPdt = m * smm1 .* x .* Tn - smm1 .* s .* s .* dTn;   % d/dtheta (s^m T) = m s^(m-1) x T - s^(m+1) dT/dx
            PoverS = smm1 .* Tn;
        end
        Br = Br + (n + 1) * rn .* G .* P;
        Bt = Bt - rn .* G .* dPdt;
        Bp = Bp - rn .* m .* H .* PoverS;
    end
end
% local spherical (r, theta, phi) -> Earth-fixed Cartesian
cp = cos(phi); sp = sin(phi);
Bitrs = [Br .* s .* cp + Bt .* x .* cp - Bp .* sp; ...
         Br .* s .* sp + Bt .* x .* sp + Bp .* cp; ...
         Br .* x - Bt .* s];
B = Bitrs;
if ~isempty(options.Rotation)
    M = options.Rotation;
    if ~(isequal(size(M, 1:2), [3 3]) && (ismatrix(M) || size(M, 3) == K))
        NonGravForces.fail(unit, "inputShape", "Rotation must be 3x3 or 3x3x%d -> got %s", K, mat2str(size(M)));
    end
    if any(~isfinite(M(:)))
        NonGravForces.fail(unit, "numericNaN", "Rotation contains NaN/Inf -> check the attitude and the Earth orientation");
    end
    B = reshape(pagemtimes(M, reshape(Bitrs, 3, 1, K)), 3, K);
end
if nargout > 1
    info = struct('itrs_T', Bitrs, 'northEastDown_T', [-Bt; Bp; -Br], 'intensity_T', vecnorm(Bitrs), ...
        'maxDegree', N, 'year', coeff.year, 'radius_m', rad);
    info.provenance = NonGravForces.provenance("env.magneticField", 'Options', struct('MaxDegree', N, 'year', coeff.year, 'rotated', ~isempty(options.Rotation)));
end
end
