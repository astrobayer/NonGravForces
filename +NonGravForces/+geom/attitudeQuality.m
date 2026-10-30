function Q = attitudeQuality(epochSca, qSca, epochImu, angleImu, gyroId, options)
%ATTITUDEQUALITY Attitude-quality diagnostic: star-camera rates against the IMU gyros, attitude against the nominal frame, 1 and 2 cpr.
%
%   Q = NonGravForces.geom.attitudeQuality(epochSca, qSca, epochImu, angleImu, gyroId)
%   Q = NonGravForces.geom.attitudeQuality(..., 'Orbit', struct('r', r, 'v', v), 'AlongVelocity', false)
%   Q = NonGravForces.geom.attitudeQuality(..., 'Period_s', 5640)
%
%   WHAT IT COMPARES. The attitude of GRACE-FO has one star-camera product (SCA1B;
%   TU Graz's operational attitude is the same numbers, CLAUDE.md section 7) and
%   one independent sensor of its CHANGE: the IMU gyros (io.readIMU1B), which
%   accumulate the angle turned about their own axes. Both are reduced to the
%   same quantity on the same intervals: over each pair of consecutive SCA1B
%   epochs [t_k, t_k+1],
%     - the star cameras give the body rotation R_k' R_k+1 (R = BODY->GCRS), whose
%       rotation vector over dt is the mean body rate omega_k (skew(omega) = R'Rdot,
%       the convention of sensor.Electrostatic);
%     - each gyro g gives (theta_g(t_k+1) - theta_g(t_k)) / dt, its accumulated
%       angle interpolated linearly to the SCA1B epochs (8 Hz samples, so the
%       interpolation spans at most 0.125 s of a rate that changes at the orbital
%       period).
%   No gyro alignment is assumed: per gyro, rate_g = a_g . omega + b_g is fitted by
%   least squares. What the fit cannot absorb -- the RESIDUAL, computed as the exact
%   orthogonal projection of the gyro rates off the span of [omega, 1] -- is attitude
%   error of either sensor that is not a constant misalignment, scale or rate bias:
%   noise, and the once- and twice-per-revolution signatures of star-camera errors
%   this diagnostic exists to expose (N-12). Its 1 and 2 cpr amplitudes are fitted
%   when the orbital period is known.
%
%   THE RESIDUAL IS THE RESULT; THE FITTED AXIS IS NOT. An Earth-pointing satellite
%   turns almost uniformly about body y at the orbital rate (1.10e-3 rad/s on
%   2020-01-01 C, the x and z rates ~1e-6), so omega is nearly constant and the
%   columns of [omega, 1] are nearly collinear: a gyro's axis, scale and bias are
%   not separable on such an arc (the independent oracle finds axes that move by
%   tens of degrees under a change of convention that changes the residual by 15 %).
%   The parameters are therefore returned as the minimum-norm solution together
%   with the condition number of the column-normalised design, and nothing should
%   be read from them unless that number is small (an attitude manoeuvre, not a
%   nominal arc).
%
%   With 'Orbit' (GCRS r, v at the SCA1B epochs) the attitude is also compared with
%   the nominal Earth-pointing frame (geom.AttitudeProvider.nominal): the rotation
%   vector of R_nominal' R, small angles about body x, y, z (roll, pitch, yaw), and
%   their 1 and 2 cpr amplitudes. The period then comes from the orbit (vis-viva
%   mean semi-major axis) unless 'Period_s' is given.
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   epochSca  1x1   struct  -      referenceSystem epoch struct, K entries (io.readSCA1B), strictly increasing
%   qSca      4xK   double  -      BODY->GCRS quaternions at epochSca
%   epochImu  1x1   struct  -      referenceSystem epoch struct, N entries (io.readIMU1B), any scale
%   angleImu  1xN   double  rad    accumulated gyro angle of each record
%   gyroId    1xN   double  -      gyro of each record
%
% NAME-VALUE ARGUMENTS
%   Name               Size  Type     Units  Description
%   ----               ----  ----     -----  -----------
%   Orbit              1x1   struct   m,m/s  .r .v, 3xK GCRS at epochSca; default [] (no nominal comparison)
%   AlongVelocity      1x1   logical  -      body +X along the velocity for the nominal frame; default false
%                                            (GRACE-C on 2020-01-01; GRACE-D flies true)
%   Period_s           1x1   double   s      orbital period for the cpr fits; default NaN (from Orbit, else none)
%   MaxGap_s           1x1   double   s      longest SCA1B interval and gyro sample spacing used; default 2
%   JumpThreshold_rad  1x1   double   rad    a gyro step larger than this is a reset or wrap: the gyro
%                                            series is split there and no interval spans it; default deg2rad(1)
%
% OUTPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   Q     1x1   struct  -      .t_s (1xM mid-interval times, s from the first SCA1B epoch), .omega_rad_s
%                              (3xM SCA1B mean body rates), .gyros (1xG ids), and per gyro in the struct
%                              array .gyro(g): .id, .rate_rad_s (1xM, NaN where not covered), .axis (3x1
%                              unit, SRF), .scale, .bias_rad_s (the minimum-norm fit; see .condition),
%                              .condition (2-norm condition number of the column-normalised design),
%                              .residual_rad_s (1xM), .rms_rad_s, .rmsRate_rad_s, .n, .cpr (1x2 amplitudes
%                              at 1 and 2 cpr, rad/s; NaN without a period); .period_s;
%                              with Orbit: .nominal.t_s (1xK), .nominal.angles_rad (3xK roll pitch yaw),
%                              .nominal.mean_rad, .nominal.rms_rad, .nominal.cpr (3x2 rad)
%
% ACCURACY
%   The rotation vector of R_k' R_k+1 is exact (atan2 form, no small-angle step); a
%   mean rate over an interval is compared with a mean rate over the same interval,
%   so no rate is interpolated. The gyro angle is interpolated linearly between 8 Hz
%   samples: for a rate changing at the orbital frequency (1.1e-3 rad/s at
%   1.1e-3 rad/s per rad) the error over 0.125 s is below 1e-10 rad. The residual is
%   an orthogonal projection (economy QR), exact to rounding however ill-conditioned
%   the fit; the parameters are minimum-norm (lsqminnorm) and carry the condition
%   number in .condition. The cpr amplitudes are least-squares sinusoids on whatever
%   span is given and are only meaningful over a span of several revolutions.
%
% ERRORS
%   NonGravForces:attitudeQuality:inputShape    sizes disagree, fewer than 2 SCA1B epochs, or an epoch is
%                                               not a referenceSystem epoch struct
%   NonGravForces:attitudeQuality:inputOrbit    Orbit is not .r/.v 3xK at the SCA1B epochs
%   NonGravForces:attitudeQuality:dataOverlap   no SCA1B interval is covered by a gyro, or a gyro has fewer
%                                               than 10 covered intervals (the fit needs them)
%   NonGravForces:AttitudeProvider:*, referenceSystem:*   propagated unwrapped
%
% REFERENCES
%   JPL D-56935 (SCA1B, IMU1B); Klinger and Mayer-Guerr (2016), Adv. Space Res. 58, 1597-1609
%   (star camera and angular accelerometer combined for GRACE attitude); HANDOVER N-12 (#69).
%
% EXAMPLE
%   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%   d = fullfile(root, 'data', 'samples', 'multiday');
%   [eS, qS] = NonGravForces.io.readSCA1B(fullfile(d, 'SCA1B_2020-01-01_C_04.txt'));
%   [eI, aI, auxI] = NonGravForces.io.readIMU1B(fullfile(d, 'IMU1B_2020-01-01_C_04.txt'));
%   Q = NonGravForces.geom.attitudeQuality(eS, qS, eI, aI, auxI.gyroId);
%   [Q.gyro.rms_rad_s]                          % what the fitted axes leave of each gyro
%
% See also NonGravForces.io.readIMU1B, NonGravForces.io.readSCA1B, NonGravForces.geom.AttitudeProvider.
%
% KEYWORDS attitude quality gyro IMU star camera SCA1B rate twice per revolution cpr
% Created 2026-10-30 by Claude (Anthropic) for Matthias Weigelt.
arguments
    epochSca (1,1) struct
    qSca double
    epochImu (1,1) struct
    angleImu double
    gyroId double
    options.Orbit = []
    options.AlongVelocity (1,1) logical = false
    options.Period_s (1,1) double = NaN
    options.MaxGap_s (1,1) double {mustBePositive} = 2
    options.JumpThreshold_rad (1,1) double {mustBePositive} = deg2rad(1)
end
unit = "attitudeQuality";
isEpoch = @(e) all(isfield(e, ["day","sec","scale"]));
if ~isEpoch(epochSca) || ~isEpoch(epochImu)
    NonGravForces.fail(unit, "inputShape", "epochSca and epochImu must be referenceSystem epoch structs -> io.readSCA1B and io.readIMU1B return them");
end
K = numel(epochSca.day);
N = numel(epochImu.day);
if K < 2 || ~isequal(size(qSca), [4 K]) || ~isequal(size(angleImu), [1 N]) || ~isequal(size(gyroId), [1 N])
    NonGravForces.fail(unit, "inputShape", "need qSca 4xK for K >= 2 SCA1B epochs and angleImu, gyroId 1xN for the N IMU epochs -> got qSca %s, angleImu %s, gyroId %s for K = %d, N = %d", ...
        mat2str(size(qSca)), mat2str(size(angleImu)), mat2str(size(gyroId)), K, N);
end
e0 = referenceSystem.daySecondToEpoch([double(epochSca.day(1)) double(epochSca.sec(1))], string(epochSca.scale));
tS = reshape(referenceSystem.secondsBetween(e0, epochSca), 1, []);
if any(diff(tS) <= 0)
    NonGravForces.fail(unit, "inputShape", "the SCA1B epochs must be strictly increasing -> sort them and drop duplicates");
end
imu = epochImu;
if string(imu.scale) ~= string(epochSca.scale)
    imu = referenceSystem.convertTimeScale(imu, string(epochSca.scale));
end
tI = reshape(referenceSystem.secondsBetween(e0, imu), 1, []);

% 1. SCA1B mean body rates over each interval: rotation vector of R_k' R_k+1 over dt
R = referenceSystem.rotationFromQuaternion(qSca);                     % 3x3xK, BODY->GCRS
Mk = pagemtimes(R(:, :, 1:end-1), 'transpose', R(:, :, 2:end), 'none');
w = [squeeze(Mk(3, 2, :) - Mk(2, 3, :))'; squeeze(Mk(1, 3, :) - Mk(3, 1, :))'; squeeze(Mk(2, 1, :) - Mk(1, 2, :))'] / 2;
sn = vecnorm(w);
cs = (squeeze(Mk(1, 1, :) + Mk(2, 2, :) + Mk(3, 3, :))' - 1) / 2;
theta = atan2(sn, cs);
f = ones(1, K - 1);
nz = sn > 0;
f(nz) = theta(nz) ./ sn(nz);
dt = diff(tS);
omega = w .* f ./ dt;
okS = dt <= options.MaxGap_s;
tMid = (tS(1:end-1) + tS(2:end)) / 2;
Q = struct('t_s', tMid, 'omega_rad_s', omega);

% 2. per gyro: accumulated angle at the SCA1B epochs, differenced over the same intervals
gyros = unique(gyroId);
Q.gyros = gyros;
G = numel(gyros);
rec = struct('id', {}, 'rate_rad_s', {}, 'axis', {}, 'scale', {}, 'bias_rad_s', {}, 'condition', {}, 'residual_rad_s', {}, ...
    'rms_rad_s', {}, 'rmsRate_rad_s', {}, 'n', {}, 'cpr', {});
for g = 1:G
    sel = gyroId == gyros(g);
    [tg, order] = sort(tI(sel));
    ag = angleImu(sel);
    ag = ag(order);
    step = diff(ag);
    gap = diff(tg);
    seg = cumsum([1, abs(step) > options.JumpThreshold_rad | gap > options.MaxGap_s | gap <= 0]);
    thA = nan(1, K);
    segA = nan(1, K);
    for s = unique(seg)
        in = seg == s;
        if nnz(in) < 2, continue, end
        cover = tS >= tg(find(in, 1)) & tS <= tg(find(in, 1, 'last'));
        thA(cover) = interp1(tg(in), ag(in), tS(cover), 'linear');
        segA(cover) = s;
    end
    sameSeg = segA(1:end-1) == segA(2:end);
    rate = (thA(2:end) - thA(1:end-1)) ./ dt;
    rate(~(sameSeg & okS)) = NaN;
    use = isfinite(rate) & all(isfinite(omega), 1);
    if nnz(use) < 10
        NonGravForces.fail(unit, "dataOverlap", "gyro %d covers %d SCA1B interval(s) -> the fit needs at least 10; pass IMU1B and SCA1B of the same span", gyros(g), nnz(use));
    end
    X = [omega(:, use)', ones(nnz(use), 1)];
    y = rate(use)';
    nrm = vecnorm(X);
    Xn = X ./ nrm;                                    % column-normalised: the condition the data offer
    [Qx, ~] = qr(Xn, 0);
    res = nan(1, K - 1);
    res(use) = (y - Qx * (Qx' * y))';                 % orthogonal projection: exact however ill-conditioned
    c = lsqminnorm(Xn, y) ./ nrm';
    a = c(1:3);
    r.id = gyros(g);
    r.rate_rad_s = rate;
    r.axis = a / norm(a);
    r.scale = norm(a);
    r.bias_rad_s = c(4);
    r.condition = cond(Xn);
    r.residual_rad_s = res;
    r.rms_rad_s = sqrt(mean(res(use) .^ 2));
    r.rmsRate_rad_s = sqrt(mean(rate(use) .^ 2));
    r.n = nnz(use);
    r.cpr = [NaN NaN];
    rec(g) = r;
end
Q.gyro = rec;

% 3. the nominal frame and the period
period = options.Period_s;
if ~isempty(options.Orbit)
    o = options.Orbit;
    if ~(isstruct(o) && all(isfield(o, ["r" "v"])) && isequal(size(o.r), [3 K]) && isequal(size(o.v), [3 K]))
        NonGravForces.fail(unit, "inputOrbit", "'Orbit' must hold .r and .v, 3x%d GCRS at the SCA1B epochs -> rotate GNV1B with io.readGNV1B and sample it at epochSca", K);
    end
    qN = NonGravForces.geom.AttitudeProvider.nominal(o.r, o.v, options.AlongVelocity);
    D = pagemtimes(referenceSystem.rotationFromQuaternion(qN), 'transpose', R, 'none');   % body in nominal
    wd = [squeeze(D(3, 2, :) - D(2, 3, :))'; squeeze(D(1, 3, :) - D(3, 1, :))'; squeeze(D(2, 1, :) - D(1, 2, :))'] / 2;
    sd = vecnorm(wd);
    th = atan2(sd, (squeeze(D(1, 1, :) + D(2, 2, :) + D(3, 3, :))' - 1) / 2);
    fd = ones(1, K);
    fd(sd > 0) = th(sd > 0) ./ sd(sd > 0);
    ang = wd .* fd;
    if isnan(period)
        gm = referenceSystem.iersConstants().gmEarth;
        aSemi = 1 ./ (2 ./ vecnorm(o.r) - vecnorm(o.v) .^ 2 / gm);
        period = 2 * pi * sqrt(mean(aSemi) ^ 3 / gm);
    end
    Q.nominal = struct('t_s', tS, 'angles_rad', ang, 'mean_rad', mean(ang, 2), ...
        'rms_rad', sqrt(mean((ang - mean(ang, 2)) .^ 2, 2)), 'cpr', localCpr(tS, ang, period));
end
Q.period_s = period;
for g = 1:G
    Q.gyro(g).cpr = localCpr(tMid, Q.gyro(g).residual_rad_s, period);
end
end

function amp = localCpr(t, y, period)
%LOCALCPR Least-squares amplitudes at 1 and 2 cycles per revolution, one row per component of y (NaN rows ignored per column).
amp = nan(size(y, 1), 2);
if ~(isfinite(period) && period > 0), return, end
w = 2 * pi / period;
for i = 1:size(y, 1)
    ok = isfinite(y(i, :));
    if nnz(ok) < 6, continue, end
    tt = t(ok)';
    X = [ones(size(tt)), cos(w * tt), sin(w * tt), cos(2 * w * tt), sin(2 * w * tt)];
    c = lsqminnorm(X, y(i, ok)');                     % a short arc leaves the sinusoids ill-determined
    amp(i, :) = [hypot(c(2), c(3)), hypot(c(4), c(5))];
end
end
