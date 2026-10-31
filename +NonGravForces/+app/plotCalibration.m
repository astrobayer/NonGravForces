function [f, files] = plotCalibration(cal, t, options)
%PLOTCALIBRATION Figure of an accelerometer calibration: scale and bias with their bands, residual ACF.
%
%   f = NonGravForces.app.plotCalibration(cal, t)
%   [f, files] = NonGravForces.app.plotCalibration(cal, t, 'Residual', info.residual, 'Visible', "off", 'ExportTo', folder)
%
%   The calibration figure of N-05 (the v3 plan's WP-20 `calibration`). One
%   row per axis (every axis carries an estimated bias):
%     left    S_ii(t), the scale the calibration carries, with its pointwise
%             1-sigma band -- a curve for a time-varying scale (N-03a), a flat
%             line for a constant one; the drift test's p-value in the title
%     middle  b_i(t), the bias curve, with the band sqrt(diag(B C B')) of its
%             basis through its covariance (jump steps included)
%     right   with 'Residual': the residual's sample ACF and, after an AR('Order')
%             whitening, its innovations' ACF, with the 95 % white-noise band.
%             The AR model is the closed-form Yule-Walker estimate
%             (tsFitAR, Method "yw") on the longest gap-free run, applied as
%             a filter to every run; tsFitAR's exact Kalman method accepts
%             gaps but took over a minute on 4000 samples, too slow to draw.
%             The first shows why the formal sigmas lie (N-03a measured the
%             bootstrap sigma at 26 times the formal one); the second whether a
%             low-order AR model would have explained the colour. An axis
%             with no more than 2 'MaxLag' finite residual samples gets an
%             empty ACF panel (its title only).
%   Everything drawn comes from the calibration's own methods (scaleAt,
%   scaleSigmaAt, biasAt, basis) and from timeSeriesPreparation.tsmodels: the
%   figure computes nothing of its own, so it cannot disagree with the object.
%   Each drawn object carries a Tag ("scale_x", "scaleBand_x", "bias_x",
%   "biasBand_x", "acf_x", "acfWhite_x"), which is how a test reads it back.
%
% INPUT ARGUMENTS
%   Name  Size  Type                            Units  Description
%   ----  ----  ----                            -----  -----------
%   cal   1x1   NonGravForces.app.Calibration   -      The calibration to draw
%   t     Kx1   int64                           ns     Epochs to draw the curves at, GPS ns
%
% NAME-VALUE ARGUMENTS
%   Name      Size  Type                        Units  Description
%   ----      ----  ----                        -----  -----------
%   Residual  1x1   timeSeriesPreparation.Arc   m/s^2  app.fitCalibration's info.residual; default none
%   Order     1x1   double                      -      AR order of the whitening; default 5
%   MaxLag    1x1   double                      -      ACF lags drawn; default 300 (samples)
%   Visible   1x1   string                      -      "on" (default) or "off"
%   ExportTo  1x1   string                      -      Folder: writes calibration.png at 200 dpi; default ""
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   f      1x1   figure  -      The figure
%   files  1xF   string  -      Files written: calibration.png in 'ExportTo', or a 1x0
%                               string without it
%
% ACCURACY
%   Not applicable: a figure draws numbers computed elsewhere, each with its
%   unit's accuracy. The bias band is the 1-sigma pointwise band of the
%   covariance the calibration carries (bootstrap or formal, as fitted); the
%   ACF band is tsAcf's 1.96 / sqrt(n), which assumes white noise.
%
% ERRORS
%   NonGravForces:plotCalibration:inputCalibration  cal is not an app.Calibration
%   NonGravForces:plotCalibration:inputEpochs       t is not a non-empty int64 vector
%   NonGravForces:plotCalibration:inputResidual     'Residual' is not a 3-column Arc
%   NonGravForces:plotCalibration:inputExport       'ExportTo' is not an existing folder
%   timeSeriesPreparation:*                         from tsAcf and tsFitAR (Yule-Walker)
%
% REFERENCES
%   The v3 plan (timeSeriesPreparation docs/inputs/HANDOVER_AccelCal_v3.md)
%   WP-20; timeSeriesPreparation docs/CONCEPT.md 10.1 (N-05).
%
% EXAMPLE
%   [obs, ref] = NonGravForces.app.syntheticPair('Days', 0.25);
%   [cal, info] = NonGravForces.app.fitCalibration(obs, ref, 'ScaleBasis', "linear", 'NBoot', 20);
%   f = NonGravForces.app.plotCalibration(cal, obs.t, 'Residual', info.residual);
%
% See also NonGravForces.app.Calibration, NonGravForces.app.fitCalibration, NonGravForces.app.plotAttribution, timeSeriesPreparation.report.
%
% KEYWORDS calibration figure plot scale bias band ACF whitening
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    cal (1,1)
    t
    options.Residual = []
    options.Order (1,1) double {mustBeInteger, mustBePositive} = 5
    options.MaxLag (1,1) double {mustBeInteger, mustBePositive} = 300
    options.Visible (1,1) string {mustBeMember(options.Visible, ["on" "off"])} = "on"
    options.ExportTo (1,1) string = ""
end
unit = "plotCalibration";
if ~isa(cal, 'NonGravForces.app.Calibration')
    NonGravForces.fail(unit, "inputCalibration", "cal must be an app.Calibration -> got a %s", class(cal));
end
if ~isa(t, 'int64') || isempty(t) || ~isvector(t)
    NonGravForces.fail(unit, "inputEpochs", "t must be a non-empty int64 vector of GPS ns -> got a %s %s", class(t), mat2str(size(t)));
end
res = options.Residual;
if ~isempty(res) && (~isa(res, 'timeSeriesPreparation.Arc') || size(res.obs, 2) ~= 3)
    NonGravForces.fail(unit, "inputResidual", "'Residual' must be app.fitCalibration's info.residual, a 3-column Arc -> got a %s", class(res));
end
if strlength(options.ExportTo) > 0 && ~isfolder(options.ExportTo)
    NonGravForces.fail(unit, "inputExport", "'ExportTo' must be an existing folder -> %s does not exist", options.ExportTo);
end

t = t(:);
hours = double(t - t(1)) / 3600e9;
S = cal.scaleAt(t); sS = cal.scaleSigmaAt(t);
b = cal.biasAt(t);
B = NonGravForces.app.Calibration.basis(cal.bias, t);
d = cal.drift();
ax = ["x" "y" "z"];
rows = 1:3;                                         % every axis carries an estimated bias
nCol = 2 + ~isempty(res);
f = figure("Visible", options.Visible, "Name", "calibration", "Color", "w", "Position", [100 100 1200 260 * numel(rows)]);
tl = tiledlayout(f, numel(rows), nCol, "TileSpacing", "compact");
title(tl, sprintf("Calibration (%s convention), %s scale, %s bias", cal.convention, ...
    cal.scaleCurve.basis, cal.bias.basis));
for i = rows
    % scale
    a1 = nexttile(tl);
    s = squeeze(S(i, i, :)); sg = squeeze(sS(i, i, :));
    if all(isfinite(sg))
        fill(a1, [hours; flipud(hours)], [s + sg; flipud(s - sg)], [0.8 0.85 1], "EdgeColor", "none", "Tag", "scaleBand_" + ax(i));
        hold(a1, "on");
    end
    plot(a1, hours, s, "b-", "LineWidth", 1.2, "Tag", "scale_" + ax(i));
    hold(a1, "off");
    ttl = "S_{" + ax(i) + ax(i) + "}";
    if ~isnan(d.pValue(i, i)), ttl = ttl + sprintf("   drift p = %.3g", d.pValue(i, i)); end
    if ~any(isfinite(sg)), ttl = ttl + "   (held)"; end
    title(a1, ttl); xlabel(a1, "hours"); ylabel(a1, "-");
    % bias
    a2 = nexttile(tl);
    bi = b(i, :).';
    C = cal.bias.cov{i};
    if ~isempty(C)
        sb = sqrt(max(sum((B * C) .* B, 2), 0));
        fill(a2, [hours; flipud(hours)], [bi + sb; flipud(bi - sb)], [1 0.85 0.8], "EdgeColor", "none", "Tag", "biasBand_" + ax(i));
        hold(a2, "on");
    end
    plot(a2, hours, bi, "r-", "LineWidth", 1.2, "Tag", "bias_" + ax(i));
    hold(a2, "off");
    title(a2, "b_" + ax(i)); xlabel(a2, "hours"); ylabel(a2, "m/s^2");
    % residual ACF before and after whitening
    if ~isempty(res)
        a3 = nexttile(tl);
        x = res.obs(:, i);
        if nnz(isfinite(x)) > 2 * options.MaxLag
            [rho, lags, conf] = timeSeriesPreparation.tsmodels.tsAcf(x, options.MaxLag);
            rhoW = timeSeriesPreparation.tsmodels.tsAcf(localWhiten(x, options.Order), options.MaxLag);
            plot(a3, lags, rho, "k-", "Tag", "acf_" + ax(i)); hold(a3, "on");
            plot(a3, lags, rhoW, "g-", "Tag", "acfWhite_" + ax(i));
            yline(a3, [conf -conf], "k:");
            hold(a3, "off");
            legend(a3, ["residual" sprintf("AR(%d) innovations", options.Order)], "Location", "northeast");
        end
        title(a3, "residual ACF, axis " + ax(i)); xlabel(a3, "lag (samples)");
    end
end
files = strings(1, 0);
if strlength(options.ExportTo) > 0
    files = fullfile(options.ExportTo, "calibration.png");
    exportgraphics(f, files, "Resolution", 200);
end
end

% =========================================================================================
function e = localWhiten(x, p)
%LOCALWHITEN Innovations of an AR(p) fitted by Yule-Walker on the longest gap-free run.
ok = isfinite(x);
edges = diff([0; ok; 0]);
s = find(edges == 1); k = find(edges == -1) - 1;
[~, longest] = max(k - s);
m = timeSeriesPreparation.tsmodels.tsFitAR(x(s(longest):k(longest)), p, Method = "yw");
e = nan(size(x));
for r = 1:numel(s)
    seg = x(s(r):k(r)) - m.mu;
    w = filter([1, -m.ar], 1, seg);
    w(1:min(p, numel(w))) = NaN;                       % the filter's start-up, not an innovation
    e(s(r):k(r)) = w;
end
end
