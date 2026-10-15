function [f, files] = plotAttribution(T, A, info, options)
%PLOTATTRIBUTION Figure of the per-revolution attribution: scales over time, scales against drivers, verdict.
%
%   f = NonGravForces.app.plotAttribution(T, A, info)
%   [f, files] = NonGravForces.app.plotAttribution(T, A, info, 'MaxDrivers', 3, 'Visible', "off", 'ExportTo', folder)
%
%   The attribution figure of N-05 (the v3 plan's WP-20 `attribution`). One
%   row per attributed parameter (A's rows):
%     first   the per-revolution scale over time with its 1-sigma bars
%             (app.perRevolution's T); a gated revolution is a gap, never a
%             point
%     then    the scale against each of its 'MaxDrivers' strongest drivers (by
%             bootstrap |coef / sigma|), with the regression's partial line
%             c0 + c_k (z_k - mean z_k) and its 1-sigma slope band from the
%             bootstrap sigma of c_k
%   and the verdict with its reason as a banner over the row. Everything drawn
%   is read from T, A and info.regression (timeSeriesPreparation.
%   regressDrivers' result): the figure fits nothing. Each drawn object
%   carries a Tag ("perRev_x", "driver_<name>_x", "fit_<name>_x",
%   "verdict_x"), which is how a test reads it back.
%
% INPUT ARGUMENTS
%   Name  Size  Type    Units  Description
%   ----  ----  ----    -----  -----------
%   T     RxN   table   -      app.perRevolution's table (tMid, scale, optional sigmaScale)
%   A     Jx6   table   -      app.attributeRevolutions' verdicts
%   info  1x1   struct  -      app.attributeRevolutions' info (.drivers, .regression)
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units  Description
%   ----        ----  ----    -----  -----------
%   MaxDrivers  1x1   double  -      Driver panels per parameter, the strongest first; default 3
%   Visible     1x1   string  -      "on" (default) or "off"
%   ExportTo    1x1   string  -      Folder: writes attribution.png at 200 dpi; default ""
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   f      1x1   figure  -      The figure
%   files  1xF   string  -      Files written
%
% ACCURACY
%   Not applicable: the figure draws regressDrivers' coefficients and sigmas
%   and perRevolution's estimates, each with its unit's accuracy. The partial
%   line shows one coefficient with the other drivers at their means; with
%   correlated drivers the scatter around it is not the regression's residual.
%
% ERRORS
%   NonGravForces:plotAttribution:inputTable   T lacks tMid or scale, or A lacks parameter or verdict
%   NonGravForces:plotAttribution:inputInfo    info lacks drivers or regression, or does not match A
%   NonGravForces:plotAttribution:inputExport  'ExportTo' is not an existing folder
%
% REFERENCES
%   The v3 plan (timeSeriesPreparation docs/inputs/HANDOVER_AccelCal_v3.md)
%   WP-20; timeSeriesPreparation docs/CONCEPT.md 6.2, 7 and 10.1 (N-05).
%
% EXAMPLE
%   D = zeros(3); D(1, 1) = 0.02;
%   [obs, ref] = NonGravForces.app.syntheticPair('Days', 2, 'ScaleDrift', D);   % >= 20 revolutions
%   T = NonGravForces.app.perRevolution(obs, ref, 'NBoot', 0);
%   [A, info] = NonGravForces.app.attributeRevolutions(T, 'NBoot', 200);
%   f = NonGravForces.app.plotAttribution(T, A, info);
%
% See also NonGravForces.app.attributeRevolutions, NonGravForces.app.perRevolution, NonGravForces.app.plotCalibration.
%
% KEYWORDS attribution figure plot per-revolution drivers verdict regression
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
arguments
    T table
    A table
    info (1,1) struct
    options.MaxDrivers (1,1) double {mustBeInteger, mustBePositive} = 3
    options.Visible (1,1) string {mustBeMember(options.Visible, ["on" "off"])} = "on"
    options.ExportTo (1,1) string = ""
end
unit = "plotAttribution";
if ~all(ismember(["tMid" "scale"], string(T.Properties.VariableNames))) || ~all(ismember(["parameter" "verdict"], string(A.Properties.VariableNames)))
    NonGravForces.fail(unit, "inputTable", "T must be app.perRevolution's table and A app.attributeRevolutions' -> check the arguments' order");
end
if ~all(isfield(info, {'drivers', 'regression'})) || numel(info.regression) ~= height(A)
    NonGravForces.fail(unit, "inputInfo", "info must be app.attributeRevolutions' second output for the same A -> it holds %d regressions for %d parameters", ...
        numel(info.regression), height(A));
end
if strlength(options.ExportTo) > 0 && ~isfolder(options.ExportTo)
    NonGravForces.fail(unit, "inputExport", "'ExportTo' must be an existing folder -> %s does not exist", options.ExportTo);
end

hours = double(T.tMid - T.tMid(1)) / 3600e9;
J = height(A);
ax = ["x" "y" "z"];
nCol = 1 + options.MaxDrivers;
f = figure("Visible", options.Visible, "Name", "attribution", "Color", "w", "Position", [100 100 300 * nCol 300 * max(J, 1)]);
tl = tiledlayout(f, max(J, 1), nCol, "TileSpacing", "compact");
title(tl, "Per-revolution attribution");
for j = 1:J
    i = find(ax == extractBetween(A.parameter(j), 3, 3));
    P = T.scale(:, i);
    a1 = nexttile(tl);
    if ismember("sigmaScale", string(T.Properties.VariableNames))
        errorbar(a1, hours, P, T.sigmaScale(:, i), "bo", "Tag", "perRev_" + ax(i));
    else
        plot(a1, hours, P, "bo", "Tag", "perRev_" + ax(i));
    end
    banner = A.parameter(j) + ": " + A.verdict(j);
    r = A.reason(j);
    if strlength(r) > 60, r = extractBefore(r, 61) + "..."; end
    if strlength(r) > 0, banner = banner + " (" + r + ")"; end
    title(a1, banner, "Tag", "verdict_" + ax(i), "Interpreter", "none");
    xlabel(a1, "hours"); ylabel(a1, "scale");
    reg = info.regression{j};
    shown = 0;
    if ~isempty(reg)
        z = abs(reg.coef(2:end) ./ reg.sigma(2:end));
        [~, order] = sort(z(:), 'descend');
        for k = order(1:min(options.MaxDrivers, numel(order))).'
            name = reg.driverNames(k);
            if startsWith(name, "amp_")
                zk = T.(name)(:, i);                % an effect amplitude is per axis, as attributed
            else
                zk = info.drivers.(name);
            end
            a = nexttile(tl);
            plot(a, zk, P, "k.", "MarkerSize", 12, "Tag", "driver_" + name + "_" + ax(i)); hold(a, "on");
            ok = isfinite(zk) & isfinite(P);
            zs = linspace(min(zk(ok)), max(zk(ok)), 50)';
            m = mean(zk(ok));
            c0 = reg.coef(1); ck = reg.coef(k + 1); sk = reg.sigma(k + 1);
            fill(a, [zs; flipud(zs)], [c0 + (ck + sk) * (zs - m); flipud(c0 + (ck - sk) * (zs - m))], [0.85 0.85 0.85], ...
                "EdgeColor", "none", "FaceAlpha", 0.6);
            plot(a, zs, c0 + ck * (zs - m), "r-", "LineWidth", 1.2, "Tag", "fit_" + name + "_" + ax(i));
            hold(a, "off");
            title(a, sprintf("%s   |z| = %.1f", name, z(k)), "Interpreter", "none"); xlabel(a, name, "Interpreter", "none");
            shown = shown + 1;
        end
    end
    for k = shown + 1:options.MaxDrivers
        axis(nexttile(tl), "off");
    end
end
files = strings(1, 0);
if strlength(options.ExportTo) > 0
    files = fullfile(options.ExportTo, "attribution.png");
    exportgraphics(f, files, "Resolution", 200);
end
end
