classdef tNgfFigures < NgfTestCase
%TNGFFIGURES Tests of the calibration figures (N-05): app.plotCalibration, app.plotAttribution.
%
%   v3 WP-20: "No precision refval (state so)". What CAN be asserted is that
%   a figure draws its source and nothing else: every tagged curve's data
%   equal the object's own values, bitwise (one authority). The figures run
%   with Visible "off" and are closed by a teardown.
%
% Created 2026-10-15 by Claude (Anthropic) for Matthias Weigelt.
    properties
        obs
        ref
        cal
        info
        T
        A
        ainfo
    end
    methods (TestClassSetup)
        function setup(tc)
            [tc.obs, tc.ref] = NonGravForces.app.syntheticPair('Days', 0.25, 'Seed', 3);
            [tc.cal, tc.info] = NonGravForces.app.fitCalibration(tc.obs, tc.ref, 'ScaleBasis', "linear", 'NBoot', 20);
            D = zeros(3); D(1, 1) = 0.02;
            [o, r] = NonGravForces.app.syntheticPair('Days', 2, 'ScaleDrift', D, 'Seed', 4);
            tc.T = NonGravForces.app.perRevolution(o, r, 'NBoot', 0);
            [tc.A, tc.ainfo] = NonGravForces.app.attributeRevolutions(tc.T, 'NBoot', 200);
        end
    end
    methods
        function f = keep(tc, f)
            %KEEP Close the figure when the test ends.
            tc.addTeardown(@close, f);
        end
    end

    % ===================================================================== contract
    methods (Test, TestTags = {'contract'})
        function plotCalibrationContract(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            [f, files] = NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'Residual', tc.info.residual, ...
                'Visible', "off", 'ExportTo', folder);
            tc.keep(f);
            tc.verifyClass(f, 'matlab.ui.Figure');
            tc.verifyEqual(f.Visible, matlab.lang.OnOffSwitchState.off);
            tc.verifyEqual(files, string(fullfile(folder, "calibration.png")));
            tc.verifyTrue(isfile(files));
            for a = ["x" "y" "z"]
                tc.verifyNotEmpty(findobj(f, 'Tag', "bias_" + a), "bias_" + a);
                tc.verifyNotEmpty(findobj(f, 'Tag', "scale_" + a), "scale_" + a);
                tc.verifyNotEmpty(findobj(f, 'Tag', "acf_" + a), "acf_" + a);
                tc.verifyNotEmpty(findobj(f, 'Tag', "acfWhite_" + a), "acfWhite_" + a);
            end
            tc.verifyEmpty(findobj(f, 'Tag', "scaleBand_y"), 'S_yy is held: no band is drawn for it');
            g = tc.keep(NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'Visible', "off"));
            tc.verifyEmpty(findobj(g, 'Tag', "acf_x"), 'without a residual there is no ACF panel');
            tc.verifyErrorId(@() NonGravForces.app.plotCalibration(5, tc.obs.t), "NonGravForces:plotCalibration:inputCalibration");
            tc.verifyErrorId(@() NonGravForces.app.plotCalibration(tc.cal, double(tc.obs.t)), "NonGravForces:plotCalibration:inputEpochs");
            tc.verifyErrorId(@() NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'Residual', 1), "NonGravForces:plotCalibration:inputResidual");
            tc.verifyErrorId(@() NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'ExportTo', fullfile(folder, "nowhere")), ...
                "NonGravForces:plotCalibration:inputExport");
        end

        function plotAttributionContract(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            [f, files] = NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'Visible', "off", 'ExportTo', folder);
            tc.keep(f);
            tc.verifyClass(f, 'matlab.ui.Figure');
            tc.verifyTrue(isfile(files));
            for j = 1:height(tc.A)
                a = extractBetween(tc.A.parameter(j), 3, 3);
                tc.verifyNotEmpty(findobj(f, 'Tag', "perRev_" + a));
                t = findall(f, 'Tag', "verdict_" + a);   % a title's text hides its handle from findobj
                tc.verifyTrue(contains(string(t.String), tc.A.verdict(j)), 'the banner names the verdict');
            end
            tc.verifyErrorId(@() NonGravForces.app.plotAttribution(tc.T(:, {'window'}), tc.A, tc.ainfo), "NonGravForces:plotAttribution:inputTable");
            bad = tc.ainfo; bad.regression = {};
            tc.verifyErrorId(@() NonGravForces.app.plotAttribution(tc.T, tc.A, bad), "NonGravForces:plotAttribution:inputInfo");
            tc.verifyErrorId(@() NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'ExportTo', fullfile(folder, "nowhere")), ...
                "NonGravForces:plotAttribution:inputExport");
        end
    end

    % ===================================================================== precision
    methods (Test, TestTags = {'precision'})
        function plotCalibrationPrecision(tc)
            % The figure draws the calibration and nothing else: bitwise its methods' values.
            f = tc.keep(NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'Residual', tc.info.residual, 'Visible', "off"));
            S = tc.cal.scaleAt(tc.obs.t); b = tc.cal.biasAt(tc.obs.t);
            ax = ["x" "y" "z"];
            for i = 1:3
                tc.verifyEqual(findobj(f, 'Tag', "scale_" + ax(i)).YData(:), squeeze(S(i, i, :)), "scale_" + ax(i));
                tc.verifyEqual(findobj(f, 'Tag', "bias_" + ax(i)).YData(:), b(i, :).', "bias_" + ax(i));
            end
            rho = timeSeriesPreparation.tsmodels.tsAcf(tc.info.residual.obs(:, 1), 300);
            tc.verifyEqual(findobj(f, 'Tag', "acf_x").YData(:), rho(:), 'the ACF is tsAcf''s, unchanged');
        end

        function plotAttributionPrecision(tc)
            % The per-revolution points are T's, and each fitted line is regressDrivers' coefficients.
            f = tc.keep(NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'Visible', "off"));
            tc.verifyEqual(findobj(f, 'Tag', "perRev_x").YData(:), tc.T.scale(:, 1));
            reg = tc.ainfo.regression{1};
            k = find(reg.driverNames == "time");
            line = findobj(f, 'Tag', "fit_time_x");
            zs = line.XData(:); z = tc.ainfo.drivers.time; ok = isfinite(z) & isfinite(tc.T.scale(:, 1));
            tc.verifyEqual(line.YData(:), reg.coef(1) + reg.coef(k + 1) * (zs - mean(z(ok))), 'the line is c0 + c_k (z - mean z)');
        end
    end

    % ===================================================================== robustness
    methods (Test, TestTags = {'robustness'})
        function plotCalibrationRobustness(tc)
            % A constant calibration with no covariance (a hand-built one) draws flat lines and no bands;
            % a residual with a gap still whitens, run by run.
            c = NonGravForces.app.Calibration(struct('scale', eye(3), 'bias', struct('basis', "constant", ...
                't0', tc.obs.t(1), 'knots', int64([]), 'coef', zeros(3, 1))));
            f = tc.keep(NonGravForces.app.plotCalibration(c, tc.obs.t, 'Visible', "off"));
            tc.verifyEqual(unique(findobj(f, 'Tag', "scale_x").YData), 1);
            tc.verifyEmpty(findobj(f, 'Tag', "biasBand_x"));
            r = tc.info.residual; y = r.obs; y(5000:5300, :) = NaN;
            g = tc.keep(NonGravForces.app.plotCalibration(tc.cal, tc.obs.t, 'Residual', r.derive(r.t, y, "test.gap"), 'Visible', "off"));
            tc.verifyTrue(all(isfinite(findobj(g, 'Tag', "acfWhite_x").YData(1:10))));
        end

        function plotAttributionRobustness(tc)
            % A gated revolution is drawn as a gap (NaN), never as a point; an unresolved parameter
            % (no regression) still gets its time panel and banner.
            T = tc.T; T.scale(3, :) = NaN;
            [Ag, ig] = NonGravForces.app.attributeRevolutions(T, 'NBoot', 100);
            f = tc.keep(NonGravForces.app.plotAttribution(T, Ag, ig, 'Visible', "off"));
            tc.verifyTrue(isnan(findobj(f, 'Tag', "perRev_x").YData(3)));
            A2 = Ag; A2.verdict(:) = "unresolved"; i2 = ig; i2.regression = cell(1, height(Ag));
            g = tc.keep(NonGravForces.app.plotAttribution(T, A2, i2, 'Visible', "off"));
            tc.verifyNotEmpty(findobj(g, 'Tag', "perRev_x"));
            tc.verifyEmpty(findobj(g, '-regexp', 'Tag', '^fit_'));
        end
    end

    % ===================================================================== speed
    methods (Test, TestTags = {'speed'})
        function plotCalibrationSpeed(tc)
            % Curves of twice the epochs: about twice the drawing, the fixed figure cost pulling it lower.
            t1 = tc.obs.t(1:2:end); t2 = tc.obs.t;
            close(NonGravForces.app.plotCalibration(tc.cal, t1, 'Visible', "off")); close(NonGravForces.app.plotCalibration(tc.cal, t2, 'Visible', "off"));
            d1 = zeros(1, 3); d2 = zeros(1, 3);
            for k = 1:3
                t0 = tic; f = NonGravForces.app.plotCalibration(tc.cal, t1, 'Visible', "off"); drawnow; d1(k) = toc(t0); close(f);
                t0 = tic; f = NonGravForces.app.plotCalibration(tc.cal, t2, 'Visible', "off"); drawnow; d2(k) = toc(t0); close(f);
            end
            tc.verifyBoundAndRecord("plotCalibration, 2x epochs / 1x", median(d2) / median(d1), 4, "-", "linear would be 2");
        end

        function plotAttributionSpeed(tc)
            % One driver panel against three: the panels dominate, so about three times, never quadratic.
            close(NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'Visible', "off", 'MaxDrivers', 1));
            d1 = zeros(1, 3); d3 = zeros(1, 3);
            for k = 1:3
                t0 = tic; f = NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'Visible', "off", 'MaxDrivers', 1); drawnow; d1(k) = toc(t0); close(f);
                t0 = tic; f = NonGravForces.app.plotAttribution(tc.T, tc.A, tc.ainfo, 'Visible', "off", 'MaxDrivers', 3); drawnow; d3(k) = toc(t0); close(f);
            end
            tc.verifyBoundAndRecord("plotAttribution, 3 driver panels / 1", median(d3) / median(d1), 6, "-", "panels are drawn once each: about 3 at most");
        end
    end
end
