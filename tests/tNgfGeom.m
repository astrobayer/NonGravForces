classdef tNgfGeom < NgfTestCase
%TNGFGEOM Tests of geom.MacroModel and geom.graceFO.
% Created 2026-08-27 by Claude Fable 5 (Anthropic).
    methods (Test, TestTags = {'contract'})
        function macroModelContract(tc)
            mm = NonGravForces.geom.graceFO("C");
            tc.verifyClass(mm, 'NonGravForces.geom.MacroModel');
            tc.verifyEqual(mm.n, 8); tc.verifySize(mm.N, [3 8]); tc.verifySize(mm.optVis, [3 8]);
            p = mm.panels; p(3).specVis = 0.5;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputEnergy');
            p = mm.panels; p(1).normal = [2;0;0];
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputNormal');
            p = rmfield(mm.panels, 'group');
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(p), 'NonGravForces:MacroModel:inputPanels');
            % the per-panel alpha is no longer read (O-23 closed by removal): a panel table carrying it is still accepted
            pa = mm.panels; [pa.alpha] = deal(0.85); tc.verifyEqual(NonGravForces.geom.MacroModel(pa).A, mm.A);
            % mass and internal shares (LP-C F10): graceFO carries the handbook mass, a share must cover every panel
            tc.verifyEqual(mm.mass_kg, NonGravForces.const().graceFO.massC); tc.verifyTrue(all(isnan(mm.internalShare)));
            ps = mm.panels; [ps.internalShare] = deal(0); ps(7).internalShare = 1;
            tc.verifyEqual(NonGravForces.geom.MacroModel(ps).internalShare, [0 0 0 0 0 0 1 0]);
            ps(1).internalShare = 0.5;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(ps), 'NonGravForces:MacroModel:inputPanels');
            ps = mm.panels; [ps.internalShare] = deal([]); ps(7).internalShare = 1;
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(ps), 'NonGravForces:MacroModel:inputPanels');
            tc.verifyErrorId(@() NonGravForces.geom.MacroModel(mm.panels, 'Mass_kg', -1), 'NonGravForces:MacroModel:inputPanels');
        end
        function graceFOContract(tc)
            tc.verifyError(@() NonGravForces.geom.graceFO("E"), 'MATLAB:validators:mustBeMember');
            d = NonGravForces.geom.graceFO("D");
            tc.verifyEqual(d.name, "GRACE-D");
            tc.verifyEqual(d.inertia(2,2), 579.759);
        end
    end
    methods (Test, TestTags = {'reference'})
        function graceFOReference(tc)
            % JPL D-56935 Table 5: total area and the nadir-plate optics.
            mm = NonGravForces.geom.graceFO("C");
            % Sum of the eight printed areas: 2*(0.9551567+3.1554792+0.2282913)+6.0711120+2.1673620.
            tc.verifyRelClose(sum(mm.A), 16.9163284, 1e-7, "graceFO total area vs Table 5 sum");
            k = find([mm.panels.name] == "Nadir");
            tc.verifyEqual(mm.optVis(:,k), [0.12; 0.68; 0.20]);
            tc.verifyEqual(mm.optIR(:,k), [0.75; 0.19; 0.06]);
            tc.verifyEqual(mm.inertia(3,3), 649.690);
        end
        function macroModelReference(tc)
            % Side normals use sin/cos 50 deg as printed (panels), renormalised to unit length in N (O-2, 2026-09-07).
            mm = NonGravForces.geom.graceFO("C");
            k = find([mm.panels.name] == "StarboardOuter");
            tc.verifyEqual(mm.panels(k).normal, [0; 0.766044; -0.642787]);
            tc.verifyEqual(norm(mm.N(:,k)), 1, 'AbsTol', 2e-16);
            tc.verifyBoundAndRecord("renormalised normal vs printed digits", norm(mm.N(:,k) - [0; 0.766044; -0.642787]), 1e-6, "-", "six printed digits: their norm is 1 - 7.3e-7 (measured 2026-09-07; the record had declared 5e-7)");
        end
    end
    methods (Test, TestTags = {'precision'})
        function graceFOPrecision(tc)
            % Non-convexity: the silhouette defect must equal 2*A_inner*c50 in Z
            % (inner and outer flank normals cancel in Y). Bound 1e-9 m^2: arithmetic.
            % c50 is taken from the packed unit normal (O-2, 2026-09-07), not from the printed 0.642787 (norm 1 - 7.3e-7).
            mm = NonGravForces.geom.graceFO("C");
            s = mm.silhouetteDefect();
            c50 = abs(mm.N(3, find([mm.panels.name] == "StarboardOuter", 1)));
            expectZ = 2 * 0.2282913 * c50 - 2 * 3.1554792 * c50 + 6.0711120 - 2.1673620;
            tc.verifyBoundAndRecord("graceFO silhouette defect z", abs(s(3) - expectZ), 1e-9, "m^2");
            tc.verifyBoundAndRecord("graceFO silhouette defect xy", norm(s(1:2)), 1e-9, "m^2");
        end
        function macroModelPrecision(tc)
            % Energy closure is enforced at 1e-6; the handbook numbers close to 1e-12.
            mm = NonGravForces.geom.graceFO("C");
            tc.verifyBoundAndRecord("MacroModel energy closure vis", max(abs(sum(mm.optVis,1)-1)), 1e-12, "-");
            tc.verifyBoundAndRecord("MacroModel energy closure ir",  max(abs(sum(mm.optIR,1)-1)),  1e-12, "-");
        end
    end
    methods (Test, TestTags = {'robustness'})
        function macroModelRobustness(tc)
            % Unknown centroids (all NaN) are accepted; they only disable torques.
            mm = NonGravForces.geom.graceFO("C"); p = mm.panels;
            for k = 1:numel(p), p(k).centroid_m = nan(3,1); end
            m2 = NonGravForces.geom.MacroModel(p);
            tc.verifyTrue(all(isnan(m2.C(:))));
        end
        function graceFORobustness(tc)
            % Non-operating array closes the energy budget by construction.
            mm = NonGravForces.geom.graceFO("C", 'SolarArrayOperating', false);
            tc.verifyLessThan(max(abs(sum(mm.optVis,1)-1)), 1e-12);
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function macroModelMetamorphic(tc)
            % Rotating every normal and centroid by R rotates the silhouette defect by R.
            mm = NonGravForces.geom.graceFO("C"); p = mm.panels;
            [Q, ~] = qr(reshape(sin(1:9), 3, 3)); R = Q * diag(sign(diag(Q))); if det(R) < 0, R(:,1) = -R(:,1); end
            for k = 1:numel(p), p(k).normal = R * p(k).normal; p(k).centroid_m = R * p(k).centroid_m; p(k).vertices_m = R * p(k).vertices_m; end
            m2 = NonGravForces.geom.MacroModel(p);
            tc.verifyBoundAndRecord("MacroModel rotation covariance", norm(m2.silhouetteDefect() - R * mm.silhouetteDefect()), 1e-12, "m^2");
        end
        function graceFOMetamorphic(tc)
            % C and D differ only in mass properties, never in geometry.
            c = NonGravForces.geom.graceFO("C"); d = NonGravForces.geom.graceFO("D");
            tc.verifyEqual(c.A, d.A); tc.verifyEqual(c.N, d.N); tc.verifyEqual(c.optVis, d.optVis);
            tc.verifyNotEqual(c.inertia, d.inertia);
        end
    end
end
