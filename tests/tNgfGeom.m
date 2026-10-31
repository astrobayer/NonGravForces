classdef tNgfGeom < NgfTestCase
%TNGFGEOM Tests of geom.MacroModel, geom.graceFO and geom.grace.
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
            % The silhouette defect is the sum of the PRINTED area vectors, which do not close on each other (a closed body,
            % convex or not, gives zero): in Z the flanks' 2*A*c50 terms plus nadir minus zenith; in Y the inner and
            % outer flank normals cancel. Bound 1e-9 m^2: arithmetic.
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
    % geom.grace (N-07, 2026-09-22): every expected number below is read from R-5, Bettadpur (2012),
    % GRACE 327-720 Rev 4.6, and cited by section; none is taken from geom.grace itself.
    methods (Test, TestTags = {'contract'})
        function graceContract(tc)
            tc.verifyError(@() NonGravForces.geom.grace("C"), 'MATLAB:validators:mustBeMember');
            b = NonGravForces.geom.grace("B");
            tc.verifyClass(b, 'NonGravForces.geom.MacroModel');
            tc.verifyEqual(b.name, "GRACE-B"); tc.verifyEqual(b.n, 8); tc.verifyTrue(b.hasVertices);
            % R-5 prints no inertia and no CoG -> proof-mass offset: unknown, never zero
            tc.verifyTrue(all(isnan(b.inertia(:)))); tc.verifyTrue(all(isnan(b.proofMassOffset)));
            tc.verifyEqual(b.centreOfMass, zeros(3, 1));
        end
    end
    methods (Test, TestTags = {'reference'})
        function graceReference(tc)
            % R-5 Sec. I.5.3.3 (launch mass) and I.5.3.5 (surface table).
            a = NonGravForces.geom.grace("A");
            tc.verifyEqual(a.mass_kg, 487.2); tc.verifyEqual(NonGravForces.geom.grace("B").mass_kg, 487.2);
            tc.verifyRelClose(sum(a.A), 2*(0.9551567+3.1554792+0.2282913)+6.0711120+2.1673620, 1e-12, "grace total area vs Sec. I.5.3.5");
            k = find([a.panels.name] == "Nadir");
            tc.verifyEqual(a.optVis(:,k), [0.12; 0.68; 0.20]); tc.verifyEqual(a.optIR(:,k), [0.75; 0.19; 0.06]);
            k = find([a.panels.name] == "PortOuter");
            tc.verifyEqual(a.panels(k).normal, [0; -0.766044; -0.642787]);
            tc.verifyEqual(a.optVis(:,k), [0.65; 0.05; 0.30]); tc.verifyEqual(a.optIR(:,k), [0.81; 0.03; 0.16]);
        end
        function graceReferenceDimensions(tc)
            % The polygons against numbers R-5 prints and geom.grace does not use:
            %  - Fig 3's 56 mm: depth of the overhang below the nadir plane (geom.grace uses the
            %    printed underside AREA, 0.2282913 m^2, and the 50 deg slant);
            %  - Sec. I.5.3.2: the orbit thrusters O11/O21 at X = -1561 mm on the rear face and the
            %    attitude thrusters A14/A22 at Z = +275 mm on the nadir face, A12/A24 at Z = -444 mm;
            %  - Fig 4's projected side area, 2.4633921 m^2, which includes the boom (0.0461901 m^2,
            %    Sec. I.5.3.5 note 4) and the dome, notch and recesses the plates do not model.
            % Dimensions are printed to the millimetre, so 1 mm (plus the half-millimetre of the
            % odd 3123 mm length) is the bound; the side-area bound is the 1 % below which the
            % unmodelled features lie (dome 46 mm x ~0.1 m, notch 38 mm x ~0.1 m).
            a = NonGravForces.geom.grace("A");
            V = [a.V{:}];
            zNadir = a.C(3, [a.panels.name] == "Nadir"); zZen = a.C(3, [a.panels.name] == "Zenith");
            tc.verifyBoundAndRecord("grace overhang depth vs Fig 3 (56 mm)", abs((max(V(3,:)) - zNadir) - 0.056), 1e-3, "m");
            tc.verifyBoundAndRecord("grace nadir plane vs thrusters A14/A22 (Z = +275 mm)", abs(zNadir - 0.275), 1e-3, "m");
            % The zenith face (720 - 275 = 445 mm, from two printed values) and the thruster (444 mm) are
            % each rounded to the millimetre, so rounding alone allows 1.5 mm between them. The thruster
            % must lie INSIDE the body; the gap is bounded at 2 mm. (Until 2026-09-22 this asserted the gap
            % at 1 mm + 1e-12 and measured exactly 1.0000e-03: a bound that passed on its tolerance only.)
            tc.verifyLessThan(zZen, -0.444, "thrusters A12/A24 must sit inside the zenith face");
            tc.verifyBoundAndRecord("grace zenith plane vs thrusters A12/A24 (Z = -444 mm)", abs(zZen + 0.444), 2e-3, "m", "measured 1 mm on CI run 35773012179; millimetre rounding of three printed values allows 1.5 mm");
            tc.verifyBoundAndRecord("grace rear face vs orbit thrusters O11/O21 (X = -1561 mm)", abs(min(V(1,:)) + 1.561), 1e-3, "m");
            side = (max(V(3,:)) - min(V(3,:))) * (max(V(1,:)) - min(V(1,:))) + 0.0461901;   % prismatic body: silhouette from the side
            tc.verifyBoundAndRecord("grace side silhouette vs Fig 4 C.o.A. area", abs(side / 2.4633921 - 1), 1e-2, "-", "measured 2.53e-3 in Python on 2026-09-22 before any MATLAB run");
        end
    end
    methods (Test, TestTags = {'precision'})
        function gracePrecision(tc)
            % Same identity as graceFOPrecision, on R-5's table: inner and outer flank normals cancel in Y,
            % the Z defect is 2*A_inner*c50 - 2*A_outer*c50 + A_nadir - A_zenith. Arithmetic: 1e-9 m^2.
            a = NonGravForces.geom.grace("A");
            s = a.silhouetteDefect();
            c50 = abs(a.N(3, find([a.panels.name] == "StarboardOuter", 1)));
            expectZ = 2 * 0.2282913 * c50 - 2 * 3.1554792 * c50 + 6.0711120 - 2.1673620;
            tc.verifyBoundAndRecord("grace silhouette defect z", abs(s(3) - expectZ), 1e-9, "m^2");
            tc.verifyBoundAndRecord("grace silhouette defect xy", norm(s(1:2)), 1e-9, "m^2");
            tc.verifyBoundAndRecord("grace energy closure", max(abs([sum(a.optVis,1), sum(a.optIR,1)] - 1)), 1e-12, "-");
        end
    end
    methods (Test, TestTags = {'robustness'})
        function graceRobustness(tc)
            % Note 2's non-operating array closes by construction (the printed 0.72/0.05/0.30 would sum to 1.07).
            a = NonGravForces.geom.grace("A", 'SolarArrayOperating', false);
            tc.verifyLessThan(max(abs(sum(a.optVis,1) - 1)), 1e-12);
            k = find([a.panels.name] == "Zenith"); tc.verifyEqual(a.optVis(:,k), [0.72; 0.05; 0.23]);
            % The polygons must not intersect: nothing occludes the nadir plate from straight below or
            % the zenith plate from straight above (GRACE-FO's symptom of a protruding plate, 2026-09-07).
            f = NonGravForces.geom.selfShadow(a, [[0; 0; 1], [0; 0; -1]], 'Samples', 20);
            tc.verifyEqual(f([a.panels.name] == "Nadir", 1), 1);
            tc.verifyEqual(f([a.panels.name] == "Zenith", 2), 1);
        end
    end
    methods (Test, TestTags = {'metamorphic'})
        function graceMetamorphic(tc)
            % A and B are one table (R-5 Sec. I.5.3.3): everything but the name is identical.
            a = NonGravForces.geom.grace("A"); b = NonGravForces.geom.grace("B");
            tc.verifyEqual(a.A, b.A); tc.verifyEqual(a.N, b.N); tc.verifyEqual(a.C, b.C);
            tc.verifyEqual(a.optVis, b.optVis); tc.verifyEqual(a.optIR, b.optIR); tc.verifyEqual(a.mass_kg, b.mass_kg);
            tc.verifyNotEqual(a.name, b.name);
        end
        function graceMetamorphicAcrossDocuments(tc)
            % R-5 Sec. I.5.3.5 and JPL D-56935 Table 5 print the same surface table. The two factories
            % transcribe them independently, so this compares two transcriptions of one fact; the
            % geometry, mass and origin differ by document and are NOT compared.
            a = NonGravForces.geom.grace("A"); c = NonGravForces.geom.graceFO("C");
            tc.verifyEqual([a.panels.name], [c.panels.name]);
            tc.verifyEqual(a.A, c.A); tc.verifyEqual(a.N, c.N);
            tc.verifyEqual(a.optVis, c.optVis); tc.verifyEqual(a.optIR, c.optIR); tc.verifyEqual(a.heatCapacity, c.heatCapacity);
            tc.verifyNotEqual(a.mass_kg, c.mass_kg);
        end
    end
end
