function [S, info] = flowShadow(macro, vrel_body, T_K, options)
%FLOWSHADOW Self-shadowing of a free-molecular gas stream: flux-weighted occlusion ratios per plate, epoch and species.
%
%   S = NonGravForces.geom.flowShadow(macro, vrel_body, T_K, 'Method', "lut", 'Lut', lut)
%   S = NonGravForces.geom.flowShadow(macro, vrel_body, T_K)                 % exact ray cast (slow)
%   S = NonGravForces.geom.flowShadow(macro, vrel_body, [], 'Lut', lut, 'Method', "lut")   % hyperthermal only
%   [S, info] = NonGravForces.geom.flowShadow(___, 'Integrals', true)
%
%   WHY THIS IS NOT geom.selfShadow ALONG v_rel. Gas molecules do not arrive as a
%   parallel beam. Relative to the body the gas is a Maxwellian drifting with
%   -v_rel, so molecules come from a spread of directions about +v_rel of angular
%   width ~1/s (s = |v_rel| / sqrt(2kT/m): ~8 for O, ~3.5 for He, ~2 for H at 500
%   km). Sentman's coefficients are integrals over that spread and over each
%   plate's front hemisphere, which is why an edge-on or even leeward plate carries
%   force. An obstruction must therefore act INSIDE those integrals, weighted by
%   the flux each direction delivers (HANDOVER LP-C, WP-C1; design approved by
%   the owner 2026-09-15).
%
%   THE RATIOS. For plate k with outward normal n, tangential stream direction t
%   (kernel.dragPanelForce's), unobstructed fraction o(w) toward direction w
%   (geom.selfShadow), and the incident molecules of one species, three
%   directional integrals over the unit sphere,
%       I_N = int (n.w)+          E(x) J3(x) dW          number flux
%       I_P = int (n.w)+^2        E(x) J4(x) dW          normal momentum
%       I_T = int (n.w)+ (-w.t)   E(x) J4(x) dW          tangential momentum
%   with x = s cos(chi), chi the angle between w and v_rel, E = exp(x^2 - s^2)
%   and J_m(x) = int_0^inf u^m exp(-(u - x)^2) du in closed form, are formed
%   with and without o(w). S.number, S.normal and S.tangential are the ratios
%   (with o) / (without o). kernel.dragPanelForce scales Sentman's re-emitted
%   normal coefficient by S.number (re-emission is proportional to the number
%   flux), the incident normal coefficient by S.normal and the tangential one by
%   S.tangential. Without o the integrals ARE Sentman's incident coefficients:
%   CN = I_N / (s pi^1.5), Cp_inc = 2 I_P / (s^2 pi^1.5), Ctau_inc = 2 I_T / (s^2 pi^1.5),
%   returned with 'Integrals' and asserted against the closed forms by the tests.
%
%   RANGES. S.number and S.normal lie in [0, 1]: their weights are non-negative.
%   S.tangential does NOT: the tangential weight (-w.t) changes sign across the
%   hemisphere, so an obstruction that removes counter-flowing molecules RAISES the
%   net tangential momentum (ratio above 1), and near head-on incidence it can
%   reverse it (negative). Both are physics, not rounding.
%
%   LIMITS THAT PIN IT. A plate with nothing in front of it gets exactly 1 in
%   all three ratios, bitwise (the same sums with o = 1). As s -> inf the
%   distribution collapses onto +v_rel and every ratio tends to o(v_rel), the
%   parallel-beam occlusion, which is S.beam and what constantCd uses.
%
%   QUADRATURE. Directions are rings of constant chi about v_rel. On each ring
%   the plate factor is A + b cos(psi - psi0), so the front-facing arc is known
%   in closed form and integrated by Gauss-Legendre on that arc alone: the kink
%   at the plate's horizon never sits inside a panel. The rings are
%   Gauss-Legendre in chi on segments split at the plate's horizon crossings
%   |pi/2 - theta_n| and pi - |pi/2 - theta_n| and at 3/s of the lightest, the
%   middle and the heaviest species, so each species' peak is resolved.
%
% INPUT ARGUMENTS
%   Name        Size  Type                           Units  Description
%   ----        ----  ----                           -----  -----------
%   macro       1x1   NonGravForces.geom.MacroModel  -      With vertices on every plate
%   vrel_body   3xK   double                         m/s    Satellite velocity relative to the gas, body frame, |v| > 0
%   T_K         1xK   double                         K      Gas kinetic temperature at the satellite, > 0; [] for the
%                                                           hyperthermal limit only (S.beam)
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type     Units  Description
%   ----        ----  ----     -----  -----------
%   Method      1x1   string   -      "raycast" (default, exact to the sampling) | "lut"
%   Samples     1x1   double   -      ray-cast samples per plate edge (default 40)
%   Lut         1x1   struct   -      a geom.selfShadow table for this macro model (required for "lut")
%   Quadrature  1x2   double   -      [nodes per chi segment, nodes per arc], default [12 16]
%   Integrals   1x1   logical  -      also return the unobstructed coefficient integrals in info (default false)
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   S        1x1   struct  -      .beam NxK (occlusion along v_rel), .number .normal .tangential NxKx8
%                                 (species in const.speciesName order; [] without T_K), .speedRatio
%                                 8xK, .kind "flow"
%   info     1x1   struct  -      .method, .quadrature, .directionsPerEpoch, .integrals (struct CN
%                                 CpInc CtauInc NxKx8, with 'Integrals'), .provenance
%
% ACCURACY
%   The speed moments are closed forms (erfc and erfcx). The directional
%   quadrature reproduces Sentman's incident closed forms to 1.6e-6 in
%   coefficient units (a ram plate's coefficient is ~2) over s in [1, 15] and
%   incidence 0-180 deg at the default [12 16], measured 2026-09-15 and recorded
%   per run by tests/tNgfGeometry.m; [8 16] gives 3e-4. Far on the lee side the
%   moments cancel (terms of order e^(-s^2) x^2 against a result of order
%   e^(-s^2)/x^4), which costs relative digits on values below 1e-6 and nothing
%   that reaches a force. A ratio inherits the occlusion's own error (geom.selfShadow: the
%   ray cast's ~2/Samples, the table's 90/n deg smearing) and none from the
%   quadrature where the occlusion is uniform, because the same sums appear
%   above and below. What is NOT modelled: molecules re-emitted from one plate
%   that strike another (multiple reflections). Against the independent
%   test-particle Monte Carlo tools/pyxcheck/tpmc_check.py the incident fluxes on
%   a plate screened by a wall agree within its standard errors (recorded per run).
%
% ERRORS
%   NonGravForces:flowShadow:inputMacro   not a MacroModel, or without vertices on every plate
%   NonGravForces:flowShadow:inputShape   vrel_body not 3xK or not finite, T_K not 1xK or [], Quadrature not positive integers
%   NonGravForces:flowShadow:inputRange   |vrel| = 0 or T_K <= 0
%   NonGravForces:flowShadow:inputLut     Method "lut" without a 'Lut'
%   NonGravForces:selfShadow:*            propagated from the occlusion look-up
%
% REFERENCES
%   Sentman L.H. (1961) LMSC-448514 (the incident momentum integrals); Bird G.A.
%   (1994) Molecular Gas Dynamics and the Direct Simulation of Gas Flows, sect. 4.2
%   (the flux of a drifting Maxwellian across a surface element); Golub G.H.,
%   Welsch J.H. (1969) Math. Comp. 23:221 (Gauss-Legendre nodes); HANDOVER LP-C WP-C1.
%
% EXAMPLE
%   mm = NonGravForces.geom.graceFO("C");
%   lut = NonGravForces.geom.selfShadow(mm, [], 'Method', "lut");
%   S = NonGravForces.geom.flowShadow(mm, [-7590; 150; 40], 950, 'Method', "lut", 'Lut', lut);
%   squeeze(S.tangential(:, 1, 2))'      % atomic oxygen: 1 on every plate nothing screens
%
% See also NonGravForces.geom.selfShadow, NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.gsiCoefficients.
%
% KEYWORDS self-shadowing free molecular flow Maxwellian occlusion drag
% Created 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    macro
    vrel_body double
    T_K double
    options.Method (1,1) string {mustBeMember(options.Method, ["raycast","lut"])} = "raycast"
    options.Samples (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.Samples, 2)} = 40
    options.Lut struct = struct([])
    options.Quadrature double = [12 16]
    options.Integrals (1,1) logical = false
end
unit = "flowShadow";
if ~isa(macro, 'NonGravForces.geom.MacroModel') || ~macro.hasVertices
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel with vertices_m on every plate -> geom.graceFO and geom.boxWing (bus) provide them");
end
if ~(ismatrix(vrel_body) && size(vrel_body, 1) == 3) || any(~isfinite(vrel_body(:)))
    NonGravForces.fail(unit, "inputShape", "vrel_body must be a finite 3xK array -> one relative velocity per column");
end
K = size(vrel_body, 2);
withSpecies = ~isequal(size(T_K), [0 0]);          % [] means "hyperthermal only"; a 1x0 array is zero epochs WITH a temperature
if withSpecies && ~(isequal(size(T_K), [1 K]) && all(isfinite(T_K)))
    NonGravForces.fail(unit, "inputShape", "T_K must be a finite 1x%d array, or [] for the hyperthermal limit -> got %s", K, mat2str(size(T_K)));
end
q = options.Quadrature;
if ~(isequal(size(q), [1 2]) && all(q == round(q)) && all(q >= 1))
    NonGravForces.fail(unit, "inputShape", "'Quadrature' must be two positive integers [chi nodes per segment, arc nodes] -> got %s", mat2str(q));
end
U = vecnorm(vrel_body);
if any(U <= 0) || (withSpecies && any(T_K <= 0))
    NonGravForces.fail(unit, "inputRange", "|vrel_body| must be > 0 and T_K > 0 -> check the inputs");
end
occ = {'Method', options.Method, 'Samples', options.Samples};
if options.Method == "lut"
    if isempty(options.Lut)
        NonGravForces.fail(unit, "inputLut", "Method lut needs the table -> pass 'Lut', geom.selfShadow(macro, [], 'Method', ""lut"")");
    end
    occ = [occ, {'Lut', options.Lut}];
end
N = macro.n;
v = vrel_body ./ U;
if K == 0
    beam = zeros(N, 0);                            % selfShadow reads an empty direction list as "return the table"
else
    beam = NonGravForces.geom.selfShadow(macro, v, occ{:});
end
S = struct('kind', "flow", 'beam', beam, 'number', [], 'normal', [], 'tangential', [], 'speedRatio', []);
integrals = struct('CN', [], 'CpInc', [], 'CtauInc', []);
if withSpecies
    c = NonGravForces.const();
    mass = c.speciesMass(:) / c.NA;                              % kg per molecule, 8x1
    s = U ./ sqrt(2 * c.kB * T_K ./ mass);                        % 8xK
    S.speedRatio = s;
    [xg, wg] = localGauss(q(1)); [xa, wa] = localGauss(q(2));
    % flow frame per epoch: e1, e2 orthonormal and perpendicular to v
    useY = abs(v(1, :)) > 0.9;
    a = zeros(3, K); a(1, ~useY) = 1; a(2, useY) = 1;
    e1 = a - sum(a .* v, 1) .* v; e1 = e1 ./ vecnorm(e1);
    e2 = cross(v, e1, 1);
    [~, iMax] = max(mass); [~, iMin] = min(mass); iMid = 2;       % Ar, H and atomic O set the peak splits
    fN = ones(N, K, 8); fP = ones(N, K, 8); fT = ones(N, K, 8);
    if options.Integrals
        integrals.CN = zeros(N, K, 8); integrals.CpInc = zeros(N, K, 8); integrals.CtauInc = zeros(N, K, 8);
    end
    block = 128;
    for k = 1:N
        % A table row that is 1 everywhere reads back exactly 1 in every direction (difference-form bilinear
        % interpolation), so every ratio of that plate is exactly 1 and the quadrature has nothing to change.
        % The integrals are still formed when they are asked for.
        if options.Method == "lut" && all(options.Lut.f(k, :) == 1) && ~options.Integrals
            continue
        end
        nk = macro.N(:, k);
        for b0 = 1:block:K
            idx = b0:min(b0 + block - 1, K); Kb = numel(idx);
            nz = (nk' * v(:, idx))'; nx = (nk' * e1(:, idx))'; ny = (nk' * e2(:, idx))';   % Kbx1
            B = hypot(nx, ny); psi0 = atan2(ny, nx); thetaN = atan2(B, nz);
            chi1 = abs(pi / 2 - thetaN); chi2 = pi - chi1;
            sb = s(:, idx);
            [bp, perm] = sort(min(pi, max(0, [zeros(Kb, 1), 3 ./ sb(iMax, :)', 3 ./ sb(iMid, :)', 3 ./ sb(iMin, :)', chi1, chi2, pi * ones(Kb, 1)])), 2);
            isHorizon = ismember(perm, [5 6]);
            lo = bp(:, 1:6); hi = bp(:, 2:7);                                             % Kbx6 segments
            % The ring integral behaves like |chi - chi_horizon|^(3/2) at a horizon crossing, which caps
            % Gauss-Legendre at algebraic convergence (measured ~n^-5). A horizon END of a segment is therefore
            % reached through a quadratic map (cubic when both ends are horizons), under which that term is a
            % polynomial; the other ends stay linear so the species peaks keep their nodes.
            tau = reshape((xg + 1) / 2, 1, 1, []);
            wT = reshape(wg / 2, 1, 1, []);
            hl = isHorizon(:, 1:6); hh = isHorizon(:, 2:7);
            map = tau + zeros(Kb, 6); dmap = ones(Kb, 6) + 0 * tau;
            m1 = hl & ~hh; m2 = ~hl & hh; m3 = hl & hh;
            map = localSet(map, m1, tau.^2);                 dmap = localSet(dmap, m1, 2 * tau);
            map = localSet(map, m2, 1 - (1 - tau).^2);       dmap = localSet(dmap, m2, 2 * (1 - tau));
            map = localSet(map, m3, 3 * tau.^2 - 2 * tau.^3); dmap = localSet(dmap, m3, 6 * tau .* (1 - tau));
            chi = reshape(lo + (hi - lo) .* map, Kb, []);                                % KbxR
            wchi = reshape((hi - lo) .* dmap .* wT, Kb, []) .* sin(chi);
            A = nz .* cos(chi); bb = B .* sin(chi);
            delta = atan2(sqrt(max(bb.^2 - A.^2, 0)), -A);                                  % half-width of the lit arc
            delta(A >= bb) = pi; delta(A <= -bb) = 0;
            u = delta .* reshape(xa, 1, 1, []);                                             % psi - psi0, KbxRxnA
            wpsi = delta .* reshape(wa, 1, 1, []);
            cu = cos(u);
            cn = max(A + bb .* cu, 0);                                                       % (n.w)+
            % directions w = sin(chi) cos(psi) e1 + sin(chi) sin(psi) e2 + cos(chi) v
            psi = psi0 + u;
            sc = sin(chi); cc = cos(chi);
            cp = sc .* cos(psi); sp = sc .* sin(psi);
            W = zeros(3, numel(cp));
            for i = 1:3
                W(i, :) = reshape(e1(i, idx)' .* cp + e2(i, idx)' .* sp + v(i, idx)' .* cc, 1, []);
            end
            o = reshape(NonGravForces.geom.selfShadow(macro, W, occ{:}, 'Plate', k), size(cp));
            l = B;
            if all(l > 0)
                tang = cc .* l - (nz ./ l) .* bb .* cu;                                     % -w.t
            else
                ll = max(l, realmin); tang = cc .* ll - (nz ./ ll) .* bb .* cu;
            end
            wc = wpsi .* cn; wco = (wpsi .* o) .* cn;
            % Accumulated node by node, NOT with sum(., 3): MATLAB's sum along a dimension changes its internal
            % order with the size of the other dimensions, and a batch of epochs must equal single epochs bitwise
            % (measured 2026-09-15: sum() differed by 1e-16 between K = 5 and K = 1).
            [SN0, SNO, SP0, SPO, ST0, STO] = deal(zeros(Kb, size(cn, 2)));
            for ia = 1:size(cn, 3)
                SN0 = SN0 + wc(:, :, ia);                 SNO = SNO + wco(:, :, ia);
                SP0 = SP0 + wc(:, :, ia) .* cn(:, :, ia); SPO = SPO + wco(:, :, ia) .* cn(:, :, ia);
                ST0 = ST0 + wc(:, :, ia) .* tang(:, :, ia); STO = STO + wco(:, :, ia) .* tang(:, :, ia);
            end
            for j = 1:8
                sj = sb(j, :)';
                x = sj .* cc;
                es = exp(-sj.^2);
                % E G0 = (sqrt(pi)/2) exp(x^2 - s^2) erfc(-x), written so that neither factor over- or underflows:
                % exp(x^2 - s^2) <= 1 since |x| <= s, and for x < 0 the product goes through erfcx.
                s2 = sj.^2 .* ones(size(x));
                esx = es .* ones(size(x));
                EG0 = zeros(size(x));
                pos = x >= 0;
                EG0(pos) = (sqrt(pi) / 2) * exp(x(pos).^2 - s2(pos)) .* erfc(-x(pos));
                EG0(~pos) = (sqrt(pi) / 2) * esx(~pos) .* erfcx(-x(~pos));
                EJ3 = es .* (x.^2 + 1) / 2 + (x.^3 + 1.5 * x) .* EG0;
                EJ4 = es .* (0.5 * x.^3 + 1.25 * x) + (x.^4 + 3 * x.^2 + 0.75) .* EG0;
                g3 = wchi .* EJ3; g4 = wchi .* EJ4;
                [IN0, INO, IP0, IPO, IT0, ITO] = deal(zeros(Kb, 1));
                for ir = 1:size(g3, 2)                    % ring by ring, for the same reason as above
                    IN0 = IN0 + g3(:, ir) .* SN0(:, ir); INO = INO + g3(:, ir) .* SNO(:, ir);
                    IP0 = IP0 + g4(:, ir) .* SP0(:, ir); IPO = IPO + g4(:, ir) .* SPO(:, ir);
                    IT0 = IT0 + g4(:, ir) .* ST0(:, ir); ITO = ITO + g4(:, ir) .* STO(:, ir);
                end
                rN = INO ./ IN0; rN(~(IN0 > 0)) = 1;
                rP = IPO ./ IP0; rP(~(IP0 > 0)) = 1;
                rT = ITO ./ IT0; bad = ~(l > 0) | ~(abs(IT0) > 0); rT(bad) = rN(bad);
                fN(k, idx, j) = rN; fP(k, idx, j) = rP; fT(k, idx, j) = rT;
                if options.Integrals
                    integrals.CN(k, idx, j) = IN0 ./ (sj * pi^1.5);
                    integrals.CpInc(k, idx, j) = 2 * IP0 ./ (sj.^2 * pi^1.5);
                    integrals.CtauInc(k, idx, j) = 2 * IT0 ./ (sj.^2 * pi^1.5);
                end
            end
        end
    end
    S.number = fN; S.normal = fP; S.tangential = fT;
end
if nargout > 1
    info.method = options.Method;
    info.quadrature = q;
    info.directionsPerEpoch = N * 6 * q(1) * q(2);
    info.integrals = integrals;
    info.provenance = NonGravForces.provenance("geom.flowShadow", 'Options', struct('method', options.Method, 'quadrature', q));
end
end

function X = localSet(X, mask, value)
% X is Kb x 6 x nG, mask Kb x 6, value 1x1xnG: overwrite the masked (epoch, segment) rows with value.
nG = size(X, 3);
M = repmat(mask, 1, 1, nG);
V = repmat(value, size(mask, 1), size(mask, 2), 1);
X(M) = V(M);
end

function [x, w] = localGauss(n)
% Gauss-Legendre nodes and weights on [-1, 1] by Golub-Welsch (symmetric tridiagonal eigenproblem).
if n == 1
    x = 0; w = 2; return
end
b = (1:n-1) ./ sqrt(4 * (1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
x = x';
w = 2 * V(1, i).^2;
end
