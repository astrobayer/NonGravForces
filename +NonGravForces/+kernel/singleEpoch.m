function [F, tau, info] = singleEpoch(macro, gas, src, Tpanel_K, options)
%SINGLEEPOCH Compiled single-epoch kernel: gas + radiation + thermal recoil, forces and torques in one call.
%
%   [F, tau] = NonGravForces.kernel.singleEpoch(macro, gas, src, Tpanel_K, 'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300)
%   [F, tau] = NonGravForces.kernel.singleEpoch(macro, gas, src, [], 'Model', "constantCd", 'Cd', 2.2)
%   [F, tau, info] = NonGravForces.kernel.singleEpoch(___, 'CdScale', 1.03, 'CrScale', 0.97)
%
%   ONE epoch, evaluated in C (mex/ngfKernel.c through
%   +kernel/private/kernelGateway). The three surface effects are computed over
%   the same panel loop and each is reduced to a force and a torque, so the
%   columns of F and tau are drag, radiation and thermal recoil in that order
%   and their sum is the total non-gravitational force on the body:
%
%       F(:,1) = cdScale * sum_k F_k^gas        (kernel.dragPanelForce)
%       F(:,2) = crScale * sum_k F_k^rad        (kernel.radiationPanelForce)
%       F(:,3) =           sum_k F_k^thermal    (kernel.thermalRecoil)
%       tau(:,j) = sum_k (c_k - r_cm) x F_k^j   (kernel.assemble)
%
%   They stay separate because Model.accel divides each by the mass and
%   Model.budget reports them one by one; summing here would force every caller
%   to re-split what it needs.
%
%   WHY IT EXISTS. Two evaluation paths are kept side by side by design
%   (HANDOVER D.1 P-11, owner ruling 2026-09-08): Model.accel evaluates K epochs
%   at once and is the path for an arc; this unit evaluates one epoch in C and is
%   the path for an integrator's right-hand side, called once per step and never
%   amortised over K. Neither is a fidelity tier of the other, neither is a
%   fallback for the other, and neither is justified by being the faster one.
%   What the pair costs is that they must AGREE, which
%   tests/tNgfMex/singleEpochEndToEndPrecision measures over the real arc: the
%   environment composed as Model.evaluate composes it, one epoch through here,
%   and the acceleration back must be Model.accel's.
%
%   WHAT IT DOES NOT DO. No environment: the density, the composition, the
%   radiation source list, the attitude, the self-shadowing and the panel
%   temperatures are inputs, computed by env and geom in MATLAB. This unit is the
%   kernel layer and nothing else. It also does no frame work: everything in and
%   out is body frame, SI, radians.
%
%   MATLAB REMAINS NORMATIVE. The C is a port of +kernel/*.m, verified against
%   them (tests/tNgfMex.m) and against the frozen independent Python mirror
%   (mex/test/ngfKernelParity.c, run by CI). Where the two ever disagree the
%   MATLAB is right and the C is wrong.
%
% INPUT ARGUMENTS
%   Name        Size   Type        Units  Description
%   ----        ----   ----        -----  -----------
%   macro       1x1    MacroModel  -      Panel model: areas, unit normals, centroids, optics
%   gas         1x1    struct      -      .rho 1x1 kg/m^3 total mass density (used by constantCd)
%                                         .species 8x1 1/m^3 number densities in const.speciesName
%                                          order, or [] for constantCd
%                                         .T_K 1x1 K gas kinetic temperature (density's T(2))
%                                         .vrel_body 3x1 m/s satellite velocity relative to the gas,
%                                          body frame, |v| > 0
%                                         .shadow Nx1 self-shadow fraction of each plate for the gas
%                                          stream in [0,1], or [] for all 1
%   src         1x1    struct      -      .dir 3xM unit vectors satellite -> source, body frame
%                                         .P Mx1 Pa radiation pressure (the eclipse factor is already
%                                          inside the Sun's)
%                                         .band 1xM string "vis" | "ir"
%                                         .shadow NxM per plate and source in [0,1], or [] for all 1
%                                         M = 0 is legal and gives no radiation force
%   Tpanel_K    Nx1    double      K      Plate temperatures (geom.thermalState), >= 0, or [] for no
%                                         thermal recoil (Options.thermal.model == "none")
%
% NAME-VALUE ARGUMENTS
%   Name              Size  Type    Units  Description
%   ----              ----  ----    -----  -----------
%   Model             1x1   string  -      "sentman" | "maxwell" | "constantCd" (required; no default)
%   Alpha             1x1   double  -      Accommodation coefficient in [0,1]; required by sentman/maxwell
%   Tw_K              1x1   double  K      Wall temperature > 0; required by sentman/maxwell
%   Cd                1x1   double  -      Drag coefficient > 0; required by constantCd
%   SpecularFraction  1x1   double  -      In [0,1]; required by maxwell
%   CdScale           1x1   double  -      Factor on the gas force (Options.drag.cdScale), default 1
%   CrScale           1x1   double  -      Factor on the radiation force (Options.srp.crScale), default 1
%   CentreOfMass      3x1   double  m      Torque reference point, default macro.centreOfMass
%
% OUTPUT ARGUMENTS
%   Name     Size   Type    Units  Description
%   ----     ----   ----    -----  -----------
%   F        3x3    double  N      Total force, body frame; columns drag | radiation | thermal
%   tau      3x3    double  N m    Torque about CentreOfMass, same column order
%   info     1x1    struct  -      .Fk 3xNx3 per-panel forces, pages in the same order
%                                  .gamma Nx1 n.v of the gas stream, .theta_rad Nx1
%                                  .speedRatio 8x1 per species (NaN for constantCd)
%                                  .Cp .Ctau Nx1 mixture coefficients per unit rho
%                                  .absorbed_W_m2 Nx1 (before CrScale), .emitted_W_m2 Nx1
%                                  .backend "mex", .kernelVersion, .provenance
%
% ACCURACY
%   The C evaluates the same closed forms as +kernel/*.m in the same association
%   order, so the difference is rounding. Measured 2026-09-08 against the frozen
%   Python mirror (mex/test/ngfKernelParity.c, the same fixtures and bounds
%   tests/tNgfMirror.m uses on the MATLAB): gsiCoefficients 8.7e-14 relative on
%   the 617 coefficients above 1e-12 -- the SAME figure to five digits that the
%   MATLAB records against that mirror, i.e. the residual is the mirror's
%   better-conditioned leeward form (HANDOVER O-20) and not a difference between
%   the two ports; dragPanelForce 4.3e-16, radiationPanelForce 2.0e-16, absorbed
%   flux 2.7e-16, thermalRecoil 1.4e-16, assemble force and torque exactly 0,
%   all relative to the case scale. The direct comparison against the MATLAB
%   units is tests/tNgfMex.m and it bounds the difference at 1e-14.
%   The composition is exact: this unit's per-panel arrays and reductions equal
%   the four C units called separately to the last bit (measured 0).
%
% ERRORS
%   NonGravForces:singleEpoch:depMissing    the MEX is not built on this machine
%   NonGravForces:singleEpoch:inputMacro    macro is not a NonGravForces.geom.MacroModel
%   NonGravForces:singleEpoch:inputModel    'Model' missing, unknown, or missing its required options
%   NonGravForces:singleEpoch:inputShape    a field of gas or src, Tpanel_K or CentreOfMass has the wrong shape
%   NonGravForces:singleEpoch:inputSources  src lacks dir/P/band, a band is not vis|ir, or a direction is not a unit vector
%   NonGravForces:singleEpoch:inputRange    |v| = 0, rho < 0, T <= 0, alpha outside [0,1], P < 0, a shadow outside [0,1]
%   NonGravForces:singleEpoch:numericNaN    NaN or Inf in an input
%   NonGravForces:kernelGateway:*           raised by the C gateway; reaching one means a contract above let something through
%
% REFERENCES
%   Sentman L.H. (1961) LMSC-448514; Doornbos E. (2011) Thermospheric density
%   and wind determination from satellite dynamics, Ch. 3 eqs 3.55-3.57;
%   Wöske F. (2019) PhD thesis eqs 5.4 and 5.14; Montenbruck O., Gill E. (2000)
%   Satellite Orbits eq. 3.73. HANDOVER Part G LP-6b (WP-17b), D.6.
%
% EXAMPLE
%   mm  = NonGravForces.geom.graceFO("C");
%   gas = struct('rho', 5.06e-14, 'species', [2e13; 8e13; 4e11; 1e10; 2e9; 5e10; 1e12; 1e12], ...
%                'T_K', 950, 'vrel_body', [7590; 120; -45], 'shadow', []);
%   src = struct('dir', [1 0; 0 0; 0 -1], 'P', [4.56e-6; 7.4e-7], 'band', ["vis" "ir"], 'shadow', []);
%   [F, tau] = NonGravForces.kernel.singleEpoch(mm, gas, src, 260 * ones(mm.n, 1), ...
%       'Model', "sentman", 'Alpha', 0.85, 'Tw_K', 300);
%   a = sum(F, 2) / 601.214;      % 3x1 non-gravitational acceleration, body frame
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.radiationPanelForce, NonGravForces.kernel.thermalRecoil, NonGravForces.kernel.assemble, NonGravForces.buildMex, NonGravForces.Model.
%
% KEYWORDS mex compiled single epoch kernel drag radiation thermal recoil integrator
% Created 2026-09-08 by Claude Opus 5 (Anthropic).
arguments
    macro
    gas (1,1) struct
    src (1,1) struct
    Tpanel_K double
    options.Model string = string.empty
    options.Alpha double = []
    options.Tw_K double = []
    options.Cd double = []
    options.SpecularFraction double = []
    options.CdScale (1,1) double = 1
    options.CrScale (1,1) double = 1
    options.CentreOfMass double = []
end
unit = "singleEpoch";
if ~isa(macro, 'NonGravForces.geom.MacroModel')
    NonGravForces.fail(unit, "inputMacro", "macro must be a NonGravForces.geom.MacroModel -> build one with geom.graceFO or geom.boxWing");
end
gw = fullfile(fileparts(mfilename('fullpath')), 'private', "kernelGateway." + mexext);
if ~isfile(gw)
    NonGravForces.fail(unit, "depMissing", "the single-epoch kernel MEX %s is not built -> run NonGravForces.buildMex once on this machine", gw);
end
N = macro.n;

% --- gas
need = ["rho","species","T_K","vrel_body","shadow"];
missing = setdiff(need, string(fieldnames(gas)));
if ~isempty(missing)
    NonGravForces.fail(unit, "inputShape", "gas lacks %s -> see the help of this function (pass [] for species and shadow where they do not apply)", strjoin(missing, ", "));
end
model = options.Model;
if isempty(model) || ~isscalar(model) || ~any(model == ["sentman","maxwell","constantCd"])
    NonGravForces.fail(unit, "inputModel", "'Model' must be sentman, maxwell or constantCd -> state the gas-surface model explicitly");
end
if ~isscalar(gas.rho) || ~isscalar(gas.T_K) || ~isequal(size(gas.vrel_body), [3 1])
    NonGravForces.fail(unit, "inputShape", "gas.rho and gas.T_K must be 1x1 and gas.vrel_body 3x1 -> one epoch per call (use Model.series for many)");
end
species = gas.species(:);
if model == "constantCd"
    species = zeros(0, 1);
elseif ~isequal(size(species), [8 1])
    NonGravForces.fail(unit, "inputShape", "gas.species must be 8x1 for %s -> env.density returns the eight MSIS species", model);
end
shadowGas = gas.shadow(:);
if ~isempty(shadowGas) && ~isequal(size(shadowGas), [N 1])
    NonGravForces.fail(unit, "inputShape", "gas.shadow must be %dx1 or [] -> geom.selfShadow for the stream direction", N);
end

% --- radiation sources
if ~all(isfield(src, {'dir','P','band'}))
    NonGravForces.fail(unit, "inputSources", "src needs fields dir (3xM), P (Mx1), band (1xM) -> use env.radiationSources and take one epoch");
end
d = src.dir;
if isempty(d)
    d = zeros(3, 0);
end
if ~ismatrix(d) || size(d, 1) ~= 3
    NonGravForces.fail(unit, "inputSources", "src.dir must be 3xM for one epoch -> got %s", mat2str(size(d)));
end
M = size(d, 2);
P = src.P(:);
band = string(src.band(:))';
if ~isequal(size(P), [M 1]) || ~isequal(size(band), [1 M])
    NonGravForces.fail(unit, "inputSources", "src.P must be %dx1 and src.band 1x%d -> got P %s, band %s", M, M, mat2str(size(P)), mat2str(size(band)));
end
if M > 0 && ~all(band == "vis" | band == "ir")
    NonGravForces.fail(unit, "inputSources", "src.band entries must be vis or ir -> check the source list");
end
if M > 0 && any(abs(vecnorm(d, 2, 1) - 1) > 1e-9)
    NonGravForces.fail(unit, "inputSources", "src.dir must hold unit vectors (1e-9) -> normalise the directions");
end
if isfield(src, 'shadow'), shadowRad = src.shadow; else, shadowRad = []; end
if ~isempty(shadowRad) && ~isequal(size(shadowRad), [N M])
    NonGravForces.fail(unit, "inputShape", "src.shadow must be %dx%d or [] -> geom.selfShadow for each source direction", N, M);
end

% --- panel temperatures
Tp = Tpanel_K(:);
if ~isempty(Tp) && ~isequal(size(Tp), [N 1])
    NonGravForces.fail(unit, "inputShape", "Tpanel_K must be %dx1 or [] -> geom.thermalState returns one column per epoch", N);
end

% --- finiteness, then ranges. Both are checked here because the compiled path
% does not go through the MATLAB units that would otherwise check them; the C
% guards behind this are a backstop, not the contract.
probe = [gas.rho; gas.T_K; gas.vrel_body; species; shadowGas; d(:); P; shadowRad(:); Tp; options.CdScale; options.CrScale];
if any(~isfinite(probe))
    NonGravForces.fail(unit, "numericNaN", "an input contains NaN/Inf -> screen the density, sources, temperatures and velocity");
end
if ~(norm(gas.vrel_body) > 0) || gas.rho < 0 || any(P < 0) || any(Tp < 0)
    NonGravForces.fail(unit, "inputRange", "|vrel_body| > 0, rho >= 0, P >= 0 and Tpanel_K >= 0 are required -> check the inputs");
end
if model ~= "constantCd" && ~(gas.T_K > 0)
    NonGravForces.fail(unit, "inputRange", "gas.T_K must be > 0 for %s -> the speed ratio divides by it", model);
end
if any(shadowGas < 0 | shadowGas > 1) || any(shadowRad(:) < 0 | shadowRad(:) > 1)
    NonGravForces.fail(unit, "inputRange", "shadow fractions must lie in [0, 1] -> check the selfShadow output");
end

% --- the model's own required options, refused rather than defaulted (rule R5)
alpha = 0; Tw = 0; Cd = 0; sigma = 0;
switch model
    case "constantCd"
        if ~isscalar(options.Cd) || ~(options.Cd > 0)
            NonGravForces.fail(unit, "inputModel", "constantCd needs a positive scalar 'Cd' -> pass 'Cd', 2.2 for a GRACE-class body");
        end
        Cd = options.Cd;
        modelCode = 2;
    otherwise
        if ~isscalar(options.Alpha) || ~isscalar(options.Tw_K)
            NonGravForces.fail(unit, "inputModel", "%s needs scalar 'Alpha' and 'Tw_K' -> pass the accommodation coefficient (kernel.accommodation) and the wall temperature", model);
        end
        if options.Alpha < 0 || options.Alpha > 1 || ~(options.Tw_K > 0)
            NonGravForces.fail(unit, "inputRange", "'Alpha' must lie in [0, 1] and 'Tw_K' must be > 0 -> got %g and %g", options.Alpha, options.Tw_K);
        end
        alpha = options.Alpha; Tw = options.Tw_K;
        if model == "maxwell"
            if ~isscalar(options.SpecularFraction) || options.SpecularFraction < 0 || options.SpecularFraction > 1
                NonGravForces.fail(unit, "inputModel", "maxwell needs a scalar 'SpecularFraction' in [0, 1] -> pass 'SpecularFraction', 0.1");
            end
            sigma = options.SpecularFraction;
            modelCode = 1;
        else
            modelCode = 0;
        end
end
if isempty(options.CentreOfMass)
    rcm = macro.centreOfMass;
else
    rcm = options.CentreOfMass;
    if ~isequal(size(rcm), [3 1])
        NonGravForces.fail(unit, "inputShape", "CentreOfMass must be 3x1 -> got %s", mat2str(size(rcm)));
    end
end
if any(~isfinite(rcm)) || any(~isfinite(macro.C(:)))
    NonGravForces.fail(unit, "numericNaN", "CentreOfMass or the panel centroids are not all finite, so no torque can be formed -> build the MacroModel with real centroid_m");
end

% --- constants travel from MATLAB, where each has exactly one owner (rule R1)
c = NonGravForces.const();
con = struct('cLight', referenceSystem.iersConstants().speedOfLight, ...
             'sigmaSB', c.sigmaSB, 'kB', c.kB, ...
             'speciesMass', (c.speciesMass(:) / c.NA));
mac = struct('A', reshape(macro.A, 1, N), 'N', macro.N, 'C', macro.C, ...
             'optVis', macro.optVis, 'optIR', macro.optIR, 'rcm', rcm);
gasIn = struct('rho', gas.rho, 'species', species, 'T_K', gas.T_K, ...
               'vrel', gas.vrel_body, 'model', modelCode, 'alpha', alpha, ...
               'Tw_K', Tw, 'Cd', Cd, 'SpecularFraction', sigma, 'shadow', shadowGas);
radIn = struct('dir', d, 'P', P, 'band', double(band == "ir"), 'shadow', shadowRad);   % 1xM, the shape the gateway checks
if M == 0
    radIn.band = zeros(1, 0);
end
[F, tau, Fk, aux] = kernelGateway(mac, gasIn, radIn, Tp, [options.CdScale; options.CrScale], con);

info = aux;
info.Fk = Fk;
info.model = model;
info.backend = "mex";
info.kernelVersion = string(mexext) + " built from mex/ngfKernel.c";
info.provenance = NonGravForces.provenance("kernel.singleEpoch", 'Options', ...
    struct('Model', model, 'Alpha', options.Alpha, 'Tw_K', options.Tw_K, 'Cd', options.Cd, ...
           'SpecularFraction', options.SpecularFraction, 'CdScale', options.CdScale, ...
           'CrScale', options.CrScale, 'M', M, 'thermal', ~isempty(Tp)));
end
