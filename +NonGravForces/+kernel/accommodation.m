function [acc, info] = accommodation(model, species, T_K, options)
%ACCOMMODATION Energy accommodation of the gas-surface interaction per epoch: fixed, Pilinski (2010) or the Walker/Mehta adsorbate mixture.
%
%   acc = NonGravForces.kernel.accommodation("fixed", species, T_K, 'Alpha', 0.85)
%   acc = NonGravForces.kernel.accommodation("pilinski2010", species, T_K)
%   acc = NonGravForces.kernel.accommodation("walker2014", species, T_K, 'SurfaceMass_amu', 263.3223)
%   [acc, info] = NonGravForces.kernel.accommodation(___)
%
%   Returns what kernel.dragPanelForce needs to form the re-emitted momentum:
%   either ONE accommodation coefficient per epoch, or a surface that is a
%   mixture of two states. The kernel never branches on the model name; it reads
%   acc.coverage, which is empty for a single state.
%
%   MODELS
%     fixed         alpha = Alpha at every epoch (Alpha scalar or 1xK). The value is
%                   the caller's statement; nothing is assumed when it is omitted.
%     pilinski2010  alpha = K n_O T / (1 + K n_O T), K = 4.98e-17 m^3 K^-1, with the
%                   atomic-oxygen number density n_O and the gas temperature T of
%                   NRLMSISE-00 (Pilinski, Argrow & Palo 2010; the constant as printed
%                   in Pilinski 2011, sect. 7.3, eqs 7.1-7.2, fitted with a clean
%                   substrate of alpha 0 to satellites between 150 and 500 km, where
%                   the fitted alpha spans 0.85-1.00). Outside that range it is an
%                   extrapolation and info says so.
%     walker2014    The adsorbate model as Paul, Sheridan, Licata & Mehta (2022, eqs
%                   10-13) state it: a fraction theta of the surface carries adsorbed
%                   atomic oxygen and accommodates fully, the rest is clean substrate,
%                     theta   = K_L P_O / (1 + K_L P_O),  P_O = n_O k T  [Pa],  K_L = 1.44e6 Pa^-1
%                     alpha_a = 1
%                     alpha_s = 3 mu / (1 + mu)^2,  mu = (sum_j x_j m_j) / m_s
%                   with x_j the mole fractions of H, He, N, N2, O and O2 (normalised
%                   over those six) and m_s the mass of a surface lattice particle; the
%                   drag of the two states is mixed with weight theta. ATTRIBUTION, as
%                   that paper gives it: the coverage law and K_L to Walker, Mehta &
%                   Koller (2014, JSR 51(5)), the substrate coefficient to Mehta et al.
%                   (2022); the name of this option follows the first. K_L's unit is not
%                   printed; Pa^-1 follows from P_O formed in SI in the WVU RSM Suite
%                   that implements it (read as Torr^-1, theta would be 0.13 instead of
%                   0.95 at n_O = 1e15 m^-3, 1000 K). m_s has no universally
%                   agreed value: 263.3223 amu (Mehta et al. 2022, as cited by Paul
%                   2022), 65 amu (Pilinski 2011, Table 11.1); the owner chose the first
%                   for the precise preset (HANDOVER Part I, 2026-09-15), and it is
%                   required here, never defaulted.
%   SESAM (Pilinski et al. 2013) is not implemented: its sticking law and shape
%   factor are in no source in hand (HANDOVER O-16).
%
% INPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   model     1x1   string  -      "fixed" | "pilinski2010" | "walker2014"
%   species   8xK   double  1/m^3  Number densities in const.speciesName order (env.density)
%   T_K       1xK   double  K      Gas kinetic temperature at the satellite (density's T(2,:))
%
% NAME-VALUE ARGUMENTS
%   Name              Size      Type    Units  Description
%   ----              ----      ----    -----  -----------
%   Alpha             1x1|1xK   double  -      Coefficient in [0, 1]; required by "fixed", refused by the others
%   SurfaceMass_amu   1x1       double  amu    Surface lattice particle mass, > 0; required by "walker2014"
%
% OUTPUT ARGUMENTS
%   Name      Size  Type    Units  Description
%   ----      ----  ----    -----  -----------
%   acc       1x1   struct  -      .model; .alpha 1xK (fixed, pilinski2010; for walker2014 the
%                                  adsorbate state's, 1); .coverage 1xK (walker2014: theta; [] for a
%                                  single state); .alphaSubstrate 1xK (walker2014; [] otherwise)
%   info      1x1   struct  -      .K, .pressureO_Pa (walker2014), .meanMass_amu (walker2014),
%                                  .insideFit (pilinski2010: alpha within the fitted 0.85-1), .provenance
%
% ACCURACY
%   Exact evaluation of the published formulas. The model error is the
%   gas-surface interaction itself: fitted accommodation spans 0.82-0.93 at
%   270-500 km (March et al. 2019), i.e. 5-10 % in the drag coefficient, and the
%   substrate mass of walker2014 alone moves the effective coefficient by ~0.1
%   at GRACE-FO altitude near solar minimum (measured with the three values in
%   circulation, 2026-09-15).
%
% ERRORS
%   NonGravForces:accommodation:inputShape       species not 8xK or T_K not 1xK with the same K
%   NonGravForces:accommodation:inputAlpha       fixed without 'Alpha', Alpha outside [0, 1] or not 1x1/1xK, or Alpha given to another model
%   NonGravForces:accommodation:inputSurfaceMass walker2014 without a positive 'SurfaceMass_amu'
%   NonGravForces:accommodation:numericNaN       NaN/Inf in Alpha, species or T_K
%   NonGravForces:accommodation:inputRange       a negative density or a non-positive temperature (models that read them)
%
% REFERENCES
%   Pilinski M.D., Argrow B.M., Palo S.E. (2010) J. Spacecraft Rockets 47(6);
%   Pilinski M.D. (2011) Dynamic gas-surface interaction modeling for satellite
%   aerodynamic computations, PhD thesis, University of Colorado, sect. 7.3 and
%   Table 11.1. Walker A., Mehta P., Koller J. (2014) J. Spacecraft Rockets 51(5)
%   1544. Paul S.N., Sheridan P.L., Licata R.J., Mehta P.M. (2022) arXiv:2210.08364v2,
%   eqs 10-13. March G. et al. (2019) J. Space Weather Space Clim. 9, A36.
%
% EXAMPLE
%   [~, sp, T] = NonGravForces.env.density("msis00", ep, lla, sw);
%   acc = NonGravForces.kernel.accommodation("walker2014", sp, T(2, :), 'SurfaceMass_amu', 263.3223);
%   [acc.coverage; acc.alphaSubstrate]
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.gsiCoefficients, NonGravForces.env.density.
%
% KEYWORDS accommodation coefficient adsorption Langmuir atomic oxygen gas-surface interaction
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic); pilinski2010 and walker2014 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["fixed","pilinski2010","walker2014"])}
    species double
    T_K double
    options.Alpha double = []
    options.SurfaceMass_amu double = []
end
unit = "accommodation";
K = size(T_K, 2);
if ~(ismatrix(species) && size(species, 1) == 8 && size(species, 2) == K && isrow(T_K) || (K == 0 && isempty(species)))
    NonGravForces.fail(unit, "inputShape", "species must be 8xK and T_K 1xK with the same K -> got species %s, T_K %s", mat2str(size(species)), mat2str(size(T_K)));
end
if any(~isfinite(species(:))) || any(~isfinite(T_K(:)))
    NonGravForces.fail(unit, "numericNaN", "species or T_K contain NaN/Inf -> screen the density output");
end
c = NonGravForces.const();
acc = struct('model', model, 'alpha', [], 'coverage', [], 'alphaSubstrate', []);
info.K = K;
switch model
    case "fixed"
        a = options.Alpha;
        if isempty(a) || ~(isscalar(a) || isequal(size(a), [1 K]))
            NonGravForces.fail(unit, "inputAlpha", "fixed needs 'Alpha' as a scalar or 1xK -> pass 'Alpha', 0.85 (the value is your statement, not a default)");
        end
        if any(~isfinite(a(:)))
            NonGravForces.fail(unit, "numericNaN", "Alpha contains NaN/Inf -> check the input");
        end
        if any(a(:) < 0 | a(:) > 1)
            NonGravForces.fail(unit, "inputAlpha", "Alpha must lie in [0, 1] -> got %g", a(find(a < 0 | a > 1, 1)));
        end
        acc.alpha = a + zeros(1, K);
    case "pilinski2010"
        localNoAlpha(unit, model, options.Alpha);
        localRange(unit, species, T_K);
        KL = 4.98e-17;                                     % m^3 K^-1, Pilinski (2011) sect. 7.3
        x = KL * species(2, :) .* T_K;                      % n_O T in m^-3 K
        acc.alpha = x ./ (1 + x);
        info.insideFit = acc.alpha >= 0.85 & acc.alpha <= 1;
    case "walker2014"
        localNoAlpha(unit, model, options.Alpha);
        ms = options.SurfaceMass_amu;
        if ~(isscalar(ms) && isfinite(ms) && ms > 0)
            NonGravForces.fail(unit, "inputSurfaceMass", "walker2014 needs a positive 'SurfaceMass_amu' -> the sources disagree (263.3223 amu Mehta et al. 2022, 65 amu Pilinski 2011); state which one you use");
        end
        localRange(unit, species, T_K);
        KL = 1.44e6;                                        % Pa^-1, Walker et al. 2014 via Paul et al. 2022 eq. 11
        PO = species(2, :) .* c.kB .* T_K;                  % Pa
        acc.coverage = KL * PO ./ (1 + KL * PO);
        acc.alpha = ones(1, K);                             % the adsorbate state accommodates fully
        six = [6 1 7 3 2 4];                                % H He N N2 O O2, the species of Paul et al. 2022 eq. 13
        m = c.speciesMass(six)' / c.NA;                     % kg per molecule, 6x1
        n6 = species(six, :);
        meanMass = sum(n6 .* m, 1) ./ max(sum(n6, 1), realmin);   % mole-fraction weighted mean molecular mass
        mu = meanMass / (ms * c.amu);
        acc.alphaSubstrate = 3 * mu ./ (1 + mu).^2;
        info.pressureO_Pa = PO;
        info.meanMass_amu = meanMass / c.amu;
end
if nargout > 1
    info.provenance = NonGravForces.provenance("kernel.accommodation", 'Options', struct('model', model, 'Alpha', options.Alpha, 'SurfaceMass_amu', options.SurfaceMass_amu));
end
end

function localNoAlpha(unit, model, a)
if ~isempty(a)
    NonGravForces.fail(unit, "inputAlpha", "model %s computes its coefficient from the environment -> remove 'Alpha', or use model fixed", model);
end
end

function localRange(unit, species, T_K)
if any(species(:) < 0) || any(T_K <= 0)
    NonGravForces.fail(unit, "inputRange", "number densities must be >= 0 and T_K > 0 -> check the density output");
end
end
