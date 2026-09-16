function [Cp, Ctau, info] = gsiCoefficients(model, s, theta, alpha, Tw_K, Tinf_K, options)
%GSICOEFFICIENTS Normal and tangential free-molecular force coefficients of a flat plate: sentman, maxwell, constantCd.
%
%   [Cp, Ctau] = NonGravForces.kernel.gsiCoefficients("sentman", s, theta, alpha, Tw_K, Tinf_K)
%   [Cp, Ctau] = NonGravForces.kernel.gsiCoefficients("maxwell", s, theta, alpha, Tw_K, Tinf_K, 'SpecularFraction', 0.1)
%   [Cp, Ctau] = NonGravForces.kernel.gsiCoefficients("constantCd", [], theta, [], [], [], 'Cd', 2.2)
%   [Cp, Ctau, info] = NonGravForces.kernel.gsiCoefficients(___)
%
%   Force per unit plate area, in the (1/2) rho U^2 convention that the drag
%   coefficient uses:
%
%       f = (1/2) rho U^2 ( Cp (-n) + Ctau t ),      g = -n.v,  t = (g + gamma n) / l,
%
%   where n is the outward plate normal, v the unit vector of the satellite's
%   velocity relative to the gas, g = -v the direction the gas streams in,
%   theta the angle between n and v (0: head-on, > pi/2: leeward),
%   gamma = cos theta = n.v and l = sin theta; t is the tangential part of the
%   gas stream. Cp > 0 pushes the plate inward, Ctau > 0 drags it along the
%   tangential stream. The drag (along g) and lift (in the n-g plane) of the
%   plate follow as CD = Cp gamma + Ctau l and CL = Cp l - Ctau gamma; both are
%   returned in info. dragPanelForce does the body-frame assembly.
%
%   MODELS
%     sentman     Free-molecular flow of a drifting Maxwellian on a plate with
%                 diffuse re-emission (Sentman 1961). With P = exp(-gamma^2 s^2)/s,
%                 Z = 1 + erf(gamma s), the momentum integrals give
%                   Cp_inc   = gamma P/sqrt(pi) + (gamma^2 + 1/(2 s^2)) Z
%                   Ctau_inc = l (P/sqrt(pi) + gamma Z)
%                 and the diffusely re-emitted molecules add, along the normal,
%                   Cp_re    = (1/2) (v_re/v_inc) (gamma sqrt(pi) Z + P),
%                   v_re/v_inc = sqrt( (1/2) [1 + alpha (2 Tw/(Tinf s^2) - 1)] )
%                 (Koppenwallner's energy accommodation: the re-emission
%                 temperature is T_i - alpha (T_i - Tw) with the energy-equivalent
%                 incident temperature T_i = m U^2/(4 k)). Cp gamma + Ctau l is
%                 Doornbos' (2011, eqs 3.55-3.56) C_D and Cp l - Ctau gamma his
%                 C_L. The theory document (§2.3) labelled that drag/lift pair
%                 Cp/Ctau until 0.10.0 and now derives both pairs; the
%                 hyperthermal limit (2 gamma^2 normal, 2 gamma l tangential,
%                 2 gamma drag) shows which is which. Species enter through s: one call per species
%                 and a density-weighted sum (dragPanelForce does this).
%     maxwell     A fraction SpecularFraction of the molecules reflects
%                 specularly (incident normal momentum doubled, no shear, no
%                 re-emission), the rest as sentman.
%     constantCd  Cp = Cd max(gamma, 0) gamma, Ctau = Cd max(gamma, 0) l: the
%                 force of every windward plate points along the gas stream (no
%                 lift), and the plate sum reproduces the textbook
%                 F = -(1/2) rho A_proj Cd |v| v with the projected area as
%                 reference (Wöske 2019 eq. 5.20). s, alpha, Tw, Tinf are ignored
%                 and may be empty.
%     cll, schamberg  refused (unsupportedModel): no closed flat-plate form
%                 and no oracle in hand (HANDOVER O-17).
%
% INPUT ARGUMENTS
%   Name     Size     Type    Units  Description
%   ----     ----     ----    -----  -----------
%   model    1x1      string  -      "sentman" | "maxwell" | "constantCd" | "cll" | "schamberg"
%   s        any      double  -      Speed ratio U / sqrt(2 k T / m) > 0 (ignored by constantCd)
%   theta    any      double  rad    Angle between the plate normal and the relative velocity v, in [0, pi]
%   alpha    any      double  -      Energy accommodation coefficient in [0, 1]
%   Tw_K     any      double  K      Wall temperature > 0
%   Tinf_K   any      double  K      Gas kinetic temperature > 0
%   All arrays must be compatible for implicit expansion (e.g. theta NxK, s 1xK).
%
% NAME-VALUE ARGUMENTS
%   Name              Size  Type    Units  Description
%   ----              ----  ----    -----  -----------
%   Cd                1x1   double  -      Drag coefficient for constantCd, > 0 (required by that model)
%   SpecularFraction  1x1   double  -      Specular fraction for maxwell, in [0, 1] (required by that model)
%   Derivatives       1x1   logical -      Also return the analytic derivatives of the coefficients with
%                                          respect to alpha and to the speed ratio s, in info (default false).
%                                          They live here, beside the closed form they differentiate, so the
%                                          two cannot drift apart (rule R1); kernel.partials assembles them.
%   Provenance        1x1   logical -      Build info.provenance when info is requested (default true); a
%                                          kernel that calls this per species passes false, because the record
%                                          costs more than the coefficients (HANDOVER LP-C, finding F12).
%
% OUTPUT ARGUMENTS
%   Name     Size     Type    Units  Description
%   ----     ----     ----    -----  -----------
%   Cp       expanded double  -      Normal coefficient (>= 0)
%   Ctau     expanded double  -      Tangential coefficient (>= 0)
%   info     1x1      struct  -      .model, .CD, .CL (same size as Cp), .CpIncident and
%                                    .CpReemitted (Cp = CpIncident + CpReemitted: the incident normal
%                                    momentum, including a specular fraction, and the diffusely
%                                    re-emitted one, which is proportional to the incident NUMBER
%                                    flux; geom.flowShadow scales the two differently), .vRatio
%                                    (re-emission speed ratio; empty for constantCd), .provenance.
%                                    With 'Derivatives': .dCp_dalpha, .dCtau_dalpha (identically zero:
%                                    alpha enters only the re-emitted normal momentum),
%                                    .dCp_ds, .dCtau_ds, .dCpIncident_ds, .dCpReemitted_ds, all the size of Cp.
%
% ACCURACY
%   Closed forms evaluated once per array (erfc and exp): rounding only, 1e-15
%   relative (Z = erfc(-gamma s) keeps leeward plates positive where
%   1 + erf(gamma s) cancels to 0 or below; measured 2026-09-07). Agreement with the numerical kinetic-theory integration in
%   tools/pyxcheck/sentman_check.py is measured by the reference test (the
%   oracle's own Simpson error is below 1e-10). The model error is the gas-
%   surface interaction itself: Sentman with alpha = 0.85 agrees with GRACE/
%   CHAMP-derived coefficients to 5-10 % (March et al. 2019); the tangential
%   term of a leeward plate (gamma < 0) is exponentially small, not zero.
%
% ERRORS
%   NonGravForces:gsiCoefficients:unsupportedModel   cll or schamberg (O-17)
%   NonGravForces:gsiCoefficients:inputShape         arrays not compatible for expansion
%   NonGravForces:gsiCoefficients:inputRange         s <= 0, theta outside [0, pi], alpha outside [0, 1], Tw or Tinf <= 0
%   NonGravForces:gsiCoefficients:inputCd            constantCd without a positive Cd
%   NonGravForces:gsiCoefficients:inputSpecular      maxwell without a SpecularFraction in [0, 1]
%   NonGravForces:gsiCoefficients:numericNaN         NaN or Inf in an input
%
% REFERENCES
%   Sentman L.H. (1961) Free molecule flow theory and its application to the
%   determination of aerodynamic forces, LMSC-448514. Doornbos E. (2011)
%   Thermospheric density and wind determination from satellite dynamics,
%   Ch. 3 (eqs 3.55-3.57). Koppenwallner G. (2008) Comment on special section:
%   new perspectives on the satellite drag environments, JSR 45(6). Wöske F.
%   (2019) PhD thesis, eq. 5.20. March et al. (2019) J. Space Weather Space Clim. 9, A36.
%
% EXAMPLE
%   [Cp, Ctau, info] = NonGravForces.kernel.gsiCoefficients("sentman", 8, 0, 1, 300, 1000);
%   [Cp, info.CD]   % head-on plate: 2.137 (2 + 1/64 incident, 0.121 re-emission)
%
% See also NonGravForces.kernel.dragPanelForce, NonGravForces.kernel.accommodation.
%
% KEYWORDS gas-surface interaction Sentman Maxwell drag coefficient free molecular
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
arguments
    model (1,1) string {mustBeMember(model, ["sentman","maxwell","constantCd","cll","schamberg"])}
    s double
    theta double
    alpha double
    Tw_K double
    Tinf_K double
    options.Cd double = []
    options.SpecularFraction double = []
    options.Derivatives (1,1) logical = false
    options.Provenance (1,1) logical = true
end
unit = "gsiCoefficients";
if any(model == ["cll","schamberg"])
    NonGravForces.fail(unit, "unsupportedModel", ...
        "model %s has no closed flat-plate form and no oracle in this toolbox (HANDOVER O-17) -> use sentman, maxwell or constantCd", model);
end
if any(~isfinite(theta(:)))
    NonGravForces.fail(unit, "numericNaN", "theta contains NaN/Inf -> screen the geometry");
end
if any(theta(:) < 0 | theta(:) > pi)
    NonGravForces.fail(unit, "inputRange", "theta must lie in [0, pi] rad -> take atan2(|n x u|, n.u), not a signed angle");
end
g = cos(theta);
l = sin(theta);
info.model = model;
if model == "constantCd"
    if isempty(options.Cd) || ~isscalar(options.Cd) || ~(options.Cd > 0) || ~isfinite(options.Cd)
        NonGravForces.fail(unit, "inputCd", "constantCd needs a positive scalar 'Cd' -> pass 'Cd', 2.2 for a GRACE-class body");
    end
    Cp = options.Cd .* max(g, 0) .* g;
    Ctau = options.Cd .* max(g, 0) .* l;
    info.vRatio = [];
    info.CpIncident = Cp;
    info.CpReemitted = zeros(size(Cp));
    if options.Derivatives
        % constantCd has neither an accommodation coefficient nor a speed ratio.
        info.dCp_dalpha = zeros(size(Cp)); info.dCtau_dalpha = zeros(size(Cp));
        info.dCp_ds = zeros(size(Cp));     info.dCtau_ds = zeros(size(Cp));
        info.dCpIncident_ds = zeros(size(Cp)); info.dCpReemitted_ds = zeros(size(Cp));
    end
else
    for v = {s, alpha, Tw_K, Tinf_K}
        if any(~isfinite(v{1}(:)))
            NonGravForces.fail(unit, "numericNaN", "s, alpha, Tw_K or Tinf_K contain NaN/Inf -> check the density and temperature inputs");
        end
    end
    try
        probe = s + theta + alpha + Tw_K + Tinf_K; %#ok<NASGU>
    catch
        NonGravForces.fail(unit, "inputShape", "s, theta, alpha, Tw_K, Tinf_K must be compatible for implicit expansion -> got %s, %s, %s, %s, %s", ...
            mat2str(size(s)), mat2str(size(theta)), mat2str(size(alpha)), mat2str(size(Tw_K)), mat2str(size(Tinf_K)));
    end
    if any(s(:) <= 0) || any(alpha(:) < 0 | alpha(:) > 1) || any(Tw_K(:) <= 0) || any(Tinf_K(:) <= 0)
        NonGravForces.fail(unit, "inputRange", "s > 0, alpha in [0, 1], Tw_K > 0, Tinf_K > 0 are required -> check the inputs");
    end
    P = exp(-(g .* s).^2) ./ s;
    Z = erfc(-g .* s);                              % = 1 + erf(g s) without cancellation on leeward plates (g s << 0)
    CpInc = g .* P / sqrt(pi) + (g.^2 + 1 ./ (2 * s.^2)) .* Z;
    CtauInc = l .* (P / sqrt(pi) + g .* Z);
    vRatio = sqrt(0.5 * ((1 - alpha) + alpha .* 2 .* Tw_K ./ (Tinf_K .* s.^2)));   % = sqrt((1/2)[1 + alpha(2Tw/(Tinf s^2) - 1)]) without the 1 - 1 cancellation at alpha = 1
    CpRe = 0.5 * vRatio .* (g * sqrt(pi) .* Z + P);
    switch model
        case "sentman"
            Cp = CpInc + CpRe;
            Ctau = CtauInc + zeros(size(Cp));
            info.CpIncident = CpInc + zeros(size(Cp));
            info.CpReemitted = CpRe + zeros(size(Cp));
        case "maxwell"
            sig = options.SpecularFraction;
            if isempty(sig) || ~isscalar(sig) || ~isfinite(sig) || sig < 0 || sig > 1
                NonGravForces.fail(unit, "inputSpecular", "maxwell needs a scalar 'SpecularFraction' in [0, 1] -> pass 'SpecularFraction', 0.1");
            end
            % the diffuse part (1 - sig) carries incident and re-emitted momentum, the specular part doubles the
            % incident normal momentum: incident (1 + sig) CpInc, re-emitted (1 - sig) CpRe
            grid = zeros(size(CpInc)) + zeros(size(CpRe));      % the broadcast shape of the two parts
            info.CpIncident = (1 + sig) * CpInc + grid;
            info.CpReemitted = (1 - sig) * CpRe + grid;
            Cp = info.CpIncident + info.CpReemitted;
            Ctau = (1 - sig) * CtauInc + zeros(size(Cp));
    end
    info.vRatio = vRatio + zeros(size(Cp));
    if options.Derivatives
        % Differentiated from the lines above, in the same order, so a change to one is visible
        % against the other. E = exp(-(g s)^2) is the exponential inside P and dZ/ds.
        E = exp(-(g .* s).^2);
        dP_ds = -E .* (2 * (g .* s).^2 + 1) ./ s.^2;
        dZ_ds = 2 * g .* E / sqrt(pi);
        dCpInc_ds = g / sqrt(pi) .* dP_ds - Z ./ s.^3 + (g.^2 + 1 ./ (2 * s.^2)) .* dZ_ds;
        dCtauInc_ds = l .* (dP_ds / sqrt(pi) + g .* dZ_ds);
        b = 2 * Tw_K ./ (Tinf_K .* s.^2);                 % the alpha-weighted term inside vRatio^2
        vSafe = max(vRatio, realmin);                     % vRatio = 0 only at alpha = 1, Tw = 0
        dv_dalpha = 0.25 * (b - 1) ./ vSafe;
        dv_ds = -alpha .* b ./ (2 * s .* vSafe);
        bracket = g * sqrt(pi) .* Z + P;
        dCpRe_dalpha = 0.5 * dv_dalpha .* bracket;
        dCpRe_ds = 0.5 * dv_ds .* bracket + 0.5 * vRatio .* (g * sqrt(pi) .* dZ_ds + dP_ds);
        switch model
            case "sentman"
                info.dCp_ds = dCpInc_ds + dCpRe_ds;
                info.dCtau_ds = dCtauInc_ds + zeros(size(Cp));
                info.dCp_dalpha = dCpRe_dalpha + zeros(size(Cp));
                info.dCpIncident_ds = dCpInc_ds + zeros(size(Cp));
                info.dCpReemitted_ds = dCpRe_ds + zeros(size(Cp));
            case "maxwell"
                info.dCpIncident_ds = (1 + sig) * dCpInc_ds + zeros(size(Cp));
                info.dCpReemitted_ds = (1 - sig) * dCpRe_ds + zeros(size(Cp));
                info.dCp_ds = info.dCpIncident_ds + info.dCpReemitted_ds;
                info.dCtau_ds = (1 - sig) * dCtauInc_ds + zeros(size(Cp));
                info.dCp_dalpha = (1 - sig) * dCpRe_dalpha + zeros(size(Cp));
        end
        % alpha reaches the tangential coefficient through nothing: the re-emitted molecules leave
        % diffusely and carry no net tangential momentum.
        info.dCtau_dalpha = zeros(size(Cp));
    end
end
info.CD = Cp .* g + Ctau .* l;
info.CL = Cp .* l - Ctau .* g;
if nargout > 2 && options.Provenance
    info.provenance = NonGravForces.provenance("kernel.gsiCoefficients", 'Options', struct('model', model, 'Cd', options.Cd, 'SpecularFraction', options.SpecularFraction));
end
end
