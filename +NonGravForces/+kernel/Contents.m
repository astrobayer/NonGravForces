% KERNEL  NonGravForces.kernel - L4 layer.
% Version 0.8.6  08-Sep-2026
%
%   gsiCoefficients      - Normal and tangential plate coefficients: sentman maxwell constantCd (cll/schamberg refuse)
%   accommodation        - Energy accommodation coefficient: fixed (sesam/langmuir refuse until O-16)
%   dragPanelForce       - Per-panel gas-momentum force in body frame
%   radiationPanelForce  - Per-panel radiation force for a source list in body frame
%   thermalRecoil        - Per-panel thermal re-emission recoil force
%   assemble             - Reductions: F = sum F_k and tau = sum (r_k - r_cm) x F_k
%   bodyTorques          - Gravity-gradient and magnetic body torques
%   partials             - Analytic partials of a w.r.t. r v Cd Cr; complex-step verifier [stub]
%   singleEpoch          - Compiled single-epoch kernel: gas radiation and thermal recoil with their reductions
