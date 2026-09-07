% KERNEL  NonGravForces.kernel - L4 layer.
% Version 0.5.1  07-Sep-2026
%
%   gsiCoefficients      - Pressure and shear coefficients per panel: sentman cll maxwell schamberg
%   accommodation        - Energy accommodation coefficient: fixed sesam langmuir
%   dragPanelForce       - Per-panel gas-momentum force in body frame
%   radiationPanelForce  - Per-panel radiation force for a source list in body frame
%   thermalRecoil        - Per-panel thermal re-emission recoil force
%   assemble             - Reductions: F = sum F_k and tau = sum (r_k - r_cm) x F_k
%   bodyTorques          - Gravity-gradient and magnetic body torques
%   partials             - Analytic partials of a w.r.t. r v Cd Cr; complex-step verifier [stub]
