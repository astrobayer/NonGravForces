% KERNEL  NonGravForces.kernel - L4 layer.
% Version 0.16.0  18-Oct-2026
%
%   gsiCoefficients      - Normal and tangential plate coefficients: sentman maxwell constantCd (cll/schamberg refuse)
%   accommodation        - Energy accommodation per epoch: fixed or Pilinski 2010 or the Walker/Mehta adsorbate mixture
%   dragPanelForce       - Per-panel gas-momentum force in body frame
%   radiationPanelForce  - Per-panel radiation force for a source list in body frame
%   thermalRecoil        - Per-panel thermal re-emission recoil force
%   assemble             - Reductions: F = sum F_k and tau = sum (r_k - r_cm) x F_k
%   bodyTorques          - Gravity-gradient and magnetic body torques
%   partials             - Analytic per-panel partials of the drag radiation and thermal forces
%   singleEpoch          - Compiled single-epoch kernel: gas radiation and thermal recoil with their reductions
