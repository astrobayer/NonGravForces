% NONGRAVFORCES  Non-gravitational forces, torques and partials for LEO satellites.
% Version 0.8.3  08-Sep-2026
%
% Single point of access
%   Model              - Configure once, evaluate many: accel | torque | series |
%                        partials | rhs | budget
%   Options            - Validated option set with presets fast/standard/precise
%
% Core (L0)
%   const              - Constants owned by this toolbox (S0, sigma_SB, species)
%   toolboxVersion     - Semantic version, one authority
%   dependencyReport   - Presence and version of the four ecosystem toolboxes
%   fail               - The one place errors are raised from (identifier
%                        taxonomy, correction hints)
%   buildMex           - Build the MEX gateways in mex/ (msis, kernel; never
%                        committed) or report their state
%   provenance         - Metadata record that travels with every result
%
% Layers (see doc)
%   io      - readers over gracefoReading, space weather, mass, CERES, IGRF
%   env     - density, wind, eclipse, Earth radiance, magnetic field, sources
%   geom    - macro models, self-shadowing, thermal state, attitude provider
%   kernel  - per-panel gas-surface and radiation forces, body torques,
%             partials, and the compiled single-epoch path
%   sensor  - measurement models: electrostatic, atom interferometer, drag-free
%   qc      - outlier, impulse and jump detection on residuals
%   app     - calibration (own linear solver), force budget, density inversion
%
% Dependencies (called by qualified name, never copied or wrapped)
%   referenceSystem >= 2.1   time, frames, ephemerides, attitude kinematics
%   gracefoReading  >= 2.1   GRACE-FO Level-1B reading
%   gjToolbox (odeSuite)     orbit integration (optional, for rhs)
