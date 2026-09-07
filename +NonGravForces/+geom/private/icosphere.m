function [verts, faces] = icosphere(nSubdivisions)
%ICOSPHERE Unit icosphere: vertices (3xV) and triangular faces (3xF) after n loop subdivisions of an icosahedron.
%
%   [verts, faces] = icosphere(n)   V = 10*4^n + 2: 12, 42, 162, 642, 2562, ...
%
%   Private helper of geom.selfShadow: the look-up table over incidence
%   directions is sampled at these vertices and interpolated barycentrically
%   inside these faces. Faces are wound consistently (outward).
%
% INPUT ARGUMENTS
%   Name            Size  Type    Units  Description
%   ----            ----  ----    -----  -----------
%   nSubdivisions   1x1   double  -      0..7 (V grows by 4 per level)
%
% OUTPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   verts   3xV   double  -      Unit vectors
%   faces   3xF   double  -      Vertex indices of each triangle, F = 20*4^n
%
% ERRORS
%   None: the caller validates the level.
%
% KEYWORDS private icosphere mesh
% Created 2026-09-07 by Claude Fable 5.1 (Anthropic).
t = (1 + sqrt(5)) / 2;
verts = [-1  1 -1  1 0 0 0 0 t t -t -t; t t -t -t -1 1 -1 1 0 0 0 0; 0 0 0 0 t t -t -t -1 1 -1 1];
verts = verts ./ vecnorm(verts);
faces = [0 11 5; 0 5 1; 0 1 7; 0 7 10; 0 10 11; 1 5 9; 5 11 4; 11 10 2; 10 7 6; 7 1 8; ...
         3 9 4; 3 4 2; 3 2 6; 3 6 8; 3 8 9; 4 9 5; 2 4 11; 6 2 10; 8 6 7; 9 8 1]' + 1;
for level = 1:nSubdivisions
    nV = size(verts, 2); nF = size(faces, 2);
    % midpoint of every edge, shared between neighbouring faces through a key
    edges = [faces([1 2], :), faces([2 3], :), faces([3 1], :)];
    edges = sort(edges, 1);
    [ue, ~, ic] = unique(edges', 'rows');
    mid = verts(:, ue(:, 1)) + verts(:, ue(:, 2));
    mid = mid ./ vecnorm(mid);
    verts = [verts, mid]; %#ok<AGROW>
    m = reshape(ic, nF, 3)' + nV;                    % midpoint index per face edge: (1,2) (2,3) (3,1)
    faces = [ [faces(1, :); m(1, :); m(3, :)], [faces(2, :); m(2, :); m(1, :)], ...
              [faces(3, :); m(3, :); m(2, :)], [m(1, :); m(2, :); m(3, :)] ];
end
end
