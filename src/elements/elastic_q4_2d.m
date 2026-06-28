function [Ke, fe] = elastic_q4_2d(coords, E, nu, beta, dc)
%ELASTIC_Q4_2D  Q4 bilinear plane-strain elastic element with chemical eigenstrain.
%
%   [Ke, fe] = elastic_q4_2d(coords, E, nu, beta, dc) returns the 8x8 element
%   stiffness matrix and the 8x1 eigenstrain force vector for a 4-node bilinear
%   quadrilateral, in PLANE STRAIN, with an in-plane chemical eigenstrain
%   eps0 = beta*dc applied to the x- and y-components (eps0_zz = 0). For fully
%   3-D isotropic Vegard swelling under plane strain, pass an effective
%   beta -> (1+nu)*beta (the standard thermal-analogue factor).
%
%   Inputs:
%     coords - (4x2) nodal coordinates [x y], counter-clockwise
%     E, nu  - Young's modulus [Pa], Poisson ratio
%     beta   - chemical expansion coefficient [m^3/mol]   (optional, default 0)
%     dc     - (4x1) or scalar concentration change c-c0 [mol/m^3] (optional)
%
%   DOF ordering: [u1x u1y u2x u2y u3x u3y u4x u4y].
%
%   Plane-strain constitutive matrix:
%     D = E/((1+nu)(1-2nu)) * [1-nu  nu     0;
%                              nu    1-nu   0;
%                              0     0   (1-2nu)/2]
%
%   See also: build_trilayer_mesh_2d, interface_element_2d, solve_mech_2d

if nargin < 4 || isempty(beta), beta = 0; end
if nargin < 5 || isempty(dc),   dc   = 0; end
if isscalar(dc), dc = dc * ones(4,1); end

% Plane-strain elasticity matrix
c  = E / ((1+nu)*(1-2*nu));
D  = c * [1-nu, nu,   0;
         nu,   1-nu, 0;
         0,    0,    (1-2*nu)/2];

% 2x2 Gauss quadrature
g  = 1/sqrt(3);
gp = [-g -g;  g -g;  g g;  -g g];
w  = [1 1 1 1];

Ke = zeros(8,8);
fe = zeros(8,1);

for q = 1:4
    xi = gp(q,1); eta = gp(q,2);
    % Bilinear shape functions and natural derivatives
    N    = 0.25*[(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
    dNdxi  = 0.25*[-(1-eta); (1-eta); (1+eta); -(1+eta)];
    dNdeta = 0.25*[-(1-xi); -(1+xi); (1+xi); (1-xi)];
    % Jacobian
    J = [dNdxi'; dNdeta'] * coords;      % 2x2
    detJ = det(J);
    if detJ <= 0
        error('elastic_q4_2d:badJacobian','Non-positive Jacobian (detJ=%g).', detJ);
    end
    dN = J \ [dNdxi'; dNdeta'];          % 2x4 : [dN/dx; dN/dy]
    % Strain-displacement matrix B (3x8)
    B = zeros(3,8);
    B(1,1:2:end) = dN(1,:);
    B(2,2:2:end) = dN(2,:);
    B(3,1:2:end) = dN(2,:);
    B(3,2:2:end) = dN(1,:);
    % Stiffness
    Ke = Ke + (B' * D * B) * detJ * w(q);
    % In-plane eigenstrain force: f = int B' D eps0 dV, eps0 = beta*dc*[1;1;0]
    % (eps_zz = 0; for 3-D isotropic swelling pass (1+nu)*beta -- see header).
    dc_q = N' * dc;
    eps0 = beta * dc_q * [1; 1; 0];
    fe = fe + (B' * D * eps0) * detJ * w(q);
end
end
