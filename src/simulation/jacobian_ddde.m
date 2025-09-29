function J = jacobian_ddde(t, inp, p, varargin)

% Dimensions of variables
nx = p.nx;

% Extract parameters
f = p.f;
h = p.h;

A = p.A;
B = p.B;
C = p.C;

% Extract states
x = inp(   1:nx );
Z = inp(nx+1:end);

% Delayed states
z = C*Z;

% Evaluate right-hand side function of delay differential equation
[~, dfx, dfz] = f(t, x, z, p, varargin{:});

% Evaluate the delay quantities
[~, dhx] = h(x, p, varargin{:});

% Jacobian of the delayed states
dz = C;

% Create Jacobian matrix
J = [dfx, dfz*dz; B*dhx, A];

% % Allocate memory
% J = zeros(numel(inp), numel(inp));
% 
% % Insert Jacobians of the right-hand side function
% J(1:nx,    1:nx ) = dfx;
% J(1:nx, nx+1:end) = dfz*dz;
% 
% % Insert Jacobians of auxiliary equations
% J(nx+1:end,    1:nx ) = B*dhx;
% J(nx+1:end, nx+1:end) = A;