function rhs = right_hand_side_ddde(t, inp, p, varargin)

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

% Evaluate the delay quantities
r = h(x, p, varargin{:});

% Evaluate the delayed states
z = C*Z;

% Allocate memory
rhs = zeros(numel(inp), 1);

% Evaluate right-hand side function of delay differential equation
rhs(1:nx) = f(t, x, z, p, varargin{:});

% Evaluate the right-hand side in the auxiliary states
rhs(nx+1:end) = A*Z + B*r;