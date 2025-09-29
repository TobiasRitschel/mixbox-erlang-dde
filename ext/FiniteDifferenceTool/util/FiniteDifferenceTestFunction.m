function [f, df, d2f] = FiniteDifferenceTestFunction(x)
% FiniteDifferenceTestFunction  Currently undocumented

% Unpack argument
x1 = x(1);
x2 = x(2);

% Evaluate function
f = [sin(x1)*sin(x2), cos(x1)*sin(x2); x1^2*x2, x1*log(x2)];

% Evaluate first order derivatives
df = [cos(x1)*sin(x2),   sin(x1)*cos(x2);
      2*x1*x2,           x1^2;
      -sin(x1)*sin(x2),  cos(x1)*cos(x2);
    log(x2), x1/x2];

% Evaluate second order derivatives
d2f(:,:,1) = [-sin(x1)*sin(x2), cos(x1)*cos(x2);
    cos(x1)*cos(x2), -sin(x1)*sin(x2)];
d2f(:,:,2) = [2*x2, 2*x1;
    2*x1, 0];
d2f(:,:,3) = [-cos(x1)*sin(x2), -sin(x1)*cos(x2);
    -sin(x1)*cos(x2), -cos(x1)*sin(x2)];
d2f(:,:,4) = [0, 1/x2;
    1/x2, -x1/x2^2];    