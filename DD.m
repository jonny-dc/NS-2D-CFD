function A = DD(n,h)
% DD(n,h)
%
% One-dimensional finite-difference derivative matrix 
% of size n times n for second derivative:
%
% This function belongs to project_main.m

e = ones(n, 1);
A = spdiags([e, -2*e, e], -1:1, n, n);
A(1, 1) = -1;
A(n, n) = -1;
A = A ./ (h^2);
