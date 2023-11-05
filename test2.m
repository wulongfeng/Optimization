% script for test2.m

%% (potentially) large 4-th order problem
%
%  f(x) = [ (1/2)x^T Ax - x^Tb ]^2 + (1/2)x^T Ax - x^T c
%
 
G = numgrid('L',32); % finite difference grid for L-shaped domain
A = delsq(G);        % finite difference matrix for L-shaped domain

n = size(A,1);
b = rand(n,1);
c = 4*ones(n,1);
f = @(x) (0.5*(x'*A*x) - b'*x)^2 + 0.5*(x'*A*x) - c'*x; % objective function
gf = @(x) 2*(0.5*(x'*A*x) - b'*x)*(A*x-b) + (A*x-c);     % gradient
g2f = @(x) A*(1+2*(0.5*(x'*A*x) - b'*x)) + 2*(A*x-b)*(A*x-b)'; % Hessian
x0 = ones(n,1);