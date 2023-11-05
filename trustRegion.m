function x=trustRegion(tol)
problem='Rosenbrock';
%problem='test2';

max_radius=1;
eta=0.2;
delta = 0.1;
k=0; % k = # iterations

% initilization
% For Rosenbrock problem
if strcmp('Rosenbrock', problem)
    x=[-1.2; 1];
    g=grad_f(x);
    func = @(x) 100*(x(2) - x(1).^2).^2 + (1-x(1)).^2;
% For test2.m
elseif strcmp('test2', problem)
    G = numgrid('L',32); % finite difference grid for L-shaped domain
    A = delsq(G);        % finite difference matrix for L-shaped domain
    n = size(A,1);
    b = rand(n,1);
    c = 4*ones(n,1);
    func = @(x) (0.5*(x'*A*x) - b'*x)^2 + 0.5*(x'*A*x) - c'*x; % objective function
    gf = @(x) 2*(0.5*(x'*A*x) - b'*x)*(A*x-b) + (A*x-c);     % gradient
    g2f = @(x) A*(1+2*(0.5*(x'*A*x) - b'*x)) + 2*(A*x-b)*(A*x-b)'; % Hessian
    x = ones(n,1);
    g=gf(x);
end

while  norm(g) > tol    
%for iter = 1:30
    k = k + 1;
    fprintf('\n\n\n');
    fprintf('Iteration %d:\n', k);
    
    % For Rosenbrock problem
    if strcmp('Rosenbrock', problem)
        g=grad_f(x);
        h=hessian_f(x);
    % For test2.m
    elseif strcmp('test2', problem)
        g=gf(x);
        h=g2f(x);
    end
    
    % update B based on Cholesky decomposition, modify block-diagonal matrix D
    B = update_B(h);

    p=dogLeg(g, B, delta);
    %p=twoDim(g, B, delta);
	%p=twoDimUpdate(g, B, delta);

    redu=func(x)-func(x+p);
    pred_redu= -(g'*p + 0.5*p'*B*p);
    rho = redu / pred_redu;
    fprintf('p = [%f, %f]\n', p(1), p(2));
    fprintf('redu= %f, pred_redu=%f, rho = %f\n', redu, pred_redu, rho);
    
    fprintf('norm chosen p = %f, delta=%f, abs diff=%f\n', norm(p), delta, abs(norm(p)-delta));
    if rho <= 0.25
        delta = 0.25 * norm(p);
    elseif rho > 0.75 && abs(norm(p)-delta) < 1e-3
        delta = min(2*delta, max_radius);
    end

    fprintf('Updated delta = %f\n', delta);

    if rho > eta
    x = x+p;
    end
    fprintf('Updated x = [%f, %f]\n', x(1), x(2));
end
end

