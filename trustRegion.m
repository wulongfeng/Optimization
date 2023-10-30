%function x=trustRegion(x, eta, max_radius, tol)
function x=trustRegion(max_radius, tol)
x=[-1.2, 1];
eta=0.15;
delta = 0.1;
g=grad_f(x);
k=0; % k = # iterations

while  norm(g) > tol    
%for iter = 1:50
    k = k + 1;
    fprintf('\n\n\n');
    fprintf('Iteration %d:\n', k);
    [g,B] = update_g_B(x);
    p=dogLeg(g, B, delta);
    redu=func(x)-func(x+p');
    pred_redu= -(g'*p + 0.5*p'*B*p);
    rho = redu / pred_redu;
    %fprintf('Iteration: %d \tp: %f\trho: %.3e\n', k, p, rho)
    fprintf('p = [%f, %f]\n', p(1), p(2));
    fprintf('redu= %f, pred_redu=%f, rho = %f\n', redu, pred_redu, rho);
    
    fprintf('norm p = %f\n', norm(p));
    if rho <= 0.25
        %fprintf('delta = %f\n', delta);
        %fprintf('p norm = %f\n', norm(p));
        delta = 0.25 * norm(p);
    %elseif rho > 0.75 && norm(p) == delta < eps
    elseif rho > 0.75 && norm(p) == delta
        delta = min(2*delta, max_radius);
    end
    fprintf('delta = %f\n', delta);


    if rho > eta
    x = x+p';
    end
    fprintf('x = [%f, %f]\n', x(1), x(2));
end
end

