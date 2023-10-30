function x=trustRegion(x, eta, max_radius, tol)

delta = 0.1;
g=grad_f(x);
h=hessian_f(x);
[L, D, P]=ldl(h);
disp(D)
h_froNorm = norm(h, 'fro');
gamma = 10^(-4*h_froNorm);
% 1*1 block or 2*2 block
is_1x1_block = all(all(abs(D - diag(diag(D))) < eps));
if is_1x1_block    
    diagonal_elements = diag(D);
    replace_indices = (diagonal_elements < gamma) | (diagonal_elements < 0);
    D(diag(replace_indices))=gamma;
else
    for i = 1:2:n
        %block = D(i:i+1, i:i+1);
        D(i,i) = abs(D(i,i+1))+1;
        D(i+1,i+1) = abs(D(i,i+1))+1;
        %disp(D)
    end
end
disp(D)
B=P'*L*D*L'*P;     
k=0; % k = # iterations


%while  norm(g) > tol    
for iter = 1:15
    k = k + 1;
    p=dogLeg(g, B, delta);
    redu=func(x)-func(x+p);
    pred_redu= -(g'*p + 0.5*p'*B*p);
    rho = redu / pred_redu;
    %fprintf('Iteration: %d \tp: %f\trho: %.3e\n', k, p, rho)
    fprintf('\n');
    fprintf('Iteration %d:\n', k);
    fprintf('p = [%f, %f]\n', p(1), p(2));
    fprintf('rho = %f\n', rho);
    fprintf('x = [%f, %f]\n', x(1), x(2));

    if rho <= 0.25
        %fprintf('delta = %f\n', delta);
        %fprintf('p norm = %f\n', norm(p));
        delta = 0.25 * norm(p);
    elseif rho > 0.75 && norm(p) == delta
        delta = min(2*delta, max_radius);
    end
    
    if rho > eta
    x = x+p;
    end
end
end

