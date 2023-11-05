function p=twoDim(g, B, radius)
eigenvalues_B = eig(B);
min_eigenvalues_B = min(eigenvalues_B);

% B is positive definite: dogleg or cauchy point
if min_eigenvalues_B > 0
    fprintf('B is positive definite \n');
    p = dogLeg(g, B, radius);
% B has zero eigenvalues but no negative eigenvalues: cauchy point
elseif min_eigenvalues_B == 0 
%if min_eigenvalues_B >=0
    fprintf('the minimal eigenvalue of B is 0 \n');
    if g'*B*g <=0
        tau=1;
    else
        tau=min(norm(g)^3/(radius*g'*B*g), 1);
    end
    p=-tau*radius/(norm(g))*g;
% B has negative eigenvalues
else
    fprintf('B has negative eigenvalue \n');
    % find a to satisfy that B + aI is positive definite
    a=lanczos(B);
    if norm((B+a*eye(size(B,1)))\g) <= radius
        % find v satisfies v'(B+aI)\g <=0
        v= find_v(B+a*eye(size(B,1)), g);
        p=-(B+a*eye(size(B,1)))\g + v;
    else
        % just add them together
        p=g+(B+a*eye(size(B,1)))\g;
    end
end

