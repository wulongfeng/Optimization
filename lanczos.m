function a=lanczos(B)
% Calculate the Lanczos decomposition
m = size(B, 1);
max_iterations = m; % Adjust as needed
Q = zeros(m, max_iterations);
T = zeros(max_iterations, max_iterations);

q = rand(m, 1); % Random initial vector
q = q / norm(q);

for k = 1:max_iterations
    Q(:, k) = q;
    z = B * q;
    alpha = q' * z;
    z = z - alpha * q;

    if k > 1
        z = z - beta * Q(:, k - 1);
    end

    T(k, k) = alpha;
    if k < max_iterations
        beta = norm(z);
        T(k + 1, k) = beta;
        T(k, k + 1) = beta;
        q = z / beta;
    end
end

% Calculate the eigenvalues of the tridiagonal matrix T
eigenvalues_T = eig(T);

% Check for negative eigenvalues in T
negative_eigenvalues = sum(eigenvalues_T < 0);

% Determine the value of a based on the negative eigenvalues
if negative_eigenvalues > 0
    a = abs(min(eigenvalues_T)) + 1; % Set a to be greater than the most negative eigenvalue
else
    a = 0; % No negative eigenvalues, no need to adjust a
end


% Calculate B + aI
% B_aI = B + a * eye(m);
% Check if B_aI is positive definite
% is_positive_definite = all(eig(B_aI) > 0);
% if is_positive_definite
%     disp('positive definite after transition');
% else
%     disp('not positive definite');
% end
end