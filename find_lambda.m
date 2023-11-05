function lambda=find_lambda(g, B, radius, tol)
% find lambda to satisfy norm((B+lambda*I)\g)==radius
[eigVecs, eigVals] = eig(B);
eigVals = diag(eigVals);
[smallestEigVal, index] = min(eigVals);
% Extract the corresponding eigenvector
smallestEigVec = eigVecs(:, index);

% from question 3
upper_bound = -smallestEigVal + (sqrt(size(B,1))*norm(g))/radius;
lower_bound = -smallestEigVal + abs(smallestEigVec'*g)/radius;
p_upper = -(B + upper_bound *eye(size(B,1)))\g;
p_upper_diff = norm(p_upper) -radius;
p_lower = -(B + lower_bound *eye(size(B,1)))\g;
p_lower_diff = norm(p_lower) -radius;

%disp('before while-loop: p_lower_diff:');
%disp(p_upper_diff);
%disp('before while-loop: p_upper_diff:');
%disp(p_lower_diff);

interval_len = abs(upper_bound-lower_bound);
while interval_len > tol
    mid = (upper_bound + lower_bound)/2;
    p_mid = -(B + mid *eye(size(B,1)))\g;
    p_mid_diff = norm(p_mid) -radius;

    if p_mid_diff * p_upper_diff >0
        upper_bound = mid;
        p_upper_diff = p_mid_diff;
    else
        lower_bound = mid;
        p_lower_diff = p_mid_diff;
    end
   interval_len = abs(upper_bound-lower_bound);
end
lambda = (upper_bound+lower_bound)/2;
disp('The final diff after finding the appropriate lambda:');
p_lambda = -(B + lambda *eye(size(B,1)))\g;
disp(abs(norm(p_lambda) -radius));
end

