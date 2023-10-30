function p=dogLeg(g, B, radius)
pb = -B\g;
pb_norm = norm(pb);
if pb_norm <= radius
    p=pb;
    return;
end

pu=-(g'*g)/(g'*B*g)*g;
pu_norm=norm(pu);
if pu_norm >= radius
    p=radius*pu/pu_norm;
    return;
end

%obtaining tau
eq = @(tau) norm(pu + (tau - 1) * (pb - pu))^2 - radius^2;
% Find the value of tau that minimizes the equation
tau = fminunc(eq, 1);
fprintf('tau = %f\n', tau);

if 0 <= tau && tau <= 1
    p = tau * pu;
    return;
elseif 1 < tau && tau <= 2
    p = pu + (tau - 1) * (pb - pu);
end
