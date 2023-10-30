function p=dogLeg(g, B, radius)
pb = -B\g;
pb_norm = norm(pb);
%disp(pb);
fprintf('pb_norm = %f\n', pb_norm);
fprintf('radius = %f\n', radius);
if pb_norm <= radius
    p=pb;
    fprintf('return pb\n');
    return;
end

pu=-(g'*g)/(g'*B*g)*g;
pu_norm=norm(pu);
%disp(pb);
fprintf('pu_norm = %f\n', pu_norm);
fprintf('radius = %f\n', radius);
if pu_norm >= radius
    p=radius*pu/pu_norm;
    fprintf('return pu\n');
    return;
end

%calcualte the value of tau
%eq = @(tau) norm(pu + (tau - 1) * (pb - pu))^2 - radius^2;
% Find the value of tau that minimizes the equation
%tau = fminunc(eq, 1);
pb_pu = pb -pu;
pb_pu_sq = pb_pu' * (pb_pu);
pu_pb_pu_sq = pu' * (pb_pu);
d = pu_pb_pu_sq^2 - pb_pu_sq*(pu'*pu-radius^2);
tau = (-pu_pb_pu_sq +sqrt(d))/pb_pu_sq +1;
fprintf('pb_pu = [%f, %f]\n', pb_pu(1), pb_pu(2));
%disp(size(pb_pu))
fprintf('pb_pu_sq = %f\n', pb_pu_sq);
%disp(size(pb_pu_sq))
fprintf('pu_pb_pu_sq = %f\n', pu_pb_pu_sq);
fprintf('d = %f\n', d);
fprintf('tau = %f\n', tau);

if 0 <= tau && tau <= 1
    p = tau * pu;
    return;
elseif 1 < tau && tau <= 2
    p = pu + (tau - 1) * (pb - pu);
end
