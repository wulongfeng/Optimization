function p=twoDimUpdate(g, B, radius)
pb=-B\g;
pu=-(g'*g)/(g'*B*g)*g;
pu_norm = norm(pu);

p1=pu/pu_norm;
p2=(pb-p1*p1'*pb)/norm(pb-p1*p1'*pb);
disp('Updated p1:');
disp(p1);
disp('Updated p2:');
disp(p2);
% p1 and p2 are orthogonal, and have been normalized
update_B = [p1'*B*p1, p1'*B*p2;
            p2'*B*p1, p2'*B*p2];
update_g = [p1'*g; p2'*g];

disp('Updated B:');
disp(update_B);
disp('Updated g:');
disp(update_g);

update_p = -update_B\update_g;
disp('Updated p:');
disp(update_p);

if norm(update_p) <= radius
    p = update_p;
else
    lambda = find_lambda(update_g, update_B, radius, 1e-3);
    disp('the value of found lambda:');
    disp(lambda);
    p = -(update_B + lambda *eye(size(B,1)))\update_g;
end
end