function v=find_v(B,g)
% Calculate the Lanczos decomposition
[V, D] = eig(B);
v = -V(:, 1);
if v'*B*g>0
    v=-v
end
end