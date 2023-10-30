function [g,B]=update_g_B(x)
g=grad_f(x);
h=hessian_f(x);
[L, D, P]=ldl(h);
%disp(D)

% calculate the value of \gamma
h_froNorm = norm(h, 'fro');
gamma = 10^(-4*h_froNorm);

% 1*1 block (diagonal) or 2*2 block-diagonal
is_1x1_block = all(all(abs(D - diag(diag(D))) < eps));
if is_1x1_block    
    diagonal_elements = diag(D);
    replace_indices = (diagonal_elements < gamma) | (diagonal_elements < 0);
    D(diag(replace_indices))=gamma;
else
    for i = 1:2:n
        %block = D(i:i+1, i:i+1);
        D(i,i) = abs(D(i,i+1))+gamma;
        D(i+1,i+1) = abs(D(i,i+1))+gamma;
        %disp(D)
    end
end
%disp(D)
B=P'*L*D*L'*P;     
end