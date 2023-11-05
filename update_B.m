function B=update_B(h)
[L, D, P]=ldl(h);
%disp(h)
%disp(D)

% calculate the value of \gamma
h_froNorm = norm(h, 'fro');
gamma = 10^(-4)*h_froNorm;
%fprintf('the value of h norm: %f, gamma: %f \n', h_froNorm, gamma);
%fprintf('diff on diagnoal elements: %.17f, error:%.17f\n', all(all(abs(D - diag(diag(D))))), eps);
%fprintf('diff < error:%.1f\n', all(all(abs(D - diag(diag(D))))) < eps);
% 1*1 block (diagonal) or 2*2 block-diagonal
is_1x1_block = all(all(abs(D - diag(diag(D))))) < eps;

if is_1x1_block
    %disp("D is diagnoal matrix")
    diagonal_elements = diag(D);
    replace_indices = (diagonal_elements < gamma) | (diagonal_elements < 0);
    D(diag(replace_indices))=gamma;
else
    %disp("D is 2*2 block-diagnoal matrix")
    n=size(D, 1);
    %fprintf('the size of D: %d\n', n);
    %disp(D);
    for i = 1:2:n
        %block = D(i:i+1, i:i+1);
        D(i,i) = abs(D(i,i+1))+gamma;
        D(i+1,i+1) = abs(D(i,i+1))+gamma;
        %disp(D)
    end
end
%disp(D)
B=P'*L*D*L'*P;     
%disp(B)
end