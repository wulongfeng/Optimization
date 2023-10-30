% Define the gradient of the objective function
function g = grad_f(x)
g(1) = 100*(2*(x(1).^2-x(2))*2*x(1)) + 2*(x(1)-1);
g(2) = 100*(-2*(x(1).^2-x(2)));
g = g';
end