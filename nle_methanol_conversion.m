function nle_methanol_conversion

% Set initial guess
xo = [0.5];

% Call fsolve
[result, fval] = fsolve(@eqns, xo);
result % Display the zero
fval % Display value of the function at the "zero"
end

function f = eqns(x)

% Define constants

K = 0.006; 
p = 1; % atm

% Define function
x = x(1);  % Let x = x(1)

f(1) = x*(3-2*x)^2/(4*(1-x)^3)-K*p^2;

end


