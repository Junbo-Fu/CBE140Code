function eigens_eye_surface_temp 

% Set initial guess
xo = [8.];

% Call fsolve
[result, fval] = fsolve(@eqns, xo);
result % Display the zero
fval % Display value of the function at the "zero"
end

function f = eqns(x)

% Define constants
B = 0.3;

% Define function
x = x(1);  % Let x = x(1)

f(1) = tan(x)+x/B;

end

