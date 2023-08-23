function [xvals]=newtoniteration_example(f,g)

% Constants and initialization

x = 1.5;   % Set starting value
nmax = 25;  % Set aximum number of iterations
eps = 1;  % Initialize error bound, eps
xvals = x;  % Initialize array of iterates
n = 0;  % Initialize iteration counter
error = 1e-5;  % Define tolerance

%  Newton iteration

while eps >= error & n<=nmax  % While error is larger than the tolerance and n <= the maximum number of iterations
       
    y = x - f(x) / g(x);  %  Newton iteration
    
    xvals = [xvals;y];  % Write next iterate in array
    eps = abs(y-x);  % Compute error
    x = y; n = n + 1;  % Update x and n
end
end
