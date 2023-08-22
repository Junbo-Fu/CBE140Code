function NH3_condenser

% 1=N2,2=H2,3=NH3,4=A
%set constants
global z K 
z=[0.2099;0.6298;0.1050;0.0553];
K=[66.7;50;0.0235;100];

% Set initial guess
VoFo = [0.5];

% Call fsolve
[VoF, fval] = fsolve(@eqns, VoFo);
VoF % Display the zero
fval % Display value of the function at the "zero"


x = z./(1 + VoF*(K-1))
y = z.*K./(1 + VoF*(K-1))
% L = 1 1bmol as basis
F=1/(1-VoF)
V=VoF*F


end

function f = eqns(VoF)

global z K

% Define function
S=0;
for i=1:4
    S = S + z(i)*(K(i)-1)/(1 + VoF*(K(i)-1));
    
end

f = S;

end