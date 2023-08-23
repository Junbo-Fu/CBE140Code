function [] = tow_example_of_couple_odes ()
% This function uses experimental rate constants as a function of
% temperature for the calculation of the concentrations of A, B, C, D, E,
% and F as a function of time, from a given reaction network.

%Thomas Dursch
%September 29th, 2009

%The reaction network is: A + B --> C --> D; A + D <---> E; 2E --> F; 
%which all occur in parallel. 

clc 
clear all

%Defining constants:
%Pre-exponential Factors:
k1p = 100000; % [=] mM^-1*s^-1
k2p = 13000; % [=] s^-1
k3p = 500; % [=] mM^-1*s^-1
k3rp = 22000; % [=] s^-1
k4p= 1.6*10^8; % [=] mM^-1*s^-1

%Actvation Energies [=] J/mol:
Ea1 = 44000;
Ea2 = 28000;
Ea3 = 21000;
Ea3r = 35000;
Ea4 = 66000;

R = 8.314; % Gas Constant [=] J/mol-K
T = 328; % Temperature at which the experiment was run [=] K

%Assume the rate constants have the Arrhenius Form: k(T) =
%kxp*exp(-Ea/(R*T)):

global k1 k2 k3 k3r k4 C

k1 = k1p*exp(-Ea1/(R*T)); % A + B --> C
k2 = k2p*exp(-Ea2/(R*T)); % C --> D
k3 = k3p*exp(-Ea3/(R*T)); % A + D --> E
k3r = k3rp*exp(-Ea3r/(R*T)); % E --> A + D
k4 = k4p*exp(-Ea4/(R*T)); % 2E --> F

%Initial Conditions [=] mMol:
initial = [100 100 0 0 0 0];
tspan = [0 10000];

%ODE Solver: 
[time, c] = ode45(@diffeq, tspan, initial); % Use ODE45 to solve the 
%system of rate equations defined by '@diffeq' over the given timespan and 
%using the intial conditions.

figure(1)
plot(time, c(:,1),'r') % Create a plot of the concentrations as a function
%of time.  

hold on % Not only plot A, but plot B, C, D, E, and F as well.

plot(time, c(:,2),'m')
plot(time, c(:,3),'c')
plot(time, c(:,4),'g')
plot(time, c(:,5),'b')
plot(time, c(:,6),'k')

% Additionally, print the values of the concentrations at t=tspan so that
% the final value is easily determined.
C(1)
C(2)
C(3)
C(4)
C(5)
C(6)

hold off

xlabel('Time [=] s')
ylabel('Concentration [=] mM')
legend('A', 'B', 'C', 'D', 'E', 'F')

end % Ending the function in order to create a new function diffeq.

function df = diffeq (t,f)
%This function is used in the ODE solver - it is the rate equations that
%will be solved by the ODE solver.

global k1 k2 k3 k3r k4 C

df = zeros(6,1); % Creates a zeros vector where the solutions to the ODEs 
%will be stored.  It is a 6 species x 1 solution vector.

%Rate Equations: f(1)=A, f(2)=B, f(3)=C, f(4)=D, f(5)=E, f(6)=F
df(1) = -k1*f(1)*f(2)-k3*f(1)*f(4)+ k3r*f(5);
df(2) = -k1*f(1)*f(2);
df(3) = k1*f(1)*f(2) - k2*f(3);
df(4) = k2*f(3) - k3*f(1)*f(4) + k3r*f(5);
df(5) = k3*f(1)*f(4) - k3r*f(5) - 2*k4*f(5)*f(5);
df(6) = k4*f(5)*f(5);

%Define f(i) as the vector C in order to report the values at time t=tspan.
C(1) = f(1);
C(2) = f(2);
C(3) = f(3);
C(4) = f(4);
C(5) = f(5);
C(6) = f(6);

end
