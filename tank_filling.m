function [] = tank_filling ()

clc
clear

% initial conditions

initialheight = [3]; % ft
trange = [0 12]; % min

% call runge kutta algorithm [ode45]

[t,h] =ode45(@diffeq,trange,initialheight);

% create output table

table = [h,t]

% create output figure

figure (1)
plot(t,h)

xlabel('t,min')
ylim([0 8])
ylabel('h,ft')
text (1,5.5,'{transient tank filling}')

end

function dhdt = diffeq (t,h)

% constants
Qin = 4; % ft3/min
A = 3.1416*4^2/4; % ft2


% differential eqns
dhdt = zeros(1,1);

dhdt(1) = Qin/A;

end






