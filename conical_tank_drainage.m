function [] = conical_tank_drainage ()

clc
clear

% initial conditions

initialheight = [5]; % ft
trange = [0 100]; % min

% call runge kutta algorithm [ode45]

[t,h] =ode45(@diffeq,trange,initialheight);

% create output table

table1 = [h,t]

% create output figure

figure (1)
plot(t,h)

xlabel('t,min')
ylim([0 8])
ylabel('h,ft')
text (20,7,'{conical tank drainage}')


V = (9/75)*3.1416*h.^3; % ft3
Vf = 0.25*(9/75)*3.1416*5^3
table2 = [V,t]

end

function dhdt = diffeq (t,h)

% constants


% differential eqns
dhdt = zeros(1,1);

dhdt(1) = -0.02*25*(2+h^2)/(9*3.1416*h^2);

end






