function [] = NH3syn_isothermalPFR_900F 

clc
clear

% initial conditions

initialconversion = [0]; 
Vrange = [0 4000]; % ft^3

% call runge kutta algorithm [ode45]

[V,x] =ode45(@diffeq,Vrange,initialconversion);

Vactual = V.*268./(x+10^(-8)); % ft^3 for 40,000 ton/yr plant
logVactual =log10(Vactual);

% create output table
table1 = [x,Vactual*10^-4]

% create output figure
figure (1)
plot(x,logVactual)

xlabel('x conversion')
ylim([-3 6])
ylabel('log10(Vactual),ft3')
text (0.1,5,'{900 F; 300 atm}')


end

function dxdV = diffeq (V,x)

% constants

k1=1.2; % lbmol/ft^3/h
beta = 0.00140; 
K = 0.00467;
P = 300; % atm
N2o = 1 % inlet N2 flow lbmol/hr

% differential eqns
dxdV = zeros(1,1);

yN=(1-x)/(4.-2.*x)
yH=3.*yN
yNH3=2.*x/(4.-2.*x)

dxdV(1)=N2o^-1*k1*(yN*yH^3-yNH3^2/(K^2*P^2))/((yNH3*yH^0.5+beta*yH^2)^1.5);



end






