function [] = NH3_pfr_reactor_inerts_hw_28_22_()

clc
clear

% initial conditions

initialconversion = 0; 
Vrange = [0 100]; % ft^3

% call runge kutta algorithm [ode45]

[V,x] =ode45(@diffeq,Vrange,initialconversion);

Vactual = V*2984./x; % ft^3 for 100,000 ton/yr plant
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

k1=0.125; % lbmol/ft^3/h
beta = 0.00140; 
K = 0.00467;
P = 300; % atm

% differential eqns
dxdV = zeros(1,1);

yN=0.225*(1-x);
yH=3*yN;
yNH3=0.45*x;

dxdV(1)=k1*(yN*yH^3-(yNH3^3/(K^2*P^2)))/((yNH3*yH^0.5+beta*yH^2)^1.5)/0.225;


end






