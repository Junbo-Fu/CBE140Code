function rising_ocean_acidity

x0=[10^(-8.1),0,0]% initial guesses for extents of reaction for reactions 3 and 4

[z] = fsolve(@eqns, x0); % solves numerically for extents of reaction for reactions 3 and 4

h_plus = 10^(-8.1)+z(1)+z(2);% [H+] based on extents of reaction
pH = -log10(h_plus)
end

function f = eqns(x0)
yco2_ppm=300;% concentration of co2 in ppm
yco2=yco2_ppm/10^6;% concentration of co2 in mole fractions
P=1;% pressure in atm
co2=yco2*P*10^(-1.47);%co2(aq);
h2co3=10^(-2.59)*co2;%[h2co3]
f(1) = (10^(-8.1)+x0(1)+x0(2))*(x0(1)-x0(2))/h2co3-10^(-3.76);%Rxn 3 equilibrium
f(2) = (10^(-8.1)+x0(1)+x0(2))*(x0(2)+x0(3))/(x0(1)-x0(2))-10^(-10.329);%Rxn 4 equilibrium
f(3) = x0(3)^.2-4.8*10^(-10); %Rxn 5 equilibrium
end
