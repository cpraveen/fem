from pylab import *

gasGam = 1.4;
gasR = 287.0;

mach_inf = 0.5;
P0_inf = 102010.0;
T0_inf = 288.6;
factor = 1.0 + 0.5*(gasGam - 1.0)*mach_inf**2;
Tem_inf = T0_inf/factor;
pre_inf = P0_inf/factor**(gasGam/(gasGam - 1.0));
rho_inf = pre_inf/(gasR*Tem_inf);
vel_inf = mach_inf*sqrt(gasGam*gasR*Tem_inf);
pre_out = 101300.0;

print("pre_inf = ", pre_inf)
print("pre_out = ", pre_out)
