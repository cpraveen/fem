% Solve using chebfun
%   -u'' + x^3 u = 1  in  (-3,3)
%        u(-3) = u(3) = 0
L = chebop(-3, 3);
L.op = @(x,u) -diff(u,2) + x^3*u;
L.lbc = 0; L.rbc = 0;
u = L\1; plot(u), grid on

x = linspace(-3,3,1000);
y = u(x);
d = [x', y'];
save('bvp1.txt','d','-ascii')
