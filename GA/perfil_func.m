function [yc, t, yu, yl] = perfil_func(a,b,x)
t = a(1)*x.^0.5 + a(2)*x + a(3)*x.^2 + a(4)*x.^3 + a(5)*x.^4;
yc = b(1)*x + b(2)*x.^2 + b(3)*x.^3 + b(4)*x.^4  + b(5)*x.^5 + b(6)*x.^6;
yu = t + yc;
yl = yc - t;