% Perfil com as equaçoes que os russos propuseram
function [xsing, ysing, trail, slopt, radius, x, yu, yl, t_c, Hte, yc, t] = perfil(a,b,number)

x = linspace(0,1,number); % pontos na corda

t = a(1)*x.^0.5 + a(2)*x + a(3)*x.^2 + a(4)*x.^3 + a(5)*x.^4;
yc = b(1)*x + b(2)*x.^2 + b(3)*x.^3 + b(4)*x.^4  + b(5)*x.^5 + b(6)*x.^6;


yu = t + yc;
yl = yc - t;
t_c = max(yu)-min(yl);
Hte = yu(end)-yl(end);

radius= (a(1)^2)/2;
ysing = yu(1);
xsing = x(1) + 0.5*radius;

upperangle = (yu(end-1)-yu(end))/(x(end)-x(end-1));
lowerangle = (yl(end-1)-yl(end))/(x(end)-x(end-1));

trail = (upperangle - lowerangle) * (180/pi);
slopt = (yc(end-1)-yc(end))/(x(end-1)-x(end));