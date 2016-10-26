function [dens, temp, vsom, visc]=atmosfera(h)

if h <= 11000 & h>-400
    temp = 288.15 - 6.5e-3*h;
    dens = 1.225*(temp/288.15)^4.256;    
elseif h>11000 & h<=20000
    temp = 216.65;
    dens = 0.3639*exp(-1.577e-4*(h-11000));    
else
    disp('Not valid altitude.');
    pause;
end    
vsom =sqrt(1.4*287*temp);
visc=0.000001458*temp^1.5*(1/(temp+110.4)); % relação de Sutherland.
