function [a, b]=trans(r0, Hte, phi, X_tcmax, t_c, theta, epsilon, X_Ycmax, Ycmax, Yt)

a5 = 1/2*(2*t_c*X_tcmax^(1/2)+2^(1/2)*r0^(1/2)*X_tcmax^(7/2)-5*2^(1/2)*r0^(1/2)*X_tcmax^(5/2)+2*phi*X_tcmax^(7/2)-Hte*X_tcmax^(7/2)-2*phi*X_tcmax^(5/2)+3*Hte*X_tcmax^(5/2)+5*2^(1/2)*r0^(1/2)*X_tcmax^2-2^(1/2)*r0^(1/2)*X_tcmax-6*t_c*X_tcmax^(3/2))/X_tcmax^(5/2)/(3*X_tcmax-3*X_tcmax^2-1+X_tcmax^3);
a4 = (1/2*2^(1/2)*r0^(1/2)/X_tcmax^(1/2)+Hte-phi-3/2*2^(1/2)*r0^(1/2)+2*a5+(2*phi+2^(1/2)*r0^(1/2)-Hte-6*a5)*X_tcmax+4*a5*X_tcmax^3)/(-1+4*X_tcmax-3*X_tcmax^2);
a3 = phi + 0.5*sqrt(2*r0) - (Hte/2) - 2*a4 - 3*a5;
a2 = (Hte/2) - sqrt(2*r0) - a3 - a4 - a5;
a1 = sqrt(2*r0);

a(1)=a1;
a(2)=a2;
a(3)=a3;
a(4)=a4;
a(5)=a5;

b6 = (-8*X_tcmax^3*X_Ycmax^3*theta-2*X_Ycmax*X_tcmax^4*theta+X_Ycmax*X_tcmax^3*theta+3*X_Ycmax^4*Yt+3*Ycmax*X_tcmax^2*X_Ycmax+4*X_tcmax^2*X_Ycmax^3*theta-X_Ycmax^3*Yt+8*X_Ycmax^4*theta*X_tcmax^3-X_Ycmax^4*X_tcmax^5*epsilon-2*X_Ycmax^2*X_tcmax^2*theta+4*X_Ycmax^2*X_tcmax^4*theta+X_Ycmax^3*theta*X_tcmax-X_Ycmax^4*X_tcmax^5*theta+X_tcmax^5*epsilon*X_Ycmax^3-X_tcmax^4*epsilon*X_Ycmax^4+2*X_tcmax^3*epsilon*X_Ycmax^4-3*X_Ycmax^4*theta*X_tcmax-X_tcmax^4*X_Ycmax^3*epsilon+X_Ycmax*X_tcmax^5*theta-3*X_Ycmax^2*X_tcmax^5*theta+X_Ycmax^2*X_tcmax^3*theta+3*X_Ycmax^3*theta*X_tcmax^5-4*X_Ycmax^4*theta*X_tcmax^4-3*X_Ycmax^5*Yt-2*Ycmax*X_tcmax^3+4*Ycmax*X_tcmax^4-2*Ycmax*X_tcmax^5+X_Ycmax^6*Yt-X_Ycmax^5*theta*X_tcmax^3+2*X_tcmax^4*epsilon*X_Ycmax^5-X_tcmax^3*epsilon*X_Ycmax^5+3*X_Ycmax^5*theta*X_tcmax-4*X_Ycmax^5*theta*X_tcmax^2+2*X_Ycmax^5*theta*X_tcmax^4-5*Ycmax*X_tcmax^4*X_Ycmax^2+4*Ycmax*X_tcmax^5*X_Ycmax-5*Ycmax*X_tcmax^4*X_Ycmax-5*Ycmax*X_tcmax^2*X_Ycmax^2-2*Ycmax*X_tcmax^3*X_Ycmax+10*Ycmax*X_tcmax^3*X_Ycmax^2-X_Ycmax^6*theta*X_tcmax^3-X_Ycmax^5*X_tcmax^2*epsilon+X_Ycmax^6*X_tcmax^2*epsilon-X_Ycmax^6*X_tcmax^3*epsilon-X_Ycmax^6*theta*X_tcmax+2*X_Ycmax^6*theta*X_tcmax^2)/X_tcmax^2/(1+X_tcmax^2-2*X_tcmax)/(X_Ycmax^3*X_tcmax^2-3*X_Ycmax^2*X_tcmax^2-X_tcmax^2+3*X_tcmax^2*X_Ycmax-6*X_Ycmax^2*X_tcmax+6*X_tcmax*X_Ycmax^3+2*X_Ycmax*X_tcmax-2*X_Ycmax^4*X_tcmax+X_Ycmax^5+3*X_Ycmax^3-X_Ycmax^2-3*X_Ycmax^4)/X_Ycmax^3;
b5 = -(-6*X_Ycmax^2*b6*X_tcmax^6+2*X_Ycmax*b6*X_tcmax^6+9*X_Ycmax^2*X_tcmax^2*theta+4*X_Ycmax*X_tcmax^4*theta-6*X_Ycmax*X_tcmax^4*b6-2*X_Ycmax*X_tcmax^3*epsilon+3*X_Ycmax^2*X_tcmax^2*epsilon-6*X_Ycmax^2*theta*X_tcmax+2*X_Ycmax*theta*X_tcmax-6*b6*X_Ycmax^5*X_tcmax^2+12*b6*X_Ycmax^5*X_tcmax^3-6*X_Ycmax^2*X_tcmax^2*b6-6*X_Ycmax*X_tcmax^3*theta+4*X_Ycmax*X_tcmax^3*b6-X_tcmax^2*theta+2*X_tcmax^3*theta+6*X_Ycmax^2*Yt-2*X_Ycmax*Yt-3*X_Ycmax^2*X_tcmax^4*theta+12*X_Ycmax^2*X_tcmax^4*b6-X_tcmax^4*theta-4*X_Ycmax^3*Yt-3*X_Ycmax^2*X_tcmax^4*epsilon+2*X_Ycmax*X_tcmax^4*epsilon+4*X_tcmax^3*X_Ycmax^3*theta-16*X_tcmax^3*X_Ycmax^3*b6-8*X_tcmax^2*X_Ycmax^3*theta+4*X_tcmax^3*X_Ycmax^3*epsilon+12*X_Ycmax^3*X_tcmax^2*b6+4*X_Ycmax^3*theta*X_tcmax-4*X_Ycmax^3*X_tcmax^2*epsilon+4*X_Ycmax^3*b6*X_tcmax^6-6*b6*X_Ycmax^5*X_tcmax^4)/X_tcmax^2/(1+X_tcmax^2-2*X_tcmax)/(2*X_tcmax-6*X_Ycmax*X_tcmax+4*X_Ycmax^2*X_tcmax-3*X_Ycmax+8*X_Ycmax^2-5*X_Ycmax^3)/X_Ycmax;
b4 = -(-2*X_tcmax^2*theta+X_tcmax^3*theta+3*X_tcmax^2*b6-4*X_tcmax^3*b6+theta*X_tcmax+2*X_tcmax^2*b5-X_tcmax^2*epsilon+b6*X_tcmax^6-Yt+b5*X_tcmax^5-3*X_tcmax^3*b5+X_tcmax^3*epsilon)/X_tcmax^2/(1+X_tcmax^2-2*X_tcmax);
b3 = (-2*theta*X_tcmax+X_tcmax^2*theta+4*X_tcmax^2*b4+5*X_tcmax^2*b5+6*X_tcmax^2*b6-X_tcmax^2*epsilon-2*b6*X_tcmax^6+2*Yt-2*b4*X_tcmax^4-2*b5*X_tcmax^5)/X_tcmax^2/(-3+2*X_tcmax);
b2 = -1/2*theta-3*b6-3/2*b3-2*b4-5/2*b5+1/2*epsilon;
b1 = theta;

b(1)=b1;
b(2)=b2;
b(3)=b3;
b(4)=b4;
b(5)=b5;
b(6)=b6;