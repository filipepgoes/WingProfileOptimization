%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNÇÃO OBJETIVO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [output_vector] = optimization_function(input_vector)

global  mach alt asa number clasa


% Root airfoil
r0_1 = input_vector(1);
Hte_1 = input_vector(2);
phi_1= input_vector(3);
X_tcmax_1 = input_vector(4);
t_c_1 = input_vector(5);
theta_1 = input_vector(6);
epsilon_1 = input_vector(7);
X_Ycmax_1 =  input_vector(8);
Ycmax_1 =  input_vector(9);
Yt_1 = input_vector(10);
epsil1 = input_vector(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break airfoil
r0_2 = input_vector(12);
Hte_2 = input_vector(13);
phi_2= input_vector(14);
X_tcmax_2 = input_vector(15);
t_c_2 = input_vector(16);
theta_2 = input_vector(17);
epsilon_2 = input_vector(18);
X_Ycmax_2 = input_vector(19);
Ycmax_2 = input_vector(20);
Yt_2 = input_vector(21);
epsil2 = input_vector(22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tip airfoil
r0_3 = input_vector(23);
Hte_3 = input_vector(24);
phi_3 = input_vector(25);
X_tcmax_3 = input_vector(26);
t_c_3 = input_vector(27);
theta_3 = input_vector(28);
epsilon_3 = input_vector(29);
X_Ycmax_3 = input_vector(30);
Ycmax_3 = input_vector(31);
Yt_3 = input_vector(32);
epsil3 = input_vector(33);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% epsil1 = 4.83;
% epsil2 = 1.69;
% epsil3 = -0.41;

% Valores obtidos junto ao 190EBD200 (ERJ195)
% epsil1 = 5.083; 
% epsil2 = 1.721;
% epsil3 = 0.184;

[a1, b1]=trans(r0_1, Hte_1, phi_1, X_tcmax_1, t_c_1, theta_1, epsilon_1, X_Ycmax_1, Ycmax_1, Yt_1);
[a2, b2]=trans(r0_2, Hte_2, phi_2, X_tcmax_2, t_c_2, theta_2, epsilon_2, X_Ycmax_2, Ycmax_2, Yt_2);
[a3, b3]=trans(r0_3, Hte_3, phi_3, X_tcmax_3, t_c_3, theta_3, epsilon_3, X_Ycmax_3, Ycmax_3, Yt_3);

[xsing1, ysing1, trail1, slopt1, radius1, x1, yu1, yl1, t_c1, Hte1, yc1, t1] = perfil(a1,b1,number.blwf);
[xsing2, ysing2, trail2, slopt2, radius2, x2, yu2, yl2, t_c2, Hte2, yc2, t2] = perfil(a2,b2,number.blwf);
[xsing3, ysing3, trail3, slopt3, radius3, x3, yu3, yl3, t_c3, Hte3, yc3, t3] = perfil(a3,b3,number.blwf);

plot(x1,yu1,'r',x1,yl1,'r',x2,yu2,'b',x2,yl2,'b',x3,yu3,'g',x3,yl3,'g');
axis equal;
pause(0.1);

[dens,temp,vsom,visc]=atmosfera(alt.cru);
V = mach.cru1*vsom;
rey.cru1 = (dens*V*asa.cma) / visc; % Reynolds para o numero de Mach de cruzeiro

[dens,temp,vsom,visc]=atmosfera(alt.cru);
V = mach.cru2*vsom;
rey.cru2 = (dens*V*asa.cma) / visc; % Reynolds para o numero de Mach de cruzeiro

[dens,temp,vsom,visc]=atmosfera(alt.estol);
V = mach.estol*vsom;
rey.estol(1) = (dens*V*asa.CHORDroot) / visc;
rey.estol(2) = (dens*V*asa.CHORDbreak) / visc;
rey.estol(3) = (dens*V*asa.CHORDtip) / visc;
rey.estol(4) = (dens*V*asa.cma) / visc;% Reynolds para o numero de Mach de estol


disp('Partindo para a condicao de stall.');
disp('Partindo para CL=0.7');
alfa=5; 
CLLA=0.7;
al_cl=1;
[time, CD, CL, CLvectorA, CDvector] = ...
    create_file(mach.estol, rey.estol(4), xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, ...
    trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLLA, al_cl);
flagA = [-CL/CD]
if flagA~=-2
    disp('Partindo para CL=1.0');
    alfa=8; 
    CLLB=1;
    al_cl=1;
    [time, CD, CL, CLvectorB, CDvector] = ...
        create_file(mach.estol, rey.estol(4), xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, ...
        trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLLB, al_cl);
    flagB = [-CL/CD]         
    if flagB~=-2
        [runflag, CLmaxw]=xfoil_optimization_function(a1, a2, a3, b1, b2, b3, CLvectorA, CLvectorB, CLLA, CLLB);
        if runflag==0 || CLmaxw<clasa
            output_vector=-2
        else
            disp('Partindo para CL=0.35');
            alfa=1; 
            CLL=0.35;
            al_cl=1;
            [time, CD1, CL1, CLvector, CDvector] = ...
                create_file(mach.cru2, rey.cru2, xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, ...
                trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLL, al_cl);
            disp('Partindo para CL=0.45');
            alfa=1; 
            CLL=0.45;
            al_cl=1;
            [time, CD2, CL2, CLvector, CDvector] = ...
                create_file(mach.cru1, rey.cru1, xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, ...
                trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLL, al_cl);
            output_vector = 0.5*[-CL1/CD1]+0.5*[-CL2/CD2]
        end
    else
        output_vector=-2   
    end
else
    output_vector=-2
end          




    