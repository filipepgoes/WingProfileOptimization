clc;
clear all;
global mach alt asa number ntr clasa


%%%%%%%%%%%%%%%%%%%%%%ALGUMAS CONDIÇÕES INICIAIS%%%%%%%%%%%%%%%%%%%%%%%%
mach.cru1=0.72;  % Condição de Mach para o cruzeiro                   %%%
mach.cru2=0.76;
mach.estol=0.30;  %Condição de Mach para a situação de estol         %%%
alt.cru=9144; % [m] Altitude de cruzeiro.                            %%%  
alt.estol= 3000;  %[m] Altitude para a qual é calculado o estol      %%%
number.blwf=40;                                                      %%%
number.xfoil=40;                                                    %%%
ntr=8;                                                               %%%
clasa=1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DADOS DA PLANTA DA ASA                                                                  %%                                                               
asa.sweep_angle1 = 35.706; % [deg] - enflechamento da asa 1.                              %%
asa.sweep_angle2 = 26.706; % [deg] - enflechamento da asa 1.                              %%
asa.dihedral_angle = 5; % [deg] - angulo de diedro.                                       %%
asa.semi_span = 13.88714; % [m] - semi-envergadura da asa (b/2).                          %%
asa.span = asa.semi_span*2; % [m] - envergadura da asa (b).                               %%
                                                                                          %%
% DADOS DAS CORDAS E DISTANCIAS                                                           %%
asa.Zroot_symmetry = 0;                                                                   %%
% CHORDroot_symmetry = 6.18966; % [m] - corda da raiz no eixo de simetria.                %%
asa.Zroot = 1.3853; % [m] - distancia do eixo de simetria ao perfil da raiz.              %%
% CHORDroot = 5.495420; % [m] - corda da raiz.                                            %%
asa.Zbreak = 5.13836; % [m] - distancia do eixo de simetria ao perfil da quebra.          %%
asa.CHORDbreak = 3.62733; % [m] - corda da quebra.                                        %%
asa.Ztip = asa.semi_span; % [m] - distancia do eixo de simetria ao perfil da ponta da asa.%%
asa.CHORDtip = 1.46042; % [m] - corda da ponta da asa.                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SIMMETRY ROOT AIRFOIL PARAMETERS
asa.XLEroot_symmetry = 0;
asa.YLEroot_symmetry = 0;
asa.THICKroot_symmetry = 1;
asa.FSECroot_symmetry = 1;

% ROOT AIRFOIL PARAMETERS
asa.XLEroot = asa.Zroot * tan(asa.sweep_angle1*(pi/180));
asa.YLEroot = asa.Zroot * tan(asa.dihedral_angle*(pi/180));
asa.THICKroot = 1;
asa.FSECroot = 1;

% BREAK AIRFOIL PARAMETERS
asa.XLEbreak = asa.Zbreak * tan(asa.sweep_angle1*(pi/180));
asa.YLEbreak = asa.Zbreak * tan(asa.dihedral_angle*(pi/180));
asa.THICKbreak = 1;
asa.FSECbreak = 1;

% leva em consideraçao bordo de fuga perpendicular a fuselagem na junçao asa-fuselagem.
asa.CHORDroot_symmetry = asa.XLEbreak + asa.CHORDbreak;
asa.CHORDroot = asa.CHORDroot_symmetry - asa.XLEroot; 

% TIP AIRFOIL PARAMETERS
asa.XLEtip = asa.XLEbreak +((asa.semi_span-asa.Zbreak)* tan(asa.sweep_angle2*(pi/180)));
asa.YLEtip = asa.semi_span * tan(asa.dihedral_angle*(pi/180));
asa.THICKtip = 1;
asa.FSECtip = 1;

% swing - area da asa.
asa.semiarea1 = ((asa.CHORDtip + asa.CHORDbreak)/2)*(asa.semi_span-asa.Zbreak);
asa.semiarea2 = ((asa.CHORDbreak + asa.CHORDroot_symmetry)/2)*(asa.Zbreak);
asa.sweep_angleLE2 = atan( ((asa.XLEtip+asa.CHORDtip)-(asa.XLEbreak+asa.CHORDbreak))/(asa.semi_span-asa.Zbreak) );
asa.semiarea3 = ((asa.CHORDbreak + (asa.CHORDroot_symmetry-asa.Zbreak*tan(asa.sweep_angleLE2)))/2)*(asa.Zbreak);
asa.swing = (asa.semiarea1 + asa.semiarea2)*2; %[m^2] area da asa.
asa.swing_ref = (asa.semiarea1 + asa.semiarea3)*2; %[m^2] area de referencia da asa (porçao trapezoidal da asa projetada no plano de simetria da aeronave).


% O calculo da corda media aerodinamica se da atraves da integracao de c^2 em relaçao a z de 0 a b/2 (semi-envergadura).
asa.integral_1 = 1/3*asa.Zbreak*asa.XLEbreak^2-2/3*asa.Zbreak*asa.XLEbreak*asa.XLEroot_symmetry+1/3*asa.Zbreak*asa.XLEroot_symmetry^2-asa.CHORDroot_symmetry*asa.Zbreak*asa.XLEbreak+asa.CHORDroot_symmetry*asa.Zbreak*asa.XLEroot_symmetry+asa.CHORDroot_symmetry^2*asa.Zbreak;
asa.integral_2 = 1/3*asa.semi_span*asa.CHORDtip^2+1/3*asa.CHORDbreak^2*asa.semi_span+1/3*asa.CHORDbreak*asa.semi_span*asa.CHORDtip-1/3*asa.Zbreak*asa.CHORDtip^2-1/3*asa.CHORDbreak^2*asa.Zbreak-1/3*asa.Zbreak*asa.CHORDtip*asa.CHORDbreak;

asa.cma = (2/asa.swing)*(asa.integral_1 + asa.integral_2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OTIMIZAÇÃO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=[0.02639;0.00208;-0.02255;0.35611;0.05000;-0.00724;-0.19506;0.74977;0.01130;0.00000;-1.42458;  0.01544;0.00239;-0.00100;0.25000;0.05000;-0.01879;-0.08037;0.70000;0.01385;0.00000;-1.68577;   0.01870;0.00300;-0.00100;0.37500;0.05583;-0.01574;-0.11546;0.70000;0.02622;0.00302;-4.13687];
x2=[0.02243;0.00124;-0.00902;0.26801;0.05073;-0.00968;-0.19130;0.81224;0.01195;-0.00020;-1.24855;0.01343;0.00249;-0.03912;0.31900;0.05016;-0.01083;-0.05942;0.72849;0.01184;0.00491;0.30124;0.01561;0.00234;-0.04734;0.30435;0.05464;-0.01710;-0.07684;0.71405;0.01209;0.00066;-2.88228];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTADOS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Root airfoil
r0_11 = x1(1); 
Hte_11 = x1(2);
phi_11= x1(3);
X_tcmax_11 = x1(4);
t_c_11 = x1(5);
theta_11 = x1(6);
epsilon_11 = x1(7);
X_Ycmax_11 =  x1(8);
Ycmax_11 =  x1(9);
Yt_11 = x1(10);
epsil11 = x1(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break airfoil
r0_21 = x1(12);
Hte_21 = x1(13);
phi_21= x1(14);
X_tcmax_21 = x1(15);
t_c_21 = x1(16);
theta_21 = x1(17);
epsilon_21 = x1(18);
X_Ycmax_21 = x1(19);
Ycmax_21 = x1(20);
Yt_21 = x1(21);
epsil21 = x1(22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tip airfoil
r0_31 = x1(23);
Hte_31 = x1(24);
phi_31 = x1(25);
X_tcmax_31 = x1(26);
t_c_31 = x1(27);
theta_31 = x1(28);
epsilon_31 = x1(29);
X_Ycmax_31 = x1(30);
Ycmax_31 = x1(31);
Yt_31 = x1(32);
epsil31 = x1(33);
% Root airfoil
r0_12 = x2(1); 
Hte_12 = x2(2);
phi_12= x2(3);
X_tcmax_12 = x2(4);
t_c_12 = x2(5);
theta_12 = x2(6);
epsilon_12 = x2(7);
X_Ycmax_12 =  x2(8);
Ycmax_12 =  x2(9);
Yt_12 = x2(10);
epsil12 = x2(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break airfoil
r0_22 = x2(12);
Hte_22 = x2(13);
phi_22= x2(14);
X_tcmax_22 = x2(15);
t_c_22 = x2(16);
theta_22 = x2(17);
epsilon_22 = x2(18);
X_Ycmax_22 = x2(19);
Ycmax_22 = x2(20);
Yt_22 = x2(21);
epsil22 = x2(22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tip airfoil
r0_32 = x2(23);
Hte_32 = x2(24);
phi_32 = x2(25);
X_tcmax_32 = x2(26);
t_c_32 = x2(27);
theta_32 = x2(28);
epsilon_32 = x2(29);
X_Ycmax_32 = x2(30);
Ycmax_32 = x2(31);
Yt_32 = x2(32);
epsil32 = x2(33);

[a11, b11]=trans(r0_11, Hte_11, phi_11, X_tcmax_11, t_c_11, theta_11, epsilon_11, X_Ycmax_11, Ycmax_11, Yt_11);
[a21, b21]=trans(r0_21, Hte_21, phi_21, X_tcmax_21, t_c_21, theta_21, epsilon_21, X_Ycmax_21, Ycmax_21, Yt_21);
[a31, b31]=trans(r0_31, Hte_31, phi_31, X_tcmax_31, t_c_31, theta_31, epsilon_31, X_Ycmax_31, Ycmax_31, Yt_31);


[xsing11, ysing11, trail11, slopt11, radius11, x11, yu11, yl11, t_c11, Hte11] = perfil(a11,b11,40); % o modelo de perfil vai mudar mais para frente.
[xsing21, ysing21, trail21, slopt21, radius21, x21, yu21, yl21, t_c21, Hte21] = perfil(a21,b21,40); % o modelo de perfil vai mudar mais para frente.
[xsing31, ysing31, trail31, slopt31, radius31, x31, yu31, yl31, t_c31, Hte31] = perfil(a31,b31,40); % o modelo de perfil vai mudar mais para frente.

[a12, b12]=trans(r0_12, Hte_12, phi_12, X_tcmax_12, t_c_12, theta_12, epsilon_12, X_Ycmax_12, Ycmax_12, Yt_12);
[a22, b22]=trans(r0_22, Hte_22, phi_22, X_tcmax_22, t_c_22, theta_22, epsilon_22, X_Ycmax_22, Ycmax_22, Yt_22);
[a32, b32]=trans(r0_32, Hte_32, phi_32, X_tcmax_32, t_c_32, theta_32, epsilon_32, X_Ycmax_32, Ycmax_32, Yt_32);


[xsing12, ysing12, trail12, slopt12, radius12, x12, yu12, yl12, t_c12, Hte12] = perfil(a12,b12,40); % o modelo de perfil vai mudar mais para frente.
[xsing22, ysing22, trail22, slopt22, radius22, x22, yu22, yl22, t_c22, Hte22] = perfil(a22,b22,40); % o modelo de perfil vai mudar mais para frente.
[xsing32, ysing32, trail32, slopt32, radius32, x32, yu32, yl32, t_c32, Hte32] = perfil(a32,b32,40); % o modelo de perfil vai mudar mais para frente.

figure;h1=plot(horzcat(x11,fliplr(x11)),horzcat(yu11,fliplr(yl11)),'k',horzcat(x12,fliplr(x12)),horzcat(yu12,fliplr(yl12)),'--r');xlim([0 1]);axis equal;ylim([-0.1 0.1]);ylabel('Espessura');xlabel('Corda');title('Perfil da Raiz');legend('PSO','GA','Location','NorthEastOutside');
figure;h2=plot(horzcat(x21,fliplr(x21)),horzcat(yu21,fliplr(yl21)),'k',horzcat(x12,fliplr(x12)),horzcat(yu12,fliplr(yl12)),'--r');xlim([0 1]);axis equal;ylim([-0.1 0.1]);ylabel('Espessura');xlabel('Corda');title('Perfil da Quebra');legend('PSO','GA','Location','NorthEastOutside');
figure;h3=plot(horzcat(x31,fliplr(x31)),horzcat(yu31,fliplr(yl31)),'k',horzcat(x12,fliplr(x12)),horzcat(yu12,fliplr(yl12)),'--r');xlim([0 1]);axis equal;ylim([-0.1 0.1]);ylabel('Espessura');xlabel('Corda');title('Perfil da Ponta');legend('PSO','GA','Location','NorthEastOutside');
