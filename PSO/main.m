function main

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

%%%%LIMITES%%%%%
r0_min      = [0.01 0.01 0.01]'; 
r0_max      = [0.03 0.03 0.03]';

Hte_min     = [0 0 0]';
Hte_max     = [0.003 0.003 0.003]';

phi_min     = [-0.21 -0.33 -0.43]';
phi_max     = [-0.001 -0.001 -0.001]';

X_tcmax_min = [0.25 0.25 0.25]';
X_tcmax_max = [0.5 0.5 0.5]';

t_c_min     = [0.05 0.05 0.05]';
t_c_max     = [0.10 0.10 0.10]';

theta_min   = [-0.03 -0.03 -0.03]';
theta_max   = [0 0 0]';

epsilon_min = [-0.3 -0.3 -0.3]';
epsilon_max = [-0.005 -0.005 -0.005]';

X_Ycmax_min = [0.7 0.7 0.7]';
X_Ycmax_max = [0.9 0.9 0.9]';

Ycmax_min   = [0.01 0.01 0.01]';
Ycmax_max   = [0.03 0.03 0.03]';

Yt_min      = [-0.005 -0.005 -0.005]';
Yt_max      = [0.005 0.005 0.005]';

epsil_min   = [-5 -5 -5]';
epsil_max   = [5 5 5]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%OPÇÕES DE OTMIZAÇÃO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dims=33;                                                       %                                                            %
min  = [r0_min(1); Hte_min(1); phi_min(1); X_tcmax_min(1); t_c_min(1); theta_min(1); epsilon_min(1); X_Ycmax_min(1); Ycmax_min(1); Yt_min(1);epsil_min(1);r0_min(2); Hte_min(2); phi_min(2); X_tcmax_min(2); t_c_min(2); theta_min(2); epsilon_min(2); X_Ycmax_min(2); Ycmax_min(2); Yt_min(2);epsil_min(2);r0_min(3); Hte_min(3); phi_min(3); X_tcmax_min(3); t_c_min(3); theta_min(3); epsilon_min(3); X_Ycmax_min(3); Ycmax_min(3); Yt_min(3);epsil_min(3)];%POSIÇÕES MÍNIMAS  dims x 1 vector
max  = [r0_max(1); Hte_max(1); phi_max(1); X_tcmax_max(1); t_c_max(1); theta_max(1); epsilon_max(1); X_Ycmax_max(1); Ycmax_max(1); Yt_max(1);epsil_max(1);r0_max(2); Hte_max(2); phi_max(2); X_tcmax_max(2); t_c_max(2); theta_max(2); epsilon_max(2); X_Ycmax_max(2); Ycmax_max(2); Yt_max(2);epsil_max(2);r0_max(3); Hte_max(3); phi_max(3); X_tcmax_max(3); t_c_max(3); theta_max(3); epsilon_max(3); X_Ycmax_max(3); Ycmax_max(3); Yt_max(3);epsil_max(3)];%POSIÇÕES MÁXIMAS  dims x 1 vector
PSOseedValue = (min+max)'/2; 

mv            = abs((max-min)/2); %VELOCIDADE MÁXIMA DOS COMPONENTES
VarRange      = [min max];        %RANGE DAS VARIÁVEIS
minmax        = 0;                % =0 MINIMIZAR  =1 MAXIMIZAR  =2 NÃO SEI
PSOparams(1)  = 35;               %MÁXIMO NÚMERO DE EPOCHS
PSOparams(2)  = 100;              %TAMANHO DA POPULAÇÃO
PSOparams(3)  = 1;                %CONSTANTE DE ACELERAÇÃO 1 - INFLUÊNCIA DO MÍNIMO LOCAL
PSOparams(4)  = 3;                %CONSTANTE DE ACELERAÇÃO 2 - INFLUÊNCIA DO MÍNIMO GLOBAL
PSOparams(5)  = 0.9;              %CONSTANTE DE INÉRCIA INICIAL
PSOparams(6)  = 0.4;              %CONSTANTE DE INÉRCIA FINAL 
PSOparams(7)  = 50;               %ITERAÇÃO EM QUE A CONSTANTE DE INÉRCIA CHEGA AO VALOR FINAL
PSOparams(8)  = 1e-8;             %VELOCIDADE DE CONVERGÊNCIA DO CRITÉRIO DE PARADA POR GRADIENTE DE ERRO
PSOparams(9)  = 35;               %NÚMERO DE ITERAÇÕES COM ERRO ABAIXO DO GRADIENTE PARA PARAR A ITERAÇÃO
PSOparams(10) = NaN;              %META DE ERRO
PSOparams(11) = 0;                %==1 PARA DETERMINAR AS POSIÇÕES INICIAIS DAS PARTÍCULAS MANUALMENTE, ==0 RANDOMICAMENTE 

functname='optimization_function';
PSOseedValue = repmat(PSOseedValue, PSOparams(2),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OTIMIZAÇÃO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,te,tr]=pso_filipe(functname,dims,mv,VarRange,minmax,PSOparams,PSOseedValue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTADOS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Root airfoil
r0_1 = x(1); 
Hte_1 = x(2);
phi_1= x(3);
X_tcmax_1 = x(4);
t_c_1 = x(5);
theta_1 = x(6);
epsilon_1 = x(7);
X_Ycmax_1 =  x(8);
Ycmax_1 =  x(9);
Yt_1 = x(10);
epsil1 = x(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break airfoil
r0_2 = x(12);
Hte_2 = x(13);
phi_2= x(14);
X_tcmax_2 = x(15);
t_c_2 = x(16);
theta_2 = x(17);
epsilon_2 = x(18);
X_Ycmax_2 = x(19);
Ycmax_2 = x(20);
Yt_2 = x(21);
epsil2 = x(22);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tip airfoil
r0_3 = x(23);
Hte_3 = x(24);
phi_3 = x(25);
X_tcmax_3 = x(26);
t_c_3 = x(27);
theta_3 = x(28);
epsilon_3 = x(29);
X_Ycmax_3 = x(30);
Ycmax_3 = x(31);
Yt_3 = x(32);
epsil3 = x(33);

[a1, b1]=trans(r0_1, Hte_1, phi_1, X_tcmax_1, t_c_1, theta_1, epsilon_1, X_Ycmax_1, Ycmax_1, Yt_1);
[a2, b2]=trans(r0_2, Hte_2, phi_2, X_tcmax_2, t_c_2, theta_2, epsilon_2, X_Ycmax_2, Ycmax_2, Yt_2);
[a3, b3]=trans(r0_3, Hte_3, phi_3, X_tcmax_3, t_c_3, theta_3, epsilon_3, X_Ycmax_3, Ycmax_3, Yt_3);


[xsing1, ysing1, trail1, slopt1, radius1, x1, yu1, yl1, t_c1, Hte1] = perfil(a1,b1,40); % o modelo de perfil vai mudar mais para frente.
[xsing2, ysing2, trail2, slopt2, radius2, x2, yu2, yl2, t_c2, Hte2] = perfil(a2,b2,40); % o modelo de perfil vai mudar mais para frente.
[xsing3, ysing3, trail3, slopt3, radius3, x3, yu3, yl3, t_c3, Hte3] = perfil(a3,b3,40); % o modelo de perfil vai mudar mais para frente.

figure;plot(x1,yu1,'k',x1,yl1,'k');axis equal;grid on;ylabel('Espessura');xlabel('Corda');
figure;plot(x2,yu2,'k',x2,yl2,'k');axis equal;grid on;ylabel('Espessura');xlabel('Corda');
figure;plot(x3,yu3,'k',x3,yl3,'k');axis equal;grid on;ylabel('Espessura');xlabel('Corda');

[dens,temp,vsom,visc]=atmosfera(alt.cru);
V = mach.cru1*vsom;
rey.cru1 = (dens*V*asa.cma) / visc;

[dens,temp,vsom,visc]=atmosfera(alt.cru);
V = mach.cru2*vsom;
rey.cru2 = (dens*V*asa.cma) / visc;

alfa=1; 
CLL=0.45;
al_cl=1;
[time1, CD1, CL1, CLvector1, CDvector1] = create_file(mach.cru1, rey.cru1, xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLL, al_cl);
disp('Copie os arquivos do BLWF que quiser.')
pause;
alfa=1; 
CLL=0.35;
al_cl=1;
[time2, CD2, CL2, CLvector2, CDvector2] = create_file(mach.cru2, rey.cru2, xsing1, ysing1, trail1, slopt1, x1, yu1, yl1, epsil1, xsing2, ysing2, trail2, slopt2, yu2, yl2, epsil2, xsing3, ysing3, trail3, slopt3, yu3, yl3, epsil3, alfa, CLL, al_cl);
disp('Copie os arquivos do BLWF que quiser.')
pause;

disp('L/D mínimo(Mach 1):')
disp(-CL1/CD1);
disp('L/D mínimo(Mach 2):')
disp(-CL2/CD2);
disp('Perfil ótimo:');
disp(x);
disp('Número de iterações:');
disp(te);
disp('Trajetória do mínimo global:');
disp(tr);