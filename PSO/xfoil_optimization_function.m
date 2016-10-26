function [runflag, CLmaxw]=xfoil_optimization_function(a1, a2, a3, b1, b2, b3, CLvectorA, CLvectorB, CLLA, CLLB)

global asa ntr mach alt number

[xsing1, ysing1, trail1, slopt1, radius1, x1, yu1, yl1, t_c1, Hte1, yc1, t1] = perfil(a1,b1,number.xfoil);
[xsing2, ysing2, trail2, slopt2, radius2, x2, yu2, yl2, t_c2, Hte2, yc2, t2] = perfil(a2,b2,number.xfoil);
[xsing3, ysing3, trail3, slopt3, radius3, x3, yu3, yl3, t_c3, Hte3, yc3, t3] = perfil(a3,b3,number.xfoil);
disp('Partindo para o xfoil');
 
[dens,temp,vsom,visc]=atmosfera(alt.estol);
V = mach.estol*vsom;
rey.estol(1) = (dens*V*asa.CHORDroot) / visc;
rey.estol(2) = (dens*V*asa.CHORDbreak) / visc;
rey.estol(3) = (dens*V*asa.CHORDtip) / visc;

[CLmax(1),alfamax(1),npontos(1)] = xfoil_create_file(mach.estol, rey.estol(1), x1, yu1, yl1);
[CLmax(2),alfamax(2),npontos(2)] = xfoil_create_file(mach.estol, rey.estol(2), x1, yu2, yl2);
[CLmax(3),alfamax(3),npontos(3)] = xfoil_create_file(mach.estol, rey.estol(3), x1, yu3, yl3);

runflag=0;
CLmaxw=0;
if any(CLmax==0) || any(npontos<2)
    runflag=0;
    CLmaxw=0
else
    dCLlocaldCL=(CLvectorB(1,:)-CLvectorA(1,:))/(CLLB-CLLA);
    Z=CLvectorB(2,:);
    [dummy,pontosenvergadura]=size(Z);
    for i=1:pontosenvergadura
        if Z(i)<=asa.Zbreak
        CLm(i)=((CLmax(2)-CLmax(1))/(asa.Zbreak-asa.Zroot_symmetry))*Z(i)+CLmax(1);
        else
        CLm(i)=((CLmax(3)-CLmax(2))/(asa.Ztip-asa.Zbreak))*(Z(i)-asa.Zbreak)+CLmax(2);
        end
    end
    for j=1:300
        CL=0.01*(j-1);
        for i=1:pontosenvergadura
            CLlocal(j,i)=dCLlocaldCL(i)*(CL-CLLB)+CLvectorB(1,i);
        end
    end
    passou=zeros(300,pontosenvergadura);
    for j=1:300
        for i=1:pontosenvergadura
            if CLlocal(j,i)>=CLm(i)
                passou(j,i)=1;
            end
        end
    end
    for j=1:300
        for i=1:pontosenvergadura
            if passou(j,i)==1
                CLmaxw=0.01*(j-1)
                runflag=1;
                break
            end
        end
        if CLmaxw~=0
            break
        end
    end
end