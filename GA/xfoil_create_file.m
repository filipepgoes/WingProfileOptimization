function [CL,alfamax,npontos] = xfoil_create_file(mach, rey, x,yu,yl)



[coly,liny]=size(yu);
for i=1:liny
   yu_perm(i)=yu(end-i+1); 
end
[colx,linx]=size(x);
for i=1:linx
   x_perm(i)=x(end-i+1);
end
fid=fopen('profile.pr','w');
fprintf(fid, 'PERFIL\n');
fprintf(fid, '%10.5f%14.10f\n',[x_perm; yu_perm]);
fprintf(fid, '%10.5f%14.10f\n',[x(2:linx); yl(2:liny)]);
fclose(fid);

alfa.inf=5;
alfa.sup=15;
alfa.div=0.5;
alfa.number=(alfa.sup-alfa.inf)/alfa.div;

fid=fopen('Entrada.key','w');

fprintf(fid, 'load profile.pr\n');
fprintf(fid, 'ppar\n');
fprintf(fid, 't\n');
fprintf(fid, '0.45\n');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, 'oper\n');
fprintf(fid, 'visc\n');
fprintf(fid, '%f\n',rey);
fprintf(fid, 'mach\n');
fprintf(fid, '%f\n',mach);
fprintf(fid, 'pacc\n');
fprintf(fid, 'profile_results.txt\n');
fprintf(fid, 'profile_results.dmp\n');
fprintf(fid, 'aseq\n'); 
fprintf(fid, '%f\n',alfa.inf);
fprintf(fid, '%f\n',alfa.sup);
fprintf(fid, '%f\n',alfa.div);
fprintf(fid, 'quit');

fclose(fid);

if exist('profile_results.txt')==2
    !del profile_results.txt
    !del profile_results.dmp
end

!process.exe starter.exe 20 C:\TG\GA\starter.exe Entrada.key

if exist('profile_results.txt')==2
        dummy = textread('profile_results.txt','%q','whitespace','\n','headerlines',12);
        [alfacalc,dummy2]=size(dummy);
        for i=1:alfacalc
            nextline=textread('profile_results.txt','%c',7,'whitespace','\n','headerlines',11+i);
            if any(nextline=='?')
                break
            end
            line=textread('profile_results.txt','%f',7,'whitespace','\n','headerlines',11+i);
            alfavector(i)=line(1);
            CLvector(i)=line(2);
        end
        if alfacalc~=0 & alfacalc~=1
            [CL,idx]=max(CLvector);
            alfamax=alfavector(idx);
            npontos=alfacalc;
        else
            CL=0;
            alfamax=0;
            npontos=0;
        end
else
    CL=0;
    alfamax=0;
    npontos=0;
end