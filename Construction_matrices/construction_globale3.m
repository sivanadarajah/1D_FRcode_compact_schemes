function [ U_final ] = construction_globale3(ksi,Np,Nbr_Elements,x,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end)
%{
    This function creates the global matrix A and then solves A*U=B
    It is only used for implicit problem
%}
global xini
global xfin
[~,Lagrange]=derivative_Lagrange(ksi,Np);
l=zeros(Np,1);
r=zeros(Np,1);
for k=1:1:Np
    l(k)=polyval(Lagrange(k,:),-1);
    r(k)=polyval(Lagrange(k,:),1);
end
%% 
A=sparse(Np*Nbr_Elements,Np*Nbr_Elements);
lenghtint=max(size(Aint));
if floor(lenghtint/Np)==5 %wide stencil
    k=2;
elseif floor(lenghtint/Np)==3 %compact stencil
    k=1;
else
    error('Problem with the matrix A');
end
F=sparse(Np*Nbr_Elements,1);
for i=1:1:Nbr_Elements*Np
    F(i)=f(x(i));
end
for i=1:1:Nbr_Elements
    row1=1+(i-1)*Np;
    row2=Np+row1-1;
    switch i
        case 1
            A(row1:row2,row1:row2+k*Np)=Aini;
        case 2
            A(row1:row2,row1-1*Np:row2+k*Np)=Aini2;
        case Nbr_Elements-1
            A(row1:row2,row1-k*Np:row2+1*Np)=Aend2;
        case Nbr_Elements
            A(row1:row2,row1-k*Np:row2)=Aend;
        otherwise
            A(row1:row2,row1-k*Np:row2+k*Np)=Aint;
    end
end

%We have built the big matrix A. Now we impose the boundary conditions
%Boundary conditions
%{
%------------------------------------------
B=coeff_u_ini(1)*A(:,1)+coeff_u_end(Np)*A(:,Np*Nbr_Elements);

for i=1:1:Np
    if i==1
        A(:,Np*(Nbr_Elements-1)+i)=A(:,Np*(Nbr_Elements-1)+i)+A(:,Np*Nbr_Elements)*coeff_u_end(i);
    elseif i==Np
        A(:,i)=A(:,i)+A(:,1)*coeff_u_ini(i);
    else
        A(:,i)=A(:,i)+A(:,1)*coeff_u_ini(i);
        A(:,Np*(Nbr_Elements-1)+i)=A(:,Np*(Nbr_Elements-1)+i)+A(:,Np*Nbr_Elements)*coeff_u_end(i);
    end
end
F_final=F-B';
F_final=F_final(2:Nbr_Elements*Np-1);
F_final=F_final';
%}
F(1:2*Np,1)=F(1:2*Np,1)+solution_steady(xini)*B_ini;
F(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=F(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_steady(xfin)*B_end;

F_final=zeros(Np*Nbr_Elements-2,1);
F_final(:,1)=F(2:Np*Nbr_Elements-1,1);
A_final=A(2:Np*Nbr_Elements-1,2:Np*Nbr_Elements-1);
A_final=sparse(A_final);
U=A_final\F_final;

U_final=zeros(1,Np*Nbr_Elements);
U_final(1)=coeff_u_ini*[solution_steady(xini);U(1:Np-1)];
U_final(Np*Nbr_Elements)=coeff_u_end*[U(Np*(Nbr_Elements-1)-2+2:Np*Nbr_Elements-2);solution_steady(xfin)];
U_final(2:Np*Nbr_Elements-1)=U;
% cond(A)
% pause();
% l'*U_final(1:Np)'-solution_steady(xini)
% pause();
% r'*U_final(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements)'-solution_steady(xfin)
% pause();

%----------------------------------------------------------


