function [sol,res]=Euler_explicit_steady3(u,Np,Nbr_Elements,x,dt,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end)
%{
       This function solves using the explicit method : Euler_explicit
%}
global xini
global xfin
sol=zeros(Np*Nbr_Elements,1);   
res=zeros(Np*Nbr_Elements,1); 
res=residual3(u,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);

sol(:,1)=u+dt*res;
sol(1)=coeff_u_ini*[solution_steady(xini);sol(2:Np)];
sol(Np*Nbr_Elements)=coeff_u_end*[sol(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_steady(xfin)];
res=max(residual3(sol,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end));
end