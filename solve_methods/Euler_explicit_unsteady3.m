function [sol]=Euler_explicit_unsteady3(u,Np,Nbr_Elements,dt,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end)
%{
       This function solves using the explicit method : Euler_explicit
%}
global xini
global xfin

rhs=-RHS3( u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );

rhs(1:2*Np,1)=rhs(1:2*Np,1)+solution_unsteady(xini,t)*B_ini;
rhs(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=rhs(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_unsteady(xfin,t)*B_end;

sol=u+dt*rhs;
sol(1)=coeff_u_ini*[solution_unsteady(xini,t);sol(2:Np)];
sol(Np*Nbr_Elements)=coeff_u_end*[sol(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_unsteady(xfin,t)];

end