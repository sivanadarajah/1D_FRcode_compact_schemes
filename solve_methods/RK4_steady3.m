function [sol,res]=RK4_steady3(u,Np,Nbr_Elements,x,dt,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end)
%{
       This function solves using the explicit method : Euler_explicit
%}
global xini
global xfin
sol=zeros(Np*Nbr_Elements,1);   
res=zeros(Np*Nbr_Elements,1); 


res1=residual3(u,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
res2=residual3(u+dt/2*res1,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
res3=residual3(u+dt/2*res2,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
res4=residual3(u+dt*res3,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
    
sol=u+dt/6*(res1+2*res2+2*res3+res4); %We solve du/dt=Delta(u)
sol(1)=coeff_u_ini*[solution_steady(xini);sol(2:Np)];
sol(Np*Nbr_Elements)=coeff_u_end*[sol(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_steady(xfin)];

res=max(residual3(sol,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end));



%{
res1=residual3(u,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
% res1(1:2*Np,1)=res1(1:2*Np,1)+solution_steady(xini)*B_ini;
% res1(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=res1(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_steady(xfin)*B_end;


res2=residual3(res1,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
% res2(1:2*Np,1)=res2(1:2*Np,1)+res1(1)*B_ini;
% res2(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=res2(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+res1(Np*Nbr_Elements)*B_end;

res3=residual3(res2,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
% res3(1:2*Np,1)=res3(1:2*Np,1)+res2(1)*B_ini;
% res3(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=res3(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+res2(Np*Nbr_Elements)*B_end;

res4=residual3(res3,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
% res4(1:2*Np,1)=res4(1:2*Np,1)+res3(1)*B_ini;
% res4(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=res4(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+res3(Np*Nbr_Elements)*B_end;

res5=residual3(res4,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
% res5(1:2*Np,1)=res5(1:2*Np,1)+res4(1)*B_ini;
% res5(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=res5(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+res4(Np*Nbr_Elements)*B_end;

sol=(u+dt*res1+dt^2*res2/2+dt^3*res3/6+dt^4*res4/24+dt^5*res5/200); %We solve du/dt=Delta(u)
sol(1)=coeff_u_ini*[solution_steady(xini);sol(2:Np)];
sol(Np*Nbr_Elements)=coeff_u_end*[sol(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_steady(xfin)];




res=max(residual3(sol,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end));
%}

end