function [sol]=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end)
global xini
global xfin
%-------------------------------------------
%Classical RK4
%-------------------------------------------
%{
k1_laplacian=-RHS3( u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k1_laplacian(1:2*Np,1)=k1_laplacian(1:2*Np,1)+solution_unsteady(xini,t-dt)*B_ini;
k1_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k1_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_unsteady(xfin,t-dt)*B_end;

K1=u+dt/2*k1_laplacian;
k2_laplacian=-RHS3( K1,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k2_laplacian(1:2*Np,1)=k2_laplacian(1:2*Np,1)+K1(1)*B_ini;
k2_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k2_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+K1(Np*Nbr_Elements)*B_end;

K2=u+dt/2*k2_laplacian;
k3_laplacian=-RHS3( K2,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k3_laplacian(1:2*Np,1)=k3_laplacian(1:2*Np,1)+K2(1)*B_ini;
k3_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k3_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+K2(Np*Nbr_Elements)*B_end;

K3=u+dt*k3_laplacian;
k4_laplacian=-RHS3( K3,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );   
k4_laplacian(1:2*Np,1)=k4_laplacian(1:2*Np,1)+k3_laplacian(1)*B_ini;
k4_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k4_laplacian(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+K3(Np*Nbr_Elements)*B_end;

sol=u+dt/6*(k1_laplacian+2*k2_laplacian+2*k3_laplacian+k4_laplacian); %We solve du/dt=Delta(u)
sol(1)=solution_unsteady(xini,t);
sol(Np*Nbr_Elements)=solution_unsteady(xfin,t);

%}


%-------------------------------------------
%Castonguay's RK54 Carpenter and KEnnedy
%-------------------------------------------

k1=-RHS3( x,t,u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k1(1:2*Np,1)=k1(1:2*Np,1)+solution_unsteady(xini,t-dt)*B_ini;
k1(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k1(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_unsteady(xfin,t-dt)*B_end;

k2=-RHS3( x,t,k1,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k2(1:2*Np,1)=k2(1:2*Np,1)+k1(1)*B_ini;
k2(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k2(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+k1(Np*Nbr_Elements)*B_end;

k3=-RHS3( x,t,k2,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k3(1:2*Np,1)=k3(1:2*Np,1)+k2(1)*B_ini;
k3(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k3(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+k2(Np*Nbr_Elements)*B_end;
k4=-RHS3( x,t,k3,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k4(1:2*Np,1)=k4(1:2*Np,1)+k3(1)*B_ini;
k4(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k4(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+k3(Np*Nbr_Elements)*B_end;

k5=-RHS3( x,t,k4,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );
k5(1:2*Np,1)=k5(1:2*Np,1)+k4(1)*B_ini;
k5(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=k5(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+k4(Np*Nbr_Elements)*B_end;


sol=(u+dt*k1+dt^2*k2/2+dt^3*k3/6+dt^4*k4/24+dt^5*k5/200);

sol(1)=coeff_u_ini*[solution_unsteady(xini,t);sol(2:Np)];%Boundary conditions
sol(Np*Nbr_Elements)=coeff_u_end*[sol(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_unsteady(xfin,t)];%Boundary conditions





end