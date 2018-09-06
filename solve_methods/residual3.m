function [ residual ] = residual3( u,Np,Nbr_Elements,x,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end )
%{
    This function calculates the residual    
%}
global xini
global xfin
residual=zeros(Np*Nbr_Elements,1);

laplacian=RHS3( u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend );%laplacian also includes the advection term and is in fact -Delta
F=zeros(Np*Nbr_Elements,1);
for i=1:1:Nbr_Elements*Np
    F(i)=f(x(i));
end
F(1:2*Np,1)=F(1:2*Np,1)+solution_steady(xini)*B_ini;
F(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)=F(Np*(Nbr_Elements-2)+1:Np*Nbr_Elements,1)+solution_steady(xfin)*B_end;
% u
% save('vecteur','u')
% pause();
% laplacian
% pause();
residual=-(laplacian-F);%We solve  : -Delta(u)=f so the residual is Delta(u)+f (except that the laplacian is in fact -laplacian
% residual
% pause();
residual(1)=0;
residual(Np*Nbr_Elements)=0;

end

