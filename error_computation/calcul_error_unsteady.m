function [ error ] = calcul_error_unsteady( U_final,Np,Nbr_Elements,M,x,t )
%{
    This function calculates the error
%}
%{
error=0;
for i=1:1:Nbr_Elements
    error=error+(((U_final(1+(i-1)*Np:Np+(i-1)*Np)-solution_unsteady(x(1+(i-1)*Np:Np+(i-1)*Np),t)')')...
       *M*(U_final(1+(i-1)*Np:Np+(i-1)*Np)-solution_unsteady(x(1+(i-1)*Np:Np+(i-1)*Np),t)'));
end
disp('M')
error=sqrt(error)

pause
%}

%{
%GL points
global xini
global xfin
error=0;
Np2=Np;
[ksi]=lgwt(Np,-1,1);
ksi=ksi';
ksi=sort(ksi);


[ksi2,W2]=lgwt(Np2,-1,1);
ksi2=ksi2';
ksi2=sort(ksi2);
[~,Lagrange]=derivative_Lagrange(ksi,Np);

%}

%GLL points
global xini
global xfin

error=0;
Np2=Np+4;
ksi=lglnodes(Np-1);
[ksi2,W2]=lglnodes(Np2-1);
[~,Lagrange]=derivative_Lagrange(ksi,Np);

%}

coord_mesh=linspace(xini,xfin,Nbr_Elements+1);
x2=zeros(1,Np2*(Nbr_Elements));
for i=1:1:Nbr_Elements
    x2(1+(i-1)*Np2:i*Np2)=coordonnees_elem(coord_mesh,i,ksi2,Np2);
end
hk=(xfin-xini)/Nbr_Elements;

Eval_sol=zeros(Np2,Np);
for k=1:1:Np
    for m=1:1:Np2
        Eval_sol(m,k)=polyval(Lagrange(k,:),ksi2(m));
    end
end


W2=diag(W2);
for i=1:1:Nbr_Elements
    u=Eval_sol*U_final(1+(i-1)*Np:Np+(i-1)*Np);
%     size(u)
%     pause();
%     size(solution_unsteady(x2(1+(i-1)*Np2:Np2+(i-1)*Np2),t))
%     pause();
    error=error+(((u'-solution_unsteady(x2(1+(i-1)*Np2:Np2+(i-1)*Np2),t)))...
       *W2*(u-solution_unsteady(x2(1+(i-1)*Np2:Np2+(i-1)*Np2),t)'));
end
disp('quadrature\n')
fprintf('Nbr Elements=%d\n',Nbr_Elements);
error=sqrt(hk/2*sum(error));

end

