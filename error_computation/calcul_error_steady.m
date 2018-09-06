function [ error ] = calcul_error_steady( U_final,Np,Nbr_Elements,M,x )
%{
    This function calculates the error
%}

%{
error=0;
for i=1:1:Nbr_Elements
    error=error+(((U_final(1+(i-1)*Np:Np+(i-1)*Np)-solution_steady(x(1+(i-1)*Np:Np+(i-1)*Np))')')...
       *M*(U_final(1+(i-1)*Np:Np+(i-1)*Np)-solution_steady(x(1+(i-1)*Np:Np+(i-1)*Np))'));
end

error=sqrt(error);

%}


%{
%GL points
global xini
global xfin
error=0;
Np2=Np+5;
[ksi]=lgwt(Np,-1,1);
ksi=ksi';
ksi=sort(ksi);


[ksi2,W2]=lgwt(Np2,-1,1);
ksi2=ksi2';
ksi2=sort(ksi2);
%}



%GLL points
global xini
global xfin

error=0;
Np2=Np+5;
ksi=lglnodes(Np-1);
[ksi2,W2]=lglnodes(Np2-1);
%}



[~,Lagrange]=derivative_Lagrange(ksi,Np);



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
M2=creationM(ksi2,Np2,hk);


W2=diag(W2);
for i=1:1:Nbr_Elements
    u=Eval_sol*U_final(1+(i-1)*Np:Np+(i-1)*Np);
%     size(u)
%     pause();
%     size(solution_steady(x(1+(i-1)*Np2:Np2+(i-1)*Np2)))
%     pause();
    error=error+(((u'-solution_steady(x2(1+(i-1)*Np2:Np2+(i-1)*Np2))))...
       *W2*(u-solution_steady(x2(1+(i-1)*Np2:Np2+(i-1)*Np2))'));
end
error=sqrt(hk/2*error);
%}
end


