function [ c ] = factor_BR2_IP( ksi,Np,invM,hk )
%{
We compute the factor enabling to switch from a BR2 scheme to a IP scheme
%}
%Method1
[~,Lagrange]=derivative_Lagrange(ksi,Np);
l=zeros(Np,1);
r=zeros(Np,1);
for k=1:1:Np
    l(k)=polyval(Lagrange(k,:),-1);
    r(k)=polyval(Lagrange(k,:),1);
end
LR=zeros(1,2*Np);
LR(1,1:Np)=r;
LR(1,Np+1:2*Np)=l;
c=1/4*LR*blkdiag(invM,invM)*LR';

%Method 2 from the article
%c=1/(2/Np^2*hk);
end

