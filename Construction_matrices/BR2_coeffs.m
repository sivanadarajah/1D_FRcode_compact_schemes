function [ Re_ini,Re_int,Re_end ] = BR2_coeffs(Np,l,r,invM,Jump_ini,Jump,Jump_end)
%{
This program creates the matrix for the BR2 formulation
%}
LR=zeros(1,2*Np);
LR(1,1:Np)=r;
LR(1,Np+1:2*Np)=l;
C=1/4*LR*blkdiag(invM,invM)*LR';
Re_ini=C*(-Jump_ini);
Re_int=C*(-Jump);
Re_end=C*-(Jump_end);
%{
Re_ini=LR*invM*Jump_ini;
Re_int=LR*invM*Jump;
Re_end=LR*invM*Jump_end;
%}
%{
invBigM2=blkdiag(invM,invM);
invBigM3=blkdiag(invM,invM,invM);
Re_ini=LR*Jump_ini*invBigM2;
Re_int=LR*Jump*invBigM3;
Re_end=LR*Jump_end*invBigM2;
%}
end

