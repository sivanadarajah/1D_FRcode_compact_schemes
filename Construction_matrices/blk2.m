function [ R ] = blk2( A,B,Np )
%----------------------------------
%Purpose : create the local matices Delta :
%The stencil is of 3 if compact or 5 if wide
%The matrices get shifted. This program handles this shift/translation

%blk2 ; 2 for 2 matrices : beginning or ending
%----------------------------------
R=zeros(2*Np,3*Np);
lenghtA=max(size(A));
lenghtB=max(size(B));
if lenghtA>lenghtB%Since means we are at the end [Qint,Qend]
    R(1:Np,1:3*Np)=A;
    R(Np+1:2*Np,Np+1:3*Np)=B;
end
if lenghtA<lenghtB%Since means we are at the beginning [Qini,Qint]
    R(1:Np,1:2*Np)=A;
    R(Np+1:2*Np,1:3*Np)=B;
end
end

