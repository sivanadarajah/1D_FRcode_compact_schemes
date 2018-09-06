function [ R ] = blk3( A,B,C,Np )
%----------------------------------
%Purpose : create the local matices Delta :
%The stencil is of 3 if compact or 5 if wide
%The matrices get shifted. This program handles this shift/translation

%blk3 ; 3 for 3 matrices : ini+1 ; int ; end-1
%----------------------------------
lenghtA=max(size(A));
lenghtC=max(size(C));
if floor(lenghtA/Np)==2 || floor(lenghtC/Np)==2
    R=zeros(3*Np,4*Np);
else
    R=zeros(3*Np,5*Np);
end

R(1:Np,1:lenghtA)=A;
%------
if floor(lenghtA/Np)==2
    R(Np+1:2*Np,1:3*Np)=B;
    R(2*Np+1:3*Np,4*Np-lenghtC+1:4*Np)=C;
elseif floor(lenghtC/Np)==2
    R(Np+1:2*Np,Np+1:4*Np)=B;
    R(2*Np+1:3*Np,4*Np-lenghtC+1:4*Np)=C;
else
    R(Np+1:2*Np,Np+1:4*Np)=B;
    R(2*Np+1:3*Np,5*Np-lenghtC+1:5*Np)=C;
end

end

