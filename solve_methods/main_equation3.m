function [ laplacian ] = main_equation3( u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend )
%{
%We calculate the laplacian --> the main equation for explit solver 
%In fact it may be Delta(u)-c*u

A=-Delta(u)-c*u!!!!!!!!!!!!!!!!!


%}       
laplacian=zeros(Np*Nbr_Elements,1);
lenghtint=max(size(Aint));
if floor(lenghtint/Np)==5 %wide stencil
    k=2;
elseif floor(lenghtint/Np)==3 %compact stencil
    k=1;
else
    error('Problem with the matrix A');
end
for i=1:1:Nbr_Elements
    row1=1+(i-1)*Np;
    row2=Np+row1-1;
    switch i
        case 1
            laplacian(row1:row2,1)=Aini(:,2:(k+1)*Np)*u(row1+1:row2+k*Np);
        case 2
            laplacian(row1:row2,1)=Aini2(:,2:(k+2)*Np)*u(row1-Np+1:row2+k*Np);
        case Nbr_Elements-1
            laplacian(row1:row2,1)=Aend2(:,1:(k+2)*Np-1)*u(row1-k*Np:row2+Np-1);
        case Nbr_Elements
            laplacian(row1:row2,1)=Aend(:,1:(k+1)*Np-1)*u(row1-k*Np:row2-1);
        otherwise
            laplacian(row1:row2,1)=Aint*u(row1-k*Np:row2+k*Np);
    end
end

end

