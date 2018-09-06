function [ gradient_lagrange,Lagrange_tot ] = derivative_Lagrange( ksi,Np )
%----------------------------------------------------------------------
%----------------------------------------------------------------------
%Purpose : Calculate the derivative of Lagrange function
%   L'1(xi1)    ...     L'Np(xi1)
%   L'1(xi2)
%   ...
%   L'1(xiN)    ...     L'Np(xiN)

%And calculate all the Lagrangian at the soluion points
%   L1_coeff_max     ...     L1_coeff_min
%   ...              ...     ...
%   LNp_coeff_max    ...     LNp_coeff_max
%----------------------------------------------------------------------
%----------------------------------------------------------------------
gradient_lagrange=zeros(Np,Np);
Lagrange_tot=zeros(Np,Np);
for k=1:1:Np
    Y=zeros(1,Np);
    Y(k)=1;
    Lagrange=lagrangepoly(ksi,Y);%lagrangepoly(x((j-1)*Np+1:j*Np),Y);%
    Lp=polyder(Lagrange);
    gradient_lagrange(:,k)=polyval(Lp,ksi);
    Lagrange_tot(k,:)=Lagrange;
end
end

