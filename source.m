function [ S ] = source( Np,Nbr_Elements,x,t )
global mu
global Y
S=zeros(Np*Nbr_Elements,1);
%S=-exp(-x)*exp(-t)*(mu+1);
%S=-exp(-t)*(atan(Y*(x-pi))-2*Y^3*mu*(x-pi)/(1+(Y*(x-pi)).^2).^2);
S=0;
end

