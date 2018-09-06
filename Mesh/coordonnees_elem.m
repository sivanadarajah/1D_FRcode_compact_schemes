function a=coordonnees_elem(x,k,ksi,Np)
%{
We compute the coordinates of the elements
%}
hk=x(k+1)-x(k);
xinik=x(k)*ones(1,Np);
a=xinik+(1+ksi)/2*hk;
end

