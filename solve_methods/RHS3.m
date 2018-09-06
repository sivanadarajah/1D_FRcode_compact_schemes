function [ rhs ] = RHS3( x,t,u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend )
global mu
global c
%{
This function  calculates the -Delta(u) !!!!! never Delta(u)
%}
Delta=main_equation3(u,Np,Nbr_Elements,Aini,Aini2,Aint,Aend2,Aend);
S=source(Np,Nbr_Elements,x,t);
rhs=Delta-S';%Delta includes also the advection term!!
end

