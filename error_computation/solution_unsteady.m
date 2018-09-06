function [ u_exacte ] = solution_unsteady( x,t )
%{
Here are some classical exact solution for unsteady problem : du/dt=Delta(u)
%}
global xini
global xfin
global c
global mu
global Y
tini=1;
%If we want to see the order  of convergence of a advective case
    %u_exacte=u0(x-c*t);
%If we want to see the order of convergence of a diffusive case
    %u_exacte=sqrt(tini/t)*exp(-x.^2/(4*mu*t));
%If we want to see a case of advection diffusion
    u_exacte=(sin(x-c*t)+cos(x-c*t))*exp(-mu*t);
    %u_exacte=exp(-x)*exp(-t);
    %u_exacte=atan(Y*(x-pi))*exp(-t);
end

