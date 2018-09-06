function [ theory ] = test_theory( Np,c_correction,kappa,hk,factor,ksi,choice_results )
%{
calculates \tau_theoric for IP 
%}
kappa=10^8;% The method does not depend on \kappa; we take the value which minimize \tau*
[~,Lagrange]=derivative_Lagrange(ksi,Np);
l=zeros(Np,1);
r=zeros(Np,1);
for k=1:1:Np
    l(k)=polyval(Lagrange(k,:),-1);
    r(k)=polyval(Lagrange(k,:),1);
end
[ g_left,~,~,~ ] = corrections_functions( ksi,Np,c_correction,kappa );
g1=g_left*r;
g_m_1=g_left*l;
if strcmp(choice_results,'IP')
    theory=1/2*2/hk*(abs(g1)-g_m_1);
elseif strcmp(choice_results,'BR2')
    %theory=1/(2*factor)*2/hk*(abs(g_left(Np))-g_left(1));%For equidistant
    %or LGL
    theory=1/(2*factor)*2/hk*(abs(g1)-g_m_1);
else
    theory=1/2*2/hk*(abs(g_left(Np))-g_left(1));
end

