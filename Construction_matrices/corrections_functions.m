function [ g_left,g_right,h_left,h_right ] = corrections_functions( ksi,Np,c_correction,kappa )
%--------------------------------------------------------------------------
%This function calculates the DERIVATIVE vector of the corrections functions
% g is the correction functions for the auxiliary variables and depends on
% kappa

%h is the correction functions for the main equation and depends on
%c_correction

%As explained, the corrections functions are of degree Np
%If we apply the formula of the work of Castonguay : p->Np-1
%--------------------------------------------------------------------------
p=Np-1;
a=factorial(2*p)/(2^(p)*(factorial(p))^2);
eta=(2*p+1)*(a*factorial(p))^2/2;

eta_kappa=kappa*eta;
G_left=(-1)^p/2*([0;LegendrePoly(p)]-(eta_kappa*[0;0;LegendrePoly(p-1)]+LegendrePoly(p+1))/(1+eta_kappa));
G_right=1/2*([0;LegendrePoly(p)]+(eta_kappa*[0;0;LegendrePoly(p-1)]+LegendrePoly(p+1))/(1+eta_kappa));


g_left=polyder(G_left);
g_left=polyval(g_left',ksi);

g_right=polyder(G_right);
g_right=polyval(g_right',ksi);

eta_c=c_correction*eta;
H_left=(-1)^p/2*([0;LegendrePoly(p)]-(eta_c*[0;0;LegendrePoly(p-1)]+LegendrePoly(p+1))/(1+eta_c));
H_right=1/2*([0;LegendrePoly(p)]+(eta_c*[0;0;LegendrePoly(p-1)]+LegendrePoly(p+1))/(1+eta_c));
h_left=polyder(H_left);
h_left=polyval(h_left',ksi);
h_right=polyder(H_right);
h_right=polyval(h_right',ksi);
end

