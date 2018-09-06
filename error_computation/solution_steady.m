function [ u_exacte ] = solution_steady( x )
%{
Here are some classical solurion to steady problem -Delta(u)=f
%}
%global Y
%Solution continuous
    %u_exacte=(1/pi^2)*cos(pi*x);
    u_exacte=sin(x);
    %u_exacte=atan(Y*(x-pi));
    %--------------------------------
    %Solution discontinuous
%     global S
%     %x1=5;
%     %x2=11;
%     x3=2.5;
%     %x4=30;
%     %x5=35;
%     u_exacte=2/pi*((atan(S*(x-x3))));
    %u_exacte=2/pi*((atan(S*(x-x1)))-(atan(S*(x-x2)))+(atan(S*(x-x3)))-(atan(S*(x-x4)))+(atan(S*(x-x5))));
%--------------------------------
%Solution of the discontinuous source
% %     global h_values
% %     global x_h_values
% %     Taille1=max(size(x_h_values));
% %     for i=1:1:Taille1
% %            if abs(x_h_values(i)-x)<10^-6
% %               u_exacte=h_values(i);
% %            end
% %     end
%--------------------------------
end

