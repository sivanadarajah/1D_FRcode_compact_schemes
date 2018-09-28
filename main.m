clear all
%close all
format long
addpath(genpath('Construction_matrices'))
addpath(genpath('error_computation'))
addpath(genpath('Mesh'))
addpath(genpath('Penalty_term'))
addpath(genpath('solve_methods'))
addpath(genpath('VN_methods'))
%------------------------------------------------
%SET YOUR PROBLEM : 
%Modify also file    
%    1)f.m (source term)
%    2)solution_steady/solution_unsteady
%------------------------------------------------

%------------------------------------------------
%{
This programm is adapted for diffusion problem
We can add an advective term for time dependent problem
%}
%------------------------------------------------




%Problem
global xini % Starting point of the 1D mesh
global xfin % Ending point of the 1D mesh
global c %velocity for the advective term
global alpha % Parameter controling the numerical flux for the advection term
global mu% Parameter controling the diffusion
global Y %Parameter used for tangential law for the second member : not necessarily used see f.m
Y=0.1;
xini=0;
xfin=2*pi;
c=0;%----Parameter controling the advection only present in time dependant problem
alpha=0;%Parameter controling the numerical flux for the advection term : (u*)={{u}}+(1-alpha)/2[|u|]
mu=1;%----Parameter controling the diffusion

pathfile=pwd;

%-----------------------------------------------
choice_analysis='TC';%'VN' (stands for Von-Neumann analysis) 'TC' (stands for classical test case)%_______%----Choose---------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
if strcmp(choice_analysis,'TC')
    choice_results='error';%'IP' we look for which value of tau yield a stable method ; 
    %'BR2' similar than 'IP' but for the BR2 scheme
    %'error' we run a classical test case knowing the solution and compute the L2-error%_______%----Choose---------------------
    %-----------------------------------------------
    %time dependency & Method resolution--> Relevant for TC
    time_dependant='Yes'; % If 'Yes' -> parabolic ; if 'No' -> Elliptic%_%-------Choose-------------
    method_resolution='Explicit';%'Explicit,'Implicit'%__________________%-------Choose-------------
    resolution_explicit='RK4';%'Euler_explicit','RK4'%________%-------Choose-------------
    %---------------------------------------------------------------------
    if strcmp(time_dependant,'Yes') && strcmp(method_resolution,'Implicit')
       error('system parabolic associated with a implicit method'); 
    end
    if strcmp(time_dependant,'No') && c~=0
       error('Advective term in a time-independant problem'); 
    end
    titre_error=strcat('\FINAL_error_IP_');
    if strcmp(choice_results,'error')
        %Parameters
        Np_array=[3,4];%number of points__________________________________________________%-------Choose-------------
        N_p=max(size(Np_array));
        for i=1:1:N_p %Loop over the polynomial degree
            Np=Np_array(i);
            titre_final=strcat(titre_error,['_p' num2str(Np-1) '.mat']);
            if Np<5
                Nbr_Elements_array=[32,64,128];%number of elements_______________________%-------Choose-------------
            else
                Nbr_Elements_array=[32];
            end
            hk_array=(xfin-xini)./Nbr_Elements_array;
            N_Elem=max(size(Nbr_Elements_array));
            N_c=4;%size of c_correctin: parameter used for the correction function of the main equation
            N_kappa=2;%size of kappa: parameter used for the correction function of the auxiliary equation
            N5=2;%size of tau
            ratio_tau=[1,1.5];%For choice_results='error', we compute \tau^{*} and we can then mutliply it byt this ratio
            result_erreur=zeros(N_c*N_kappa+2,N5*(N_Elem+(N_Elem-1)+1)+2);%Table of results
            L2_erreur=zeros(N_c*N_kappa,N5*(N_Elem));
            dt_max_array=zeros(N_c*N_kappa,N5*(N_Elem));
        %--------------------------------------------------------------------    
            choice_solution_points='LGL';%'Equidistant','LGL','GL';%_%-------Choose-------------
        %--------------------------------------------------------------------
            if strcmp(choice_solution_points,'Equidistant')
               ksi=linspace(-1,1,Np); 
            end
            if strcmp(choice_solution_points,'LGL')
               [ksi,W]=lglnodes(Np-1);
            end
            if strcmp(choice_solution_points,'GL')
               ksi=lgwt(Np,-1,1);
               ksi=ksi';
               ksi=sort(ksi);
            end
            %Corrections functions        
                p=Np-1;
                ap=factorial(2*p)/(2^p*(factorial(p)^2));
                cDG=0;
                cSD=2*p/((2*p+1)*(p+1)*(ap*factorial(p))^2);
                cHU=2*(p+1)/((2*p+1)*p*(ap*factorial(p))^2);
                switch Np
                    case 2
                        c_plus=8.706;%Value obtained after a choice_analysis='VN' analysis
                    case 3
                        c_plus=0.186;
                    case 4
                        c_plus=3.67*10^-3;
                    case 5
                        c_plus=4.79*10^-5;
                    case 6
                        c_plus=4.24*10^-7;
                end
                c_array=[cDG,cSD,cHU,c_plus];
                %extension_c=10.^linspace(-8,8,20);
                %c_array=[c_array,extension_c];
                kappa_array=[cDG,c_plus];
                %extension_kappa=10.^linspace(-8,8,20);
                %kappa_array=[kappa_array,extension_kappa];
            for j=1:1:N_Elem %Loop over the number of elements
                Nbr_Elements=Nbr_Elements_array(j);%Nbr_Elements_array(1);%
                hk=hk_array(j);%hk_array(1);%
                M=creationM(ksi,Np,hk);
                invM=inv(M);
                for k=1:1:N_c %Loop over the number of main correction functions
                    c_correction=c_array(k);
                    for m=1:1:N_kappa %Loop over the number of auxiliary correction functions
                        kappa=kappa_array(m);%kappa_array(m);%kappa_array(j);%kappa;% Nbr_Elements=Nbr_Elements_array(j);%
            %--------------------------------------------------------------------        
                        choice_scheme='compact';%Scheme : IP, LDG,..%___________________%-------Choose-------------
            %-------------------------------------------------------------------- 
            %---------------------------------------------------------------------  
            %CHOICE of TAU & BETA        
                        beta=0;%________________________________________________%-------Choose----------------  
                        %_________________________________________________%-------Choose----------------
                        s=0;
                        if strcmp(choice_scheme,'compact') && beta~=0
                            error('Choice of Compact scheme but \beta different from 0');
                        end
                        if strcmp(choice_scheme,'wide') && s~=0
                            error('Choice of Wide scheme with a BR2 penalty term');
                        end

                        factor=factor_BR2_IP( ksi,Np,invM,hk );
                        %theory(j,i,k)=test_theory( Np,c_correction,kappa,hk,factor,choice_solution_points,choice_results );
                        penalty_term=test_theory( Np,c_correction,kappa,hk,factor,ksi,choice_results );
                        %tau=theory(j,i,k);
                        %tau_values=linspace(0,100,N5);
                        N5r=max(size(ratio_tau));
                        if N5r~=N5
                            error('prelocatted value of size of tau erronous');
                        end
                        %tau=0;
    %                     if tau~=0 && s~=0
    %                         error('BR2 and IP influence');
    %                     end
                        %factor_array(i,j)=factor;
                        for n=1:1:N5
                            %theory=test_theory( Np,c_correction,kappa,hk,factor,ksi,choice_results );
                            %tau=0;
                            if strcmp(choice_scheme,'compact')
                                tau=ratio_tau(n)*penalty_term;
                            else
                                tau=ratio_tau(n);
                            end
                                [L2_erreur(m+N_kappa*(k-1),j+N_Elem*(n-1)),~,~,dt_max_array(m+N_kappa*(k-1),j+N_Elem*(n-1))]=resolution(Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
                                    choice_analysis,choice_results,time_dependant,method_resolution,resolution_explicit,choice_scheme);
                                L2_erreur
                                
                                %error_array
                            result_erreur(1,3+(n-1)*(N_Elem+(N_Elem-1)+1))=ratio_tau(n);
                            result_erreur(2,3+(N_Elem+(N_Elem-1)+1)*(n-1):3+(N_Elem+(N_Elem-1)+1)*(n-1)+N_Elem-1)=Nbr_Elements_array;
                        end
                    end
                     result_erreur(1+(k-1)*N_kappa+2:k*N_kappa+2,2)=kappa_array;
                    result_erreur((k-1)*N_kappa+3,1)=c_correction;
                end
            end
%----------------------------------------------------------------------------------------------------------------------------------------
        %Results
        order_convergence=zeros(N_c*N_kappa,N5*(N_Elem-1));
        for n=1:1:N5
            for k=1:1:N_c
                for m=1:1:N_kappa
                    for j=1:1:N_Elem-1
                        order_convergence(m+N_kappa*(k-1),j+(N_Elem-1)*(n-1))=log10(L2_erreur(m+N_kappa*(k-1),j+1+N_Elem*(n-1))/L2_erreur(m+N_kappa*(k-1),j+N_Elem*(n-1)))/log10(hk_array(j+1)/hk_array(j));
                    end
                end
            end
        end
        for n=1:1:N5
            result_erreur(3:2+N_kappa*N_c,3+(N_Elem+(N_Elem-1)+1)*(n-1):3+(N_Elem+(N_Elem-1)+1)*(n-1)+N_Elem-1)=L2_erreur(:,1+(n-1)*N_Elem:N_Elem*n);
            result_erreur(3:2+N_kappa*N_c,3+N_Elem+(N_Elem+(N_Elem-1)+1)*(n-1):3+N_Elem+(N_Elem+(N_Elem-1)+1)*(n-1)+N_Elem-1-1)=order_convergence(:,1+(n-1)*(N_Elem-1):1+(N_Elem-1)-1+(n-1)*(N_Elem-1));
            if Np<5
                result_erreur(3:2+N_kappa*N_c,2+(N_Elem+(N_Elem-1)+1)*n)=dt_max_array(:,1+(n-1)*N_Elem);
            else
                result_erreur(3:2+N_kappa*N_c,2+(N_Elem+(N_Elem-1)+1)*n)=dt_max_array(:,3+(n-1)*N_Elem);
            end
        end
        titre_final=strcat(pathfile,titre_final);
        save(titre_final,'result_erreur');
        end
    end

    
  %FIND TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL 
  %TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL TAU NUMERICAL 
    if strcmp(choice_results,'IP')|| strcmp(choice_results,'BR2')
        %Parameters
        Np_array=[3,4];%Np=3;%________________________________________________________%-------Choose-------------
        N_p=max(size(Np_array));%size(kappa_array);%
        N_c=5;
        N_kappa=1;
        error_array=zeros(N_c,N_kappa*N_p);
        tau_array=zeros(N_c,N_kappa*N_p);
        s_array=zeros(N_c,N_kappa*N_p);
        penalty_term=zeros(N_c,N_kappa*N_p);
        result_tau_array=zeros(2+N_c,1+N_kappa*N_p);
        for i=1:1:N_p
            Np=Np_array(i);
            Nbr_Elements=32;
            hk=(xfin-xini)/Nbr_Elements;
        %--------------------------------------------------------------------    
            choice_solution_points='LGL';%'Equidistant','LGL','GL';%_%-------Choose-------------
        %--------------------------------------------------------------------
            if strcmp(choice_solution_points,'Equidistant')
               ksi=linspace(-1,1,Np); 
            end
            if strcmp(choice_solution_points,'LGL')
               [ksi,W]=lglnodes(Np-1);
            end
            if strcmp(choice_solution_points,'GL')
               ksi=lgwt(Np,-1,1);
               ksi=ksi';
               ksi=sort(ksi);
            end
            %Corrections functions        
            p=Np-1;
            ap=factorial(2*p)/(2^p*(factorial(p)^2));
            cDG=0;
            cSD=2*p/((2*p+1)*(p+1)*(ap*factorial(p))^2);
            cHU=2*(p+1)/((2*p+1)*p*(ap*factorial(p))^2);
            switch Np
                case 2
                    c_plus=8.706;
                case 3
                    c_plus=0.186;
                case 4
                    c_plus=3.67*10^-3;
                case 5
                    c_plus=4.79*10^-5;
                case 6
                    c_plus=4.24*10^-7;
            end
            c_array=[cDG,cSD,cHU,c_plus,10^5];%1,10,100,1000];
            N_c_check=max(size(c_array));
            kappa_array=[cDG];%,cSD,cHU,c_plus,10^5];
            N_kappa_check=max(size(kappa_array));
            if N_c_check~=N_c
                error('prelocated size of c array false')
            end
            if N_kappa_check~=N_kappa
                error('prelocated size of kappa array false')
            end
            M=creationM(ksi,Np,hk);
            invM=inv(M);
            for k=1:1:N_c
                c_correction=c_array(k);
                for m=1:1:N_kappa
                    kappa=kappa_array(m);%kappa_array(j);%kappa;% Nbr_Elements=Nbr_Elements_array(j);%
        %--------------------------------------------------------------------        
                    choice_scheme='compact';%Scheme : IP, LDG,..%___________________%-------Choose-------------
        %-------------------------------------------------------------------- 
        %---------------------------------------------------------------------  
        %CHOICE of TAU & BETA        
                    beta=0;%________________________________________________%-------Choose----------------  
                    %_________________________________________________%-------Choose----------------
                    tau=0;
                    s=0;
%                     factor=factor_BR2_IP( ksi,Np,invM );
%                     theory(k,m+N_kappa*(i-1))=test_theory( Np,c_correction,kappa,hk,factor,choice_solution_points,choice_results );
                    if strcmp(choice_results,'IP')
                            [error_array(k,m),tau_array(k,m+N_kappa*(i-1)),~,~]=resolution(Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
                                                choice_analysis,choice_results,time_dependant,method_resolution,resolution_explicit,choice_scheme);
                            tau_array
                    end
                    if strcmp(choice_results,'BR2')
                            [error_array(k,m),tau_array(k,m),s_array(k,m),~]=resolution(Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
                                                choice_analysis,choice_results,time_dependant,method_resolution,resolution_explicit,choice_scheme);
                            s_array
                    end
                result_tau_array(2+k,1+m+N_kappa*(i-1))=tau_array(k,m+N_kappa*(i-1));
                end
            end
            result_tau_array(1,2+(i-1)*N_kappa)=Np;
            result_tau_array(2,2+(i-1)*N_kappa:2+i*N_kappa-1)=kappa_array;
        end
        result_tau_array(3:3+N_c-1,1)=c_array;
        save('Implicit_value_tau_numerical.mat','result_tau_array');
    end
end


%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------
% VON NEUMANN ANALYSIS VON NEUMANN ANALYSIS VON NEUMANN ANALYSIS
if strcmp(choice_analysis,'VN')
    Nbr_Elements=0;
    hk=1;
    choice_results='time';%through this simulation we can
    %1)retrieve the values of c+ : consider a pure advective equation
    %2)Find the maximal time step of the diffusion problem for instance
%-------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------------------------------------
% FIND THE MAXIMAL TIME STEP
    if strcmp(choice_results,'time')
        titre_ini='Final Von neumann analysis p=';
        Np_array=[3];%________________________________________________________%-------Choose-------------
        N_p=max(size(Np_array));
        N_c=2;%Elem
        N_kappa=24;%N_c
        N_tau=1;%N_kappa
        result_time=zeros(N_c*N_kappa+2,N_p*N_tau+2);
        dt_array=zeros(N_c*N_tau,N_p*N_tau);
        for i=1:1:N_p
            Np=Np_array(i);
            titre_huynh=strcat(titre_ini,sprintf('%d',Np-1));
        %--------------------------------------------------------------------    
            choice_solution_points='LGL';%'Equidistant','LGL','GL';%_%-------Choose-------------
        %--------------------------------------------------------------------
            if strcmp(choice_solution_points,'Equidistant')
               ksi=linspace(-1,1,Np); 
            end
            if strcmp(choice_solution_points,'LGL')
               ksi=lglnodes(Np-1); 
            end
            if strcmp(choice_solution_points,'GL')
               ksi=lgwt(Np,1,-1);
               ksi=ksi';
            end
            M=creationM(ksi,Np,hk);
            invM=inv(M);
            p=Np-1;
            ap=factorial(2*p)/(2^p*(factorial(p)^2));
            cDG=0;
            cSD=2*p/((2*p+1)*(p+1)*(ap*factorial(p))^2);
            cHU=2*(p+1)/((2*p+1)*p*(ap*factorial(p))^2);
            switch Np
                case 2
                    c_plus=8.706;
                case 3
                    c_plus=0.186;
                case 4
                    c_plus=3.67*10^-3;
                case 5
                    c_plus=4.79*10^-5;
                case 6
                    c_plus=4.24*10^-7;
            end
            c_array=[cDG,c_plus];%,cSD,cHU,c_plus];
            %c_array=sort(c_array);
            %c_array=sort(c_array);
            %c_array=[cDG,cSD,cHU,c_plus];
            kappa_array=[cDG,cSD,cHU,c_plus];
            extension_kappa=10.^linspace(-8,8,20);
            kappa_array=[kappa_array,extension_kappa];
            N2r=max(size(c_array));
            for j=1:1:N_c
                c_correction=c_array(j);% c_correction is the correction functions for the main equation
                for k=1:1:N_kappa
                    kappa=kappa_array(k);%kappa is the corrections functions for the auxiliary variable
                    factor=factor_BR2_IP( ksi,Np,invM );
                    tau_theory=test_theory( Np,c_correction,kappa,hk,factor,ksi,'IP' );
                    tau_array=[1];%,tau_theory*1.1,tau_theory*1.5];%
                    N4r=max(size(tau_array));
                    if N4r~=N_tau
                        error('size of N4 (tau) is not the same as the prelocated one');
                    end
                    for m=1:1:N_tau
%--------------------------------------------------------------------        
                        choice_scheme='wide';%Scheme : IP, LDG,..%___________________%-------Choose-------------
%--------------------------------------------------------------------  
                    %Compact scheme : WE WILL DO A SWIPE OF TAU BETA         
                        beta=-1/2;%________________________________________________%-------Choose----------------
                        s=0;
                        tau=tau_array(m);
                        if strcmp(choice_scheme,'compact') && beta~=0
                            error('Choice of Compact scheme but \beta different from 0');
                        end
%                         [~,~,~,dt_array((j-1)*N_kappa+k,(i-1)*N_tau+m)]=resolution_huynh_plot(Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
%                         choice_analysis,choice_results,0,0,0,choice_scheme);
                        [~,~,~,dt_array((j-1)*N_kappa+k,(i-1)*N_tau+m)]=resolution(Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
                        choice_analysis,choice_results,0,0,0,choice_scheme);
                    end
                    result_time(2,(i-1)*N_tau+3:1:i*N_tau+2)=tau_array;
                end
                result_time((j-1)*N_kappa+3:1:j*N_kappa+2,2)=kappa_array;
                result_time((j-1)*N_kappa+3,1)=c_correction;
            end
        result_time(1,(i-1)*N_tau+3)=Np;
        end
        result_time(3:N_c*N_kappa+2,3:N_p*N_tau+2)=dt_array;
        save(titre_huynh,'result_time');
    end
end