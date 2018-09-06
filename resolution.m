function [ erreur,tau_result,s_result,t_result ] = resolution( Np,ksi,Nbr_Elements,hk,M,invM,beta,tau,s,c_correction,kappa,...
        choice_analysis,choice_results,time_dependant,method_resolution,resolution_explicit,choice_scheme )
%{
PURPOSE : Calculates the error of the solution for any problem
Main structure
We create the local matrices that will be used for all cases --> Uniformity
of all cases
    I-Von Neumann analysis
        1)Stability for tau --> removed for this GitHub
        2)Plot graphs --> removed for this GitHub
        3)Find the maximum time-steps
    II-Classical test case
        1)Implicit
            -->function construction_globale
        2)Explicit
            a)Time independant
                i) Euler method
                    --> function euler_explicit steady
                    --> function residual
                    --> function RHS
                ii)RK4 method
                    --> function RK4 steady
                    --> function residual
                    --> function RHS        
            b)Time dependant
                i)Error analysis
                    A)dt constant
                        i) Euler method unsteady
                            --> function euler_explicit unsteady
                            --> function RHS
                        ii)RK4 method unsteady
                            --> function RK4 unsteady
                            --> function RHS
                    B)Look for dt_max
                        i) Euler method unsteady
                            --> function euler_explicit unsteady
                            --> function RHS
                        ii)RK4 method unsteady
                            --> function RK4 unsteady
                            --> function RHS
                ii)Limit for tau
                    A) Euler method unsteady
                        --> function euler_explicit unsteady
                        --> function RHS
                    B)RK4 method unsteady
                        --> function RK4 unsteady
                        --> function RHS    
                iii)Limit for BR2
                    A) Euler method unsteady
                        --> function euler_explicit unsteady
                        --> function RHS
                    B)RK4 method unsteady
                        --> function RK4 unsteady
                        --> function RHS 
%}
global xini
global xfin
global mu
global c
global Y
%Initialization of all outputs to 0.
erreur=0;
tau_result=0;
s_result=0;
t_result=0;
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%% I-Von-Neumann analysis
if strcmp(choice_analysis,'VN')
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%   3)Dt maximum
    if strcmp(choice_results,'time')
        %We construct the Laplacian
        [~,~,~,~,Aint,~,~,~,~]=local_matrices3(Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis);
        dtini=10;%High value to be unstable
        D_dt=10;
        for i=1:1:7
            dt=dtini;
            D_dt=D_dt/10;%
            test_stable=0;
            number_iterations=1;
            while number_iterations<12 && test_stable~=1
                [lambda_max]=VN_dt_maximum( Np,Aint,dt,choice_scheme );
                if lambda_max>=1 %lambda_max>=1
                    test_stable=0;
                    fprintf('Unstable for dt=%f ; lambda=%f \n',dt,lambda_max);
                    dt=dt-D_dt;
                else
                    test_stable=1;
                    fprintf('--------------------\n');
                    fprintf('Stable for dt=%f\n',dt);
                    fprintf('--------------------\n');
                end
                number_iterations=number_iterations+1;
            end
            if test_stable==1
                dtini=dt+D_dt;
            else
                dtini=dtini/10;
                fprintf('\n Problem : We do not have the correct order of dt yet\n \n');
            end
        end
        t_result=dt;
    end
end
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%% II-Test case analysis
%We retrieve the boundary conditions, the matrices with the Laplacian
[coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end]=local_matrices3(Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis);
if strcmp(choice_analysis,'TC')
    coord_mesh=linspace(xini,xfin,Nbr_Elements+1);%coordinates of the mesh
    x=zeros(1,Np*(Nbr_Elements));
    for i=1:1:Nbr_Elements
        x(1+(i-1)*Np:i*Np)=coordonnees_elem(coord_mesh,i,ksi,Np);
    end
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
    %1)Implicit
    if strcmp(method_resolution,'Implicit')%Steady problem
        if strcmp(choice_results,'IP')%We find the minimum value of tau which gives a compact scheme.
            s=0;
            dtau=100;
            tau=-1000;
            for i=1:1:3
                dtau=dtau/10;
                test_final_stability=1;
                while test_final_stability==1
                    %We retrieve the boundary conditions, the matrices with the Laplacian
                    [coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end]=local_matrices3(Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis);
                    %We retrieve the solution
                    u=construction_globale3(ksi,Np,Nbr_Elements,x,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);%Solveur
                    erreur=calcul_error_steady(u',Np,Nbr_Elements,M,x)
                    pause();
                    if abs(erreur)>10^-2
                        fprintf('Unstable tau =%f\n',tau);
                        tau=tau+dtau;
                    else
                        test_final_stability=0;
                    end
                end
                fprintf('--------------------\n');
                fprintf('Stable for tau=%f\n',tau);
                fprintf('--------------------\n');
                tau=tau-dtau;
            end
            tau_result=tau;%We retrieve the value of tau which gives a stable method
        end
        if strcmp(choice_results,'error')
            %We immediatly compute the error
            u=construction_globale3(ksi,Np,Nbr_Elements,x,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
            erreur=calcul_error_steady(u',Np,Nbr_Elements,M,x);
        end
    end
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
    %2)Explicit
    if strcmp(method_resolution,'Explicit')%We use a explicit method : time step is required
        %Time parameter for diffusion
        xmin = min(abs(x(1)-x(2)));
        CFL=0.5;
        dt = CFL/(mu)*xmin^2;
        dt = .5*dt;%Time step suggested for a diffusive problem Hesthaven and Warburton
        if c~=0
            %For the time step convective
            xmin = min(abs(x(1)-x(2)));
            CFL=0.75;
            dt_adv = CFL/(c)*xmin;
            dt_adv = .5*dt_adv;
            dt=min(dt,dt_adv);
        end
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
        %2.a)Time independant
        if strcmp(time_dependant,'No')%Steady solution
            u=zeros(Np*Nbr_Elements,1);%Guess of the solution
            u(1)=coeff_u_ini*[solution_steady(xini);u(2:Np)];
            u(Np*Nbr_Elements)=coeff_u_end*[u(Np*(Nbr_Elements-1)+1:Np*Nbr_Elements-1);solution_steady(xfin)];
            res=1;
            Nite=0;
            %-----------------------------------------------
            %2.a.i)Euler explicit
            if strcmp(resolution_explicit,'Euler_explicit')            
                while res>10^-10
                    %We compute the solution with the Euler algorithm
                    [u,res]=Euler_explicit_steady3(u,Np,Nbr_Elements,x,dt,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                    res
                    Nite=Nite+1;
    %                 clf;
    %                 axis([xini xfin -1.5 1.5])
    %                 hold on
    %                 plot(x,u,'+-','LineWidth',1.2)
    %                 hold on
    %                 %plot(x,u0(x-dt*i*c))
    %                 plot(x,solution_exacte(x));
    %                 title(sprintf('NDG at t = sec Equidistant points, Np=%d',Np))
    %                 xlabel('x')
    %                 ylabel('solution')
    %                 legend('NDG','Exact solution')
    %                 pause(0.0000001)
                    if max(abs(u))>1000
                        error('Unstable');
                    else
                    end
                end
                plot(x,u,'+-')
                hold on
                plot(x,solution_steady(x),'*-')
                erreur=calcul_error_steady(u,Np,Nbr_Elements,M,x);
                erreur
            end
            %-----------------------------------------------
            %2.a.ii)RK4
            if strcmp(resolution_explicit,'RK4')            
                while res>10^-10
                    %We compute the solution with the Runge Kutta algorithm
                    [u,res]=RK4_steady3(u,Np,Nbr_Elements,x,dt,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                    res
                    Nite=Nite+1;
                    %If we want to see the evolution of the solution :
                    %uncomment the following lines
    %                 clf;
    %                 axis([xini xfin -1.5 1.5])
    %                 hold on
    %                 plot(x,u,'+-','LineWidth',1.2)
    %                 hold on
    %                 %plot(x,u0(x-dt*i*c))
    %                 plot(x,solution_exacte(x));
    %                 title(sprintf('NDG at t = sec Equidistant points, Np=%d',Np))
    %                 xlabel('x')
    %                 ylabel('solution')
    %                 legend('NDG','Exact solution')
    %                 pause(0.0000001)
                    if max(abs(u))>2
                        error('Unstable');
                    else
                    end
                end
                erreur=calcul_error_steady(u,Np,Nbr_Elements,M,x);
                erreur
            end
        end
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%---------------------------------------------------------------------------------------------
        %2.b)Time dependant
        if strcmp(time_dependant,'Yes')%Unsteady problems
            tini=0;%Initial time
            tfin=2;%Ending time
            Nt = ceil((tfin-tini)/dt);
            dt = (tfin-tini)/Nt;


    %If we want to study the rate of convergence of a advection problem
            %u0=5*1/(3*sqrt(2*pi))*exp(-(x-20).^2/(2*3^2));%
    %If we want to study the rate of convergence of a diffusive problem
            %u0=exp(-x.^2/(4*mu*tini));
    %If we want to study the case of a advection diffusion problem
            u0=sin(x)+cos(x);
            %u0=exp(-x);
            %u0=atan(Y*(x-pi));

%---------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
%           i)Result for the error
            if strcmp(choice_results,'error')%We want to compute the error of the method
                choice_temps='dt_max';%
                'dt_max'%We also want to find the maximal time step possible
                %'fix', we use the time step suggested by Hesthaven and
                %Warburton
                if strcmp(choice_temps,'fix')%We use a constant time step
                    %-----------------------------------------------
                    %2.b.i.i)Euler explicit
                    if strcmp(resolution_explicit,'Euler_explicit')
                        u=u0';
%                         figure('units','normalized','outerposition',[0 0 1 1]),
%                         set(gca, 'Position', [0.03 0.03 0.94 0.94])
%                         hold on
                        for i=1:1:Nt+1
                            t=tini+i*dt;
                            u=Euler_explicit_unsteady3(u,Np,Nbr_Elements,dt,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
%                             clf;
%                             axis([xini xfin -1.5 1.5])
%                             hold on
%                             plot(x,u,'+-','LineWidth',1.2)
%                             hold on
%                             %plot(x,u0(x-dt*i*c))
%                             plot(x,solution_unsteady(x,t));
%                             title(sprintf('NDG at t = %d sec Equidistant points, Np=%d',t,Np))
%                             xlabel('x')
%                             ylabel('solution')
%                             legend('NDG','Exact solution')
%                             pause(0.001)
                            if max(abs(u))>3
                                error('unstable');
                            else
                            end
                        end
                        erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
                        erreur
                    end
                    %-----------------------------------------------
                    %2.b.i.ii)RK4
                    if strcmp(resolution_explicit,'RK4')
                        u=u0';
%                         figure('units','normalized','outerposition',[0 0 1 1]),
%                         set(gca, 'Position', [0.03 0.03 0.94 0.94])
%                         hold on
                        for i=1:1:Nt
                            t=tini+i*dt;
                            u=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
%                             clf;
%                             axis([xini xfin -2 2])
%                             hold on
%                             plot(x,u,'+-','LineWidth',1.2)
%                             hold on
%                             %plot(x,u0(x-dt*i*c))
%                             plot(x,solution_unsteady(x,t),'+-');
%                             title(sprintf('NDG at t = %d sec Equidistant points, Np=%d',t,Np))
%                             xlabel('x')
%                             ylabel('solution')
%                             legend('NDG','Exact solution')
%                             pause(0.000001)
                            if max(abs(u))>3
                                fprintf('Unstable for tau=%f \n',tau);
                                return
                            else
                            end
                        end
                        erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
                        erreur
                    end
                end
%------------------------------------------------------------------------------------------                
                if strcmp(choice_temps,'dt_max')%We try to find the maximal step
                    dt_ini=1;
                    D_dt=1;
                    nbr_chiffre=0;
%------------------------------------------------------------------------------------------
                    %2.b.i.i)Euler explicit
                    if strcmp(resolution_explicit,'Euler_explicit')
                        while nbr_chiffre<2
                            if nbr_chiffre==0
                                limite_ite=11;
                            else
                                limite_ite=12;
                            end
                            dt=dt_ini;
                            D_dt=D_dt/10;
                            test_stabilite=0;
                            number_iterations=1;
                            while number_iterations<limite_ite && test_stabilite==0
                                t=tini;
                                test_stable=0;
                                u=u0';
                                while t<tfin && test_stable==0
                                    t=t+dt;
                                    u=Euler_explicit_unsteady3(u,Np,Nbr_Elements,dt,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                                    if max(abs(u))>2
                                        test_stable=1;
                                        fprintf('Unstable for dt=%f ; u=%d \n',dt,max(abs(u)));
                                    else
                                    end
                                end
                                if test_stable==1
                                    dt=dt-D_dt;
                                    %pause();
                                    %fprintf('\n Problem : We do not have the correct order of dt yet\n \n');
                                    number_iterations=number_iterations+1;
                                    t=tini;
                                else
                                    dt_ini=dt+D_dt;
                                    test_stabilite=1;
                                end
                            end
                            if test_stabilite==0
                                dt_ini=dt_ini/10;
                                %pause();
                                fprintf('\n Problem : We do not have the correct order of dt yet\n \n');
                            else
                                fprintf('Stable for dt=%d',dt);
                                dt_ini=dt+D_dt;
%                                 pause();
                                nbr_chiffre=nbr_chiffre+1;
                            end    
                        end
                        erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
                        erreur
                        t_result=dt
                        pause();
                    end
%------------------------------------------------------------------------------------------
                    %2.b.i.ii)RK4
                    if strcmp(resolution_explicit,'RK4')
                        while nbr_chiffre<2
                            if nbr_chiffre==0
                                limite_ite=11;
                            else
                                limite_ite=12;
                            end
                            dt=dt_ini;
                            D_dt=D_dt/10;
                            test_stabilite=0;
                            number_iterations=1;
                            while number_iterations<limite_ite && test_stabilite==0
                                t=tini;
                                test_stable=0;
                                u=u0';
                                while t<tfin && test_stable==0
                                    t=t+dt;
                                    u=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                                    if max(abs(u))>2
                                        test_stable=1;
                                        fprintf('Unstable for dt=%f ; u=%d \n',dt,max(abs(u)));
                                    else
                                    end
                                end
                                if test_stable==1
                                    dt=dt-D_dt;
                                    %pause();
                                    %fprintf('\n Problem : We do not have the correct order of dt yet\n \n');
                                    number_iterations=number_iterations+1;
                                    t=tini;
                                else
                                    dt_ini=dt+D_dt;
                                    test_stabilite=1;
                                end
                            end
                            if test_stabilite==0
                                dt_ini=dt_ini/10;
                                %pause();
                                fprintf('\n Problem : We do not have the correct order of dt yet\n \n');
                            else
                                fprintf('Stable for dt=%d\n',dt);
                                dt_ini=dt+D_dt;
%                                 pause();
                                nbr_chiffre=nbr_chiffre+1;
                            end    
                        end
%                         erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
%                         erreur
                        t_result=dt;
                    end
                    %We compute the error for tfin=1 with the maximal time
                    %step found previously
                    tini=0;
                    tfin=1;%100 100 is used for purely advective PDE
                    Nt = ceil((tfin-tini)/dt); %dt has proved to yield a stable scheme; we re-launch the simulation with this dt to compute the error.
                    dt = (tfin-tini)/Nt;

                    if strcmp(resolution_explicit,'RK4')
                        u=u0';
%                         figure('units','normalized','outerposition',[0 0 1 1]),
%                         set(gca, 'Position', [0.03 0.03 0.94 0.94])
%                         hold on
                        for i=1:1:Nt
                            t=tini+i*dt;
                            u=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
%                             clf;
%                             axis([xini xfin -1.5 1.5])
%                             hold on
%                             plot(x,u,'+-','LineWidth',1.2)
%                             hold on
%                             %plot(x,u0(x-dt*i*c))
%                             plot(x,solution_unsteady(x,t));
%                             title(sprintf('NDG at t = %d sec Equidistant points, Np=%d',t,Np))
%                             xlabel('x')
%                             ylabel('solution')
%                             legend('NDG','Exact solution')
%                             pause(0.000001)
                            if max(abs(u))>2
                                fprintf('Unstable for tau=%f \n',tau);
                                return
                            else
                            end
                        end
                        erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
                        erreur
                    end
                end      
            end           
%------------------------------------------------------------------------------------------            
%-------------------------------------------------------------------------
%           i)Result for the tau ; we compute the minimum value of tau
%           yielding a stable method
            if strcmp(choice_results,'IP')
                s=0;
                disp('1')
                dtau=10;
                tau=-100;
                for i=1:1:3
                    dtau=dtau/10;
                    test_final_stability=1;
                    while test_final_stability==1
                        [coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end]=local_matrices3(Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis);
                        test_stable=0;
                        iterations_temps=1;
                        u=u0';              
                        while iterations_temps<Nt && test_stable~=1
                            t=tini+iterations_temps*dt;
                            if strcmp(resolution_explicit,'RK4')
                                    u=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                            end
                            if strcmp(resolution_explicit,'Euler_explicit')
                                    u=Euler_explicit_unsteady3(u,Np,Nbr_Elements,dt,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                            end
                            if max(abs(u))>1.5
                                test_stable=1;
                                fprintf('Unstable tau =%f\n',tau);
                            else
                            end
                            iterations_temps=iterations_temps+1;
                        end
                        if test_stable==0
                            test_final_stability=0;
                            fprintf('--------------------\n');
                            fprintf('Stable for tau=%f\n',tau);
                            fprintf('--------------------\n');
                        else
                            tau=tau+dtau;
                        end
                    end
                    tau=tau-dtau
                end
                tau_result=tau+dtau;
                erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
            end
%------------------------------------------------------------------------------------------            
%-----------------------------------------------------------------------------------------------------------------------------------
%Results for BR2
            if strcmp(choice_results,'BR2')%We find the minimal value of s, the penalty term of the BR2 scheme, ensuring a stable method
                tau=0;
                disp('BR2')
                ds=10;
                s=-100;
                for i=1:1:3
                    ds=ds/10;
                    test_final_stability=1;
                    while test_final_stability==1
                        [coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end]=local_matrices3(Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis); 
                        test_stable=0;
                        iterations_temps=1;
                        u=u0';              
                        while iterations_temps<Nt && test_stable~=1
                            t=tini+iterations_temps*dt;
                            if strcmp(resolution_explicit,'RK4')
                                    u=RK4_unsteady3(u,Np,Nbr_Elements,dt,x,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                            end
                            if strcmp(resolution_explicit,'Euler_explicit')
                                    u=Euler_explicit_unsteady3(u,Np,Nbr_Elements,dt,t,coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end);
                            end
                            if max(abs(u))>1.5
                                test_stable=1;
                                fprintf('Unstable s=%f\n',s);
                            else
                            end
                            iterations_temps=iterations_temps+1;
                        end
                        if test_stable==0
                            test_final_stability=0;
                            fprintf('--------------------\n');
                            fprintf('Stable for s=%f\n',s);
                            fprintf('--------------------\n');
                        else
                            s=s+ds;
                        end
                    end
                    s=s-ds;
                end
                s_result=s+ds;
                erreur=calcul_error_unsteady(u,Np,Nbr_Elements,M,x,t);
            end
        end
    end
%------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------
end

