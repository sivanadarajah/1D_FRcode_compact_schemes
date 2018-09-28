close all
clear all
filenameIP='FINAL_error_IP__p2';
filenameLDG='FINAL_error_LDG__p2';
resultIP=struct2array(load(strcat(filenameIP,'.mat')));
resultLDG=struct2array(load(strcat(filenameLDG,'.mat')));
[Nx,Ny]=size(resultLDG);
N2=3;%Number of element 
N5=2;%Number of parameter tau
N3=4;%Number of parameter c
N4=24;%Number of parameter kappa
new_resultIP=zeros(Nx-22*4,Ny-4+2*N5);
new_resultLDG=zeros(Nx,Ny+2*N5);
N2IP=3;
N2LDG=3;
for n=1:1:N5
    if n==1
        N2=N2-1;
        new_resultIP(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2)=resultIP(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2);
        N2=N2+1;
        new_resultLDG(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2)=resultLDG(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2);
    else
        N2=N2-1;
        new_resultIP(:,3+3+(N2+(N2-1)+1)*(n-1)-1:3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-3)=resultIP(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2);
        N2=N2+1;
        new_resultLDG(:,3+3+(N2+(N2-1)+1)*(n-1)-1:3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-3)=resultLDG(:,3+(N2+(N2-1)+1)*(n-1):3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-2);
    end
end

for n=1:1:N5%Loop over tau
    for k=1:1:N3%loop over c
            N2=N2IP;
            OOA_arrayIP=resultIP(2+1+2*(k-1),3+N2+(N2+(N2-1)+1)*(n-1):3+N2+(N2+(N2-1)+1)*(n-1)+N2-1-1);
            OOA_meanIP=mean(OOA_arrayIP);
            OOA_polyIP=polyfit(1:1:N2-1,OOA_arrayIP,0);
            R_2IP=sum((OOA_arrayIP-OOA_polyIP).^2)/sum((OOA_arrayIP-OOA_meanIP).^2);
        for m=1:1:N4%loop over kappa
            %For LDG
            N2=N2LDG;
            OOA_arrayLDG=resultLDG(2+m+N4*(k-1),3+N2+(N2+(N2-1)+1)*(n-1):3+N2+(N2+(N2-1)+1)*(n-1)+N2-1-1);
            OOA_meanLDG=mean(OOA_arrayLDG);
            OOA_polyLDG=polyfit(1:1:N2-1,OOA_arrayLDG,0);
            R_2LDG=sum((OOA_arrayLDG-OOA_polyLDG).^2)/sum((OOA_arrayLDG-OOA_meanLDG).^2);
            if n==1
                N2=N2IP;
                new_resultIP(2+m+2*(k-1),3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1)=OOA_polyIP;
                new_resultIP(2+m+2*(k-1),3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1+1)=R_2IP;
                %For LDG
                N2=N2LDG;
                new_resultLDG(2+m+N4*(k-1),3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1)=OOA_polyLDG;
                new_resultLDG(2+m+N4*(k-1),3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1+1)=R_2LDG;
            else
                N2=N2IP;
                new_resultIP(2+m+2*(k-1),3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-1)=OOA_polyIP;
                new_resultIP(2+m+2*(k-1),3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1)=R_2IP;
                %For LDG
                N2=N2LDG;
                new_resultLDG(2+m+N4*(k-1),3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1-1)=OOA_polyLDG;
                new_resultLDG(2+m+N4*(k-1),3+3+(N2+(N2-1)+1)*(n-1)+N2+(N2-1)+1)=R_2LDG;
            end
        end
    end
end

%We have calculated the "mean" OOA through the different meshes, now we
%will semilogx this OOA for the different 'c'



N2=3;
title_ini='\bf Order of convergence along the parameter $\bf \kappa$ for p=3 with c=';
for i=1:1:N3 %loop over parameter 'c'
    figure('Units','inches',...
'Position',[1 1 3 3],...
'PaperPositionMode','auto');
    Np=3;
    p=Np-1;
    ap=factorial(2*p)/(2^p*(factorial(p)^2));
    cDG=0;
    cSD=2*p/((2*p+1)*(p+1)*(ap*factorial(p))^2);
    cHU=2*(p+1)/((2*p+1)*p*(ap*factorial(p))^2);
    switch Np
        case 3
            c_plus=0.186;
        case 4
            c_plus=3.67*10^-3;
        case 5
            c_plus=4.79*10^-5;
        case 6
            c_plus=4.24*10^-7;
    end
%---------------------------------------------------------
    switch i
        case 1
            titre=strcat(title_ini,'\bf $ c_{DG}$');
        case 2
            titre=strcat(title_ini,' $\bf c_{SD}$');
        case 3
            titre=strcat(title_ini,' $\bf c_{HU}$');
        case 4
            titre=strcat(title_ini,' $\bf c_{+}$');
    end
%---------------------------------------------------------
        couleurLDG=[1 0 0];%red
        couleur=[0 1 0];%green
        couleurIP=[0 0 1];%blue
        for k=1:1:N5%Loop over tau
            m=0;
            fin=0;
            while fin==0
                m=m+1;
                switch m
                    case 1
                        kappa=cDG;
                        if k==1
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)-1);
                        else
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-3);
                        end
                        aLDG=plot(log10(kappa),OOALDG,'Marker','s','MarkerSize',6,'MarkerEdgeColor',couleurLDG,'MarkerFaceColor',[0 0 0]);
                        hold on
                    case 2
                        kappa=cSD;
                        if k==1
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)-1);
                        else
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-3);
                        end
                        bLDG=plot(log10(kappa),OOALDG,'Marker','v','MarkerSize',6,'MarkerEdgeColor',couleurLDG,'MarkerFaceColor',[0 0 0]);
                        hold on
                    case 3
                        kappa=cHU;
                        if k==1
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)-1);
                        else
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-3);
                        end
                        eLDG=plot(log10(kappa),OOALDG,'Marker','o','MarkerSize',6,'MarkerEdgeColor',couleurLDG,'MarkerFaceColor',[0 0 0]);
                        hold on
                    case 4
                        kappa=c_plus;
                        if k==1
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)-1);
                        else
                            OOALDG=new_resultLDG(2+m+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-3);
                        end
                        dLDG=plot(log10(kappa),OOALDG,'Marker','diamond','MarkerSize',6,'MarkerEdgeColor',couleurLDG,'MarkerFaceColor',[0 0 0]);
                        hold on
                    otherwise 
                        switch k
                            case 1
                                kappa_array=[cDG;cSD;cHU;c_plus];
                                Legend1IP = strcat('IP ; \tau = \tau_{theory}');
                                Legend1LDG = strcat('LDG ; \tau = 0');
                                OOA_arrayIP=new_resultIP(3+2*(i-1),3+(N2-1+(N2-2)+1)*(k-1)+N2-1+(N2-2)-1);
                                OOA_arrayLDG=new_resultLDG(3+N4*(i-1):3+N4-1+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)-1);
                                %OOA_array=new_result(3+N4*(i-1):3+N4-1+N4*(i-1),3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1);
                                kappa_array=[kappa_array;resultLDG(3+4:3+4+N4-4-1,2)];
                                %kappa_array=result(3:3+N2-4-1,2);
                                A=zeros(N4,3);
                                A(:,1)=kappa_array;
                                A(:,3)=OOA_arrayLDG;
                                B=sortrows(A);
                                kappa_array=B(:,1);
                                OOA_arrayLDG=B(:,3);
                                h1IP=plot(log10([10^-8,10^8]),[OOA_arrayIP,OOA_arrayIP],'--','Color',couleurIP,'LineWidth',1.0);
                                hold on
                                h1LDG=plot(log10(kappa_array),OOA_arrayLDG,'+--','Color',couleurLDG,'LineWidth',1.0);
                            case 2
                                kappa_array=[cDG;cSD;cHU;c_plus];
                                Legend2IP = strcat('IP ; \tau = 1.5\tau_{theory}');
                                Legend2LDG = strcat('LDG ; \tau = 0.1');
                                OOA_arrayIP=new_resultIP(3+2*(i-1),3+3+(N2-1+(N2-2)+1)*(k-1)+N2-1+(N2-2)+1-3);
                                OOA_arrayLDG=new_resultLDG(3+N4*(i-1):3+N4-1+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-3);
                                %OOA_array=new_result(3+N4*(i-1):3+N4-1+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-1);
                                kappa_array=[kappa_array;resultLDG(3+4:3+4+N4-4-1,2)];
                                %kappa_array=result(3:3+N2-4-1,2);
                                A=zeros(N4,3);
                                A(:,1)=kappa_array;
                                A(:,3)=OOA_arrayLDG;
                                B=sortrows(A);
                                kappa_array=B(:,1);
                                OOA_arrayLDG=B(:,3);
                                h2IP=plot(log10([10^-8,10^8]),[OOA_arrayIP,OOA_arrayIP],'-.','Color',couleurIP,'LineWidth',1.0);
                                h2LDG=plot(log10(kappa_array),OOA_arrayLDG,':','Color',couleurLDG,'LineWidth',1.0);
                            case 3
                                kappa_array=[cDG;cSD;cHU;c_plus];
                                Legend3 = strcat(' with \tau = 1.5\tau_{theory}');
                                OOA_arrayIP=new_resultIP(3+N4*(i-1):3+N4-1+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-1);
                                %OOA_array=new_result(3+N4*(i-1):3+N4-1+N4*(i-1),3+3+(N2+(N2-1)+1)*(k-1)+N2+(N2-1)+1-1);
                                kappa_array=[kappa_array;resultIP(3+4:3+4+N4-4-1,2)];
                                A=zeros(N4,2);
                                A(:,1)=kappa_array;
                                A(:,2)=OOA_arrayIP;
                                B=sortrows(A);
                                kappa_array=B(:,1);
                                OOA_arrayIP=B(:,2);
                                %kappa_array=result(3:3+N2-4-1,2);
                                h3=plot(log10(kappa_array),OOA_arrayIP,'-','Color',couleur,'LineWidth',1.0);
                        end
                        fin=1;       
                end    
            end
        end
        set(gca,'gridlinestyle',':')
        grid on
        pbaspect([1 1 1])
        ylim([0 Np+1])
        h = findobj(gca,'Type','line');
        legend([h1IP;h2IP;h1LDG;h2LDG],{Legend1IP,Legend2IP,Legend1LDG,Legend2LDG},'Fontsize',8,'Location','southwest','FontName','Arial')
        %title(titre,'Fontsize',20,'Interpreter','LaTex')
        ticks = -8:2:8;
        labels = arrayfun(@(x)sprintf('10^{%i}', x), ticks, 'uni', 0);
        set(gca, 'xtick', ticks, 'XTickLabels', labels)
        xt = get(gca, 'XTick');
        set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',8,'FontName','Arial')
        xlabel('$\kappa$','interpreter','latex','Fontsize',10,'FontName','Arial')
        ylabel('OOA','Fontsize',8,'FontName','Arial')
        xlim([-8 8])
        hold off;
end