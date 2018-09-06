function [ coeff_u_ini,coeff_u_end,Aini,Aini2,Aint,Aend2,Aend,B_ini,B_end ] = local_matrices3( Np,ksi,hk,invM,beta,tau,s,c_correction,kappa,choice_scheme,choice_analysis )
global c
global alpha
global mu
global xini
global xfin
%{
This is the most important function of the codes :
 it creates the local matrices for any problem : diffusion,
 diffusion/advection,...
Every test uses this function
%}
%% Corections functions && Derivative operator
[g_left,g_right,h_left,h_right]=corrections_functions(ksi,Np,c_correction,kappa);%for both the axuiliary variable and main equation
[D_operator,Lagrange]=derivative_Lagrange(ksi,Np);
l=zeros(Np,1);
r=zeros(Np,1);
for k=1:1:Np
    l(k)=polyval(Lagrange(k,:),-1);
    r(k)=polyval(Lagrange(k,:),1);
end
%% Creation local matrices
C_auxiliare=sparse(Np,Np);
C_main=sparse(Np,Np);
Jump_ini=sparse(Np,2*Np);%For the first cell ; only the second cell is involved
Jump=sparse(Np,3*Np);%Both neighbour to cell K are involved
Jump_end=sparse(Np,2*Np);%Only the second_to_last participate to the last cell
Me_ini=sparse(Np,2*Np);%Average operator
Me_loc=sparse(Np,3*Np);
Me_end=sparse(Np,2*Np);
I_boundary=sparse(Np,Np);

%-----------------------------------------
%-----------------------------------------
Top=zeros(1,Np);
Bot=zeros(1,Np);
Top(1)=1;
Bot(Np)=1;

C_auxiliaire=g_left'*Top+g_right'*Bot;
C_main=h_left'*Top+h_right'*Bot;


%-----------------------------------------
%-----------------------------------------
Jump_ini(1,1:Np)=zeros(1,Np);%l';%  %Initial conditions --> No jump at the beginning
Jump_ini(Np,1:Np)=r';
Jump_ini(Np,Np+1:2*Np)=-l';%Participation of the second_element on the first cell
%------------------------------------------
Jump(1,1:Np)=r';%Participation of the left neighbor
Jump(1,Np+1:2*Np)=-l';
Jump(Np,Np+1:2*Np)=r';
Jump(Np,2*Np+1:3*Np)=-l';%Participation of the right neighbor
%-------------------------------------------
Jump_end(1,1:Np)=r';
Jump_end(1,Np+1:2*Np)=-l';
Jump_end(Np,Np+1:2*Np)=zeros(1,Np);%-r';%   %Initial conditions --> No jump at the end
%-----------------------------------------
%-----------------------------------------
Me_ini(1,1:Np)=2*l';
Me_ini(Np,1:Np)=r';
Me_ini(Np,Np+1:2*Np)=l';%Participation of the second_element on the first cell
%------------------------------------------
Me_loc(1,1:Np)=r';%Participation of the left neighbor
Me_loc(1,Np+1:2*Np)=l';
Me_loc(Np,Np+1:2*Np)=r';
Me_loc(Np,2*Np+1:3*Np)=l';%Participation of the right neighbor
%-------------------------------------------
Me_end(1,1:Np)=r';
Me_end(1,Np+1:2*Np)=l';
Me_end(Np,Np+1:2*Np)=2*r';%Initial conditions
%-----------------------------------------
%-----------------------------------------
I_boundary(1,1:Np)=l';
I_boundary(Np,1:Np)=r';
if s~=0
    [Re_ini,Re_int,Re_end]=BR2_coeffs(Np,l,r,invM,Jump_ini,Jump,Jump_end);
else
    Re_ini=zeros(Np,2*Np);
    Re_int=zeros(Np,3*Np);
    Re_end=zeros(Np,2*Np);
end
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
if strcmp(choice_analysis,'TC')
%-------------------------
    %We compute the auxiliary equation
    Qini=zeros(Np,2*Np);
    Qint=zeros(Np,3*Np);
    Qend=zeros(Np,2*Np);
    %Qini
    Qini(1:Np,1:Np)=D_operator-C_auxiliaire*I_boundary;%we get Q=d/dx -u *g_left-u *g_right
    Qini=Qini+C_auxiliaire*(1/2*Me_ini-beta*Jump_ini);%we get Q=d/dx (u*-u) *g_left+(u*-u) *g_right
    %Qint
    Qint(1:Np,Np+1:2*Np)=D_operator-C_auxiliaire*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
    Qint=Qint+C_auxiliaire*(1/2*Me_loc-beta*Jump);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
    %Qend
    Qend(1:Np,Np+1:2*Np)=D_operator-C_auxiliaire*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
    Qend=Qend+C_auxiliaire*(1/2*Me_end-beta*Jump_end);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right

    %-------------------------------------------------------------------------------------------
    if strcmp(choice_scheme,'compact')
        %We compute the main equation for a compact scheme
        Deltaini=sparse(Np,2*Np);
        Deltaint=sparse(Np,3*Np);
        Deltaend=sparse(Np,2*Np);
        %Aini
        Deltaini=(D_operator-C_main*I_boundary)*Qini;%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaini=Deltaini+C_main*(1/2*Me_ini*blkdiag(D_operator,D_operator)-tau*hk/2*Jump_ini+s*hk/2*Re_ini);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{d/dx}}-tau[| |]
        Deltaini=mu*4/hk^2*Deltaini;
        %Aint
        Deltaint=(D_operator-C_main*I_boundary)*Qint;%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaint=Deltaint+C_main*(1/2*Me_loc*blkdiag(D_operator,D_operator,D_operator)-tau*hk/2*Jump+s*hk/2*Re_int);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{d/dx}}-tau[| |]
        Deltaint=mu*4/hk^2*Deltaint;
        Deltaini2=Deltaint;
        Deltaend2=Deltaint;
        %Aend
        Deltaend=(D_operator-C_main*I_boundary)*Qend;%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaend=Deltaend+C_main*(1/2*Me_end*blkdiag(D_operator,D_operator)-tau*hk/2*Jump_end+s*hk/2*Re_end);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{d/dx}}-tau[| |]    
        Deltaend=mu*4/hk^2*Deltaend;
    end
    if strcmp(choice_scheme,'wide')
        %We compute the main equation for a wide stencil (LDG)
        Deltaini=sparse(Np,3*Np);
        Deltaini2=sparse(Np,4*Np);
        Deltaint=sparse(Np,5*Np);
        Deltaend2=sparse(Np,4*Np);
        Deltaend=sparse(Np,3*Np);
        %Aini
        Deltaini(1:Np,1:2*Np)=(D_operator-C_main*I_boundary)*Qini+C_main*(-hk/2*tau*Jump_ini);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaini=Deltaini+C_main*((1/2*Me_ini+beta*Jump_ini)*blk2(Qini,Qint,Np));%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{Q}}+beta[|Q|]-tau[| |]
        Deltaini=mu*4/hk^2*Deltaini;
        %Aini2
        Deltaini2(1:Np,1:3*Np)=(D_operator-C_main*I_boundary)*Qint+C_main*(-hk/2*tau*Jump);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaini2=Deltaini2+C_main*((1/2*Me_loc+beta*Jump)*blk3(Qini,Qint,Qint,Np));%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{Q}}+beta[|Q|]-tau[| |]
        Deltaini2=mu*4/hk^2*Deltaini2;
        %Aint
        Deltaint(1:Np,Np+1:4*Np)=(D_operator-C_main*I_boundary)*Qint+C_main*(-hk/2*tau*Jump);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaint=Deltaint+C_main*(1/2*Me_loc+beta*Jump)*blk3(Qint,Qint,Qint,Np);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{Q}}+beta[|Q|]-tau[| |]
        Deltaint=mu*4/hk^2*Deltaint;
        %Aend2
        Deltaend2(1:Np,Np+1:4*Np)=(D_operator-C_main*I_boundary)*Qint+C_main*(-hk/2*tau*Jump);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaend2=Deltaend2+C_main*(1/2*Me_loc+beta*Jump)*blk3(Qint,Qint,Qend,Np);%we get Q=d/dx (u*-u)*g_left(u*-u)*g_right
        Deltaend2=mu*4/hk^2*Deltaend2;
        %Aend
        Deltaend(1:Np,Np+1:3*Np)=(D_operator-C_main*I_boundary)*Qend+C_main*(-hk/2*tau*Jump_end);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaend=Deltaend+C_main*(1/2*Me_end+beta*Jump_end)*blk2(Qint,Qend,Np);%we get Q=d/dx (u*-u)*g_left(u*-u)*g_right
        Deltaend=mu*4/hk^2*Deltaend;
    end
    Aini=-Deltaini;
    Aini2=-Deltaini2;
    Aint=-Deltaint;
    Aend2=-Deltaend2;
    Aend=-Deltaend;
    if c~=0
        %We compute an advective term
        Convini=zeros(Np,2*Np);
        Convint=zeros(Np,3*Np);
        Convend=zeros(Np,2*Np);
        %Qini
        Convini(1:Np,1:Np)=D_operator-C_main*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
        Convini=Convini+C_main*(1/2*Me_ini+(1-alpha)/2*Jump_ini);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
        Convini=c*2/hk*Convini;
        %Qint
        Convint(1:Np,Np+1:2*Np)=D_operator-C_main*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
        Convint=Convint+C_main*(1/2*Me_loc+(1-alpha)/2*Jump);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
        Convint=c*2/hk*Convint;
        %Qend
        Convend(1:Np,Np+1:2*Np)=D_operator-C_main*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
        Convend=Convend+C_main*(1/2*Me_end+(1-alpha)/2*Jump_end);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
        Convend=c*2/hk*Convend;
        lenghtint=max(size(Aint));
        if floor(lenghtint/Np)==5 %wide stencil
            k=2;
        elseif floor(lenghtint/Np)==3 %compact stencil
            k=1;
        else
            error('Problem with the matrix A');
        end
        Aini(1:Np,1:2*Np)=Aini(1:Np,1:2*Np)+Convini;
        Aini2(1:Np,1:3*Np)=Aini2(1:Np,1:3*Np)+Convint;
        Aint(1:Np,1+(k-1)*Np:(k+2)*Np)=Aint(1:Np,1+(k-1)*Np:(k+2)*Np)+Convint;
        Aend2(1:Np,1+(k-1)*Np:(k+2)*Np)=Aend2(1:Np,1+(k-1)*Np:(k+2)*Np)+Convint;
        Aend(1:Np,1+(k-1)*Np:(k+1)*Np)=Aend(1:Np,1+(k-1)*Np:(k+1)*Np)+Convend;
    end
%----------------------------------------------------------------------------------
%Boundary conditions
    B_ini=sparse(2*Np,1);
    B_end=sparse(2*Np,1);
    coeff_u_ini=zeros(1,Np);
    coeff_u_end=zeros(1,Np);
    for i=1:1:Np
        if i==1
            coeff_u_ini(i)=1/l(1);
            coeff_u_end(i)=-r(i)/r(Np);
        elseif i==Np
            coeff_u_ini(i)=-l(i)/l(1);
            coeff_u_end(i)=1/r(Np);
        else
            coeff_u_ini(i)=-l(i)/l(1);
            coeff_u_end(i)=-r(i)/r(Np);
        end
    end
    if strcmp(choice_scheme,'compact')
        B_ini(1:Np,1)=B_ini(1:Np,1)-Aini(:,1)*coeff_u_ini(1);
        B_ini(Np+1:2*Np,1)=B_ini(Np+1:2*Np,1)-Aini2(:,1)*coeff_u_ini(1);
        B_end(1:Np,1)=B_end(1:Np,1)-Aend2(:,3*Np)*coeff_u_end(Np);
        B_end(Np+1:2*Np,1)=B_end(Np+1:2*Np,1)-Aend(:,2*Np)*coeff_u_end(Np);        
        for i=1:1:Np
            if i==1
                Aend2(:,2*Np+i)=Aend2(:,2*Np+i)+Aend2(:,3*Np)*coeff_u_end(i);
                Aend(:,Np+i)=Aend(:,Np+i)+Aend(:,2*Np)*coeff_u_end(i);
            elseif i==Np
                Aini(:,i)=Aini(:,i)+Aini(:,1)*coeff_u_ini(i);
                Aini2(:,i)=Aini2(:,i)+Aini2(:,1)*coeff_u_ini(i);
            else
                Aend2(:,2*Np+i)=Aend2(:,2*Np+i)+Aend2(:,3*Np)*coeff_u_end(i);
                Aend(:,Np+i)=Aend(:,Np+i)+Aend(:,2*Np)*coeff_u_end(i);
                Aini(:,i)=Aini(:,i)+Aini(:,1)*coeff_u_ini(i);
                Aini2(:,i)=Aini2(:,i)+Aini2(:,1)*coeff_u_ini(i);
            end
        end
    end
    if strcmp(choice_scheme,'wide')
        B_ini(1:Np,1)=B_ini(1:Np,1)-Aini(:,1)*coeff_u_ini(1);
        B_ini(Np+1:2*Np,1)=B_ini(Np+1:2*Np,1)-Aini2(:,1)*coeff_u_ini(1);
        B_end(1:Np,1)=B_end(1:Np,1)-Aend2(:,4*Np)*coeff_u_end(Np);
        B_end(Np+1:2*Np,1)=B_end(Np+1:2*Np,1)-Aend(:,3*Np)*coeff_u_end(Np);        
        for i=1:1:Np
            if i==1
                Aend2(:,3*Np+i)=Aend2(:,3*Np+i)+Aend2(:,4*Np)*coeff_u_end(i);
                Aend(:,2*Np+i)=Aend(:,2*Np+i)+Aend(:,3*Np)*coeff_u_end(i);
            elseif i==Np
                Aini(:,i)=Aini(:,i)+Aini(:,1)*coeff_u_ini(i);
                Aini2(:,i)=Aini2(:,i)+Aini2(:,1)*coeff_u_ini(i);
            else
                Aend2(:,3*Np+i)=Aend2(:,3*Np+i)+Aend2(:,4*Np)*coeff_u_end(i);
                Aend(:,2*Np+i)=Aend(:,2*Np+i)+Aend(:,3*Np)*coeff_u_end(i);
                Aini(:,i)=Aini(:,i)+Aini(:,1)*coeff_u_ini(i);
                Aini2(:,i)=Aini2(:,i)+Aini2(:,1)*coeff_u_ini(i);
            end
        end
    end
end
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
%-------------------------------------------------------------------------------------------
if strcmp(choice_analysis,'VN')
    Qint=zeros(Np,3*Np);
    Qint(1:Np,Np+1:2*Np)=D_operator-C_auxiliaire*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
    Qint=Qint+C_auxiliaire*(1/2*Me_loc-beta*Jump);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
    if strcmp(choice_scheme,'compact')
        Deltaint=zeros(Np,3*Np);
        Deltaint=(D_operator-C_main*I_boundary)*Qint;%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaint=Deltaint+C_main*(1/2*Me_loc*blkdiag(D_operator,D_operator,D_operator)-tau*hk/2*Jump);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{d/dx}}-tau[| |]
        Deltaint=mu*4/hk^2*Deltaint;
    end
    if strcmp(choice_scheme,'wide')
        Deltaint=zeros(Np,5*Np);
        Deltaint(1:Np,Np+1:4*Np)=(D_operator-C_main*I_boundary)*Qint+C_main*(-hk/2*tau*Jump);%we get A=dQ/dx+ -q *g_left-q *g_right
        Deltaint=Deltaint+C_main*(1/2*Me_loc+beta*Jump)*blk3(Qint,Qint,Qint,Np);%we get A=dQ/dx+(q*-q)*g_left+(q*-q)*g_rightt // q*={{Q}}+beta[|Q|]-tau[| |]
        Deltaint=mu*4/hk^2*Deltaint;
    end
    Aint=Deltaint;
    Aini=0;
    Aini2=0;
    Aend2=0;
    Aend=0;
    coeff_u_ini=0;
    coeff_u_end=0;
    B_ini=0;
    B_end=0;
    if c~=0
        Convint=zeros(Np,3*Np);
        Convint(1:Np,Np+1:2*Np)=D_operator-C_main*I_boundary;%we get Q=d/dx -u*g_left-u*g_right
        Convint=Convint+C_main*(1/2*Me_loc+(1-alpha)/2*Jump);%we get Q=d/dx (u*-u)*g_left+(u*-u)*g_right
        Convint=c*2/hk*Convint;
        lenghtint=max(size(Aint));
        if floor(lenghtint/Np)==5 %wide stencil
            k=2;
        elseif floor(lenghtint/Np)==3 %compact stencil
            k=1;
        else
            error('Problem with the matrix A');
        end
        Aint(1:Np,1+(k-1)*Np:(k+2)*Np)=Aint(1:Np,1+(k-1)*Np:(k+2)*Np)-Convint;
    end
end

end

