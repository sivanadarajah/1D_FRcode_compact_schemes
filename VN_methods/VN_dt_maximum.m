function [ lambda_max ] = VN_dt_maximum( Np,A,dt,choice_scheme )
%{
This function has for purpose to calculate the maximum time step in order
to be stable for a given value of tau
We need to precise a time solver (it is explictly RK4 5stage 4th order)
With this solver we can calculate the matrix M and its spectra if the
maximum of its eigenvalue is greater than one then the time step is too
large
%}
w_array=0+pi/500:pi/500:pi-pi/500;
N3=max(size(w_array));
lambda_array=zeros(Np,N3);
%-------------------------------------------------------------------------
%Compact
if strcmp(choice_scheme,'compact')
    C_m_1_final=A(:,1:Np);
    C0_final=A(:,Np+1:2*Np);
    C_p_1_final=A(:,2*Np+1:3*Np);
    for i=1:1:N3
        w=w_array(i);
        S=exp(-1i*w)*C_m_1_final+C0_final+...
            exp(1i*w)*C_p_1_final;
        M=RK54_VN(S,dt,Np);
        [~,Dia]=eig(M);
        %These are the different eigenvalues for the different corrections
        %functions
        for j=1:1:Np
            lambda(j)=Dia(j,j);
        end
        %We take a set of values for w and calculate the eigenvalues
        % We go all over the spectrum [0,2\pi]
        for j=1:1:Np
            lambda_array(j,i)=abs(double(lambda(j)));
        end
    end
[lambda_max]=max(max(abs(lambda_array)))
end
%-------------------------------------------------------------------------
%Wide stencil
if strcmp(choice_scheme,'wide')
    C_m_2_final=A(:,1:Np);
    C_m_1_final=A(:,Np+1:2*Np);
    C0_final=A(:,2*Np+1:3*Np);
    C_p_1_final=A(:,3*Np+1:4*Np);
    C_p_2_final=A(:,4*Np+1:5*Np);
    for i=1:1:N3
        w=w_array(i);
        S=exp(-1i*2*w)*C_m_2_final+exp(-1i*w)*C_m_1_final+C0_final+...
    exp(1i*w)*C_p_1_final+exp(1i*2*w)*C_p_2_final;
        M=RK54_VN(S,dt,Np);
        [~,Dia]=eig(M);
        %These are the different eigenvalues for the different corrections
        %functions
        for j=1:1:Np
            lambda(j)=Dia(j,j);
        end
        %We take a set of values for w and calculate the eigenvalues
        % We go all over the spectrum [0,2\pi]
        for j=1:1:Np
            lambda_array(j,i)=abs(double(lambda(j)));
            %lambda_array(j,i)=max(real(double(lambda(j))));
        end
    end
[lambda_max]=max(max(abs(lambda_array)));
end

end

