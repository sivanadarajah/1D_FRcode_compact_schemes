function M=creationM(ksi,Np,hk)

M=zeros(Np,Np);
for i=1:1:Np
    Y_i=zeros(1,Np);
    Y_i(i)=1;
    Lagrange_i=lagrangepoly(ksi,Y_i);
    for j=1:1:Np
        Y_j=zeros(1,Np);
        Y_j(j)=1;
        Lagrange_j=lagrangepoly(ksi,Y_j);
        P=conv(Lagrange_i,Lagrange_j);
        P_int=polyint(P,0);
        M(i,j)=diff(polyval(P_int,[-1 1]));
    end
end
M=M*hk/2;
end

