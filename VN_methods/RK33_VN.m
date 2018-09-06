function [ M ] = RK33_VN( S,dt,Np )

M=diag(ones(1,Np))+dt*S+(dt*S)^2/2+(dt*S)^3/6;

%M=diag(ones(1,Np))+dt*S+(dt*S)^2/2+(dt*S)^3/6;

%M=diag(ones(1,Np))
end

