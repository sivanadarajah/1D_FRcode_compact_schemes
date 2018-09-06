function [ M ] = RK54_VN( S,dt,Np )

M=diag(ones(1,Np))+dt*S+(dt*S)^2/2+(dt*S)^3/6+(dt*S)^4/24+(dt*S)^5/200;

end

