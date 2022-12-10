clear all;close all

HW_Wave1D_CN_Pade(4,4,1,64,0.01);

function HW_Wave1D_CN_Pade(L,Tmax,c,N,dt)
% function NR_Wave1D_CN_Pade(L,Tmax,c,N,dt)
% This script simulates the 1D Wave equation with periodic boundary conditions.
% Iterative CN in time is used with a fourth-order Pade method in space.

dx=L/N; IterSteps=2; t=0; x=(-N/2:N/2-1)'*dx; q=exp(-x.^2/0.1); v=0*q;
y = [q;v];
E = zeros(2*N,2*N);
F = E; 
C = eye(N); D = -2.4*eye(N);
for i = 1:N-1
    C(i,i+1) = 0.1;  D(i,i+1) = 1.2; 
    C(i+1,i) = 0.1;  D(i,i+1) = 1.2; 
end
D = D/dx^2;
E(1:N,1:N) = eye(N); E(N+1:2*N,N+1:2*N) = C;
F(1:N,N+1:2*N) = eye(N); F(N+1:2*N,1:N) = D;
A = E/dt - 1/2*F ; B = E/dt +1/2*F; 

for n = 1:Tmax/dt
    y = A\(B*y);
    t=t+dt; q=y(1:N); v=(N+1:2*N); NR_PlotXY(x,q,t,-L/2,L/2,-0.2,1.2);
end

end 
% function NR_Wave1D_ItCN_Pade
