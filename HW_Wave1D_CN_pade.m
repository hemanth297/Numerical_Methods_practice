clear all; close all

L = 4; N= 64; Tmax = 4;c = 1; dt = 0.025;
NR_Wave1D_CN_Pade_interleavingqv(L,Tmax,c,N,dt);

function NR_Wave1D_CN_Pade_interleavingqv(L,Tmax,c,N,dt)

dx=L/N; IterSteps=2; t=0; x=(-N/2:N/2-1)'*dx; q=exp(-x.^2/0.1); v=0*q;
y(1:2:2*N) = q; y(2:2:2*N) = v; y = y'; %% interleaqving q,v
E = eye(2*N,2*N); F = zeros(2*N,2*N);
for i = 1:N
 F(2*i-1, 2*i) = 1;
 F(2*i,2*i-1) = -2.4/dx^2;
end
for i = 1:N-1
    E(2*i, 2*(i+1)) = 0.1;
    E(2*(i+1), 2*i) = 0.1;
    F(2*i, 2*i+1) = 1.2/dx^2;
    F(2*(i+1), 2*i-1) = 1.2/dx^2;
end
%%%% making matrix circulent
%E(2,2*N) = E(2,4); E(2*N,2) = E(2,4);  % use only for periodic BC
%F(2*N,1) = F(2,3); F(2,2*N-1) = F(2,3); % use only for periodic BC

F = c*F;
A = E/dt - 1/2*F ; B = E/dt +1/2*F; 
for n = 1:Tmax/dt
    y = A\(B*y);
    t=t+dt; q=y(1:2:2*N); v=(2:2:2*N); PlotXY_comp(x,q,t,-L/2,L/2,-1.2,1.2);
end


end % function NR_Wave1D_ItCN_Pade


function PlotXY_comp(x,y,t,xmin,xmax,ymin,ymax)
plot(x,y,'bo',x,exp(-x.^2/0.1),'r--');  xlabel('x');  ylabel('q'); legend("live plot","initial plot")
title(sprintf('Time = %5.2f',t)); axis([xmin xmax ymin ymax]); pause(0.0001);
end % function PlotXY