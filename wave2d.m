%
% 2D Wave Equation
% Euler Explicit Scheme
% Boundary u=0
%
close all
clear all

dx=0.01; dy=0.01; dt=0.001;
x=0:dx:1; y=0:dy:1; lx=length(x); ly=length(y);

%f= @(x,y) max(0,0.1-2*(x+y-5).^2);
%u=f
%ub=f;

u=zeros(lx,ly);
ub=zeros(lx,ly);
un=zeros(lx,ly);
u(50:55,50:55)=4;
ub(50:55,50:55)=4;
u=u';
ub=ub';

figure(1)
imagesc(x,y,u)
set(gca,'ydir','normal')
    xlabel('x')
    ylabel('y')
    colorbar

figure(2)
 surf(x,y,u);     % 3-D surface plot
    zlim([0,5])
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title('Numerical solution of 2-D wave equation')
%    rotate3d on
%    pause()
    close all
for tt=1:500; %timestep
    for ii=2:ly-1;
        for jj=2:lx-1;
     un(jj,ii)=2*u(jj,ii)-ub(jj,ii)+(dt^2/dx^2)*(u(jj+1,ii)+u(jj-1,ii)...
        -2*u(jj,ii))+(dt^2/dy^2)*(u(jj,ii+1)+u(jj,ii-1)-2*u(jj,ii));
        end
    end
%    figure(3)
%    imagesc(x,y,un)
%    set(gca,'ydir','normal')
%    xlabel('x')
%    ylabel('y')
%    colorbar
%    hold off
    figure(4)
    surf(x,y, un);     % 3-D surface plot
    view([-37,80]);
    zlim([-1,5])
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title('Numerical solution of 2-D wave equation')
%    rotate3d on
%    hold off
   drawnow
    ub=u;
    u=un;
end