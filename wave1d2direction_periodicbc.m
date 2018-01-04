%---------------------------------------------------------------
% 
%  difference scheme to solve the wave equation
%     u_{tt}(x,t) -c^2 u_{xx}(x,t) = 0,        xl < x < xr, 0 < t < tf
%     u(x, 0) = f(x),                 xl < x < xr
%     periodic boundary condition
% 
% The analytic solution is:
%        u(x,t)= f(x-t)=e^(-50(x-t-0.3)^2)
%---------------------------------------------------------------
clear all;                  
close all;

xl=0; xr=2;                 % x domain [xl,xr]
M = 200;                    % M: number of division for x
dx = (xr-xl) / M;           % dx: mesh size
tf = 1.5;                   % final simulation time
Nt = 2000;                  % Nt: number of time steps
dt = tf/Nt;                 
c1 = 50;                     % parameter for the solution
c  = 1;
mu = c.^2.*dt/dx; 

% Evaluate the initial conditions
x = xl : dx : xr;              % generate the grid point
f = exp(-c1*(x-1).^2);        % dimension f(1:M+1) 
g = zeros(1,M+1);
figure(2)
plot(x,f,'b-');
ylim([-0.5,1.5])
xlabel('x')
ylabel('u')
%pause()
% store the solution at all grid points for all time steps
u = zeros(Nt,M+1);   
udata=f;
tdata=0;

% Find the approximate solution at each time step
for n = 1:Nt
    t = n*dt;                    % current time
    if n==1               % first time step
        for m=2:M
       u(n,m)=f(m);
        end
	u(n,1) = u(n,M);       % the left-end point
	u(n,M+1) = u(n,2);     % the right-end poin
    elseif n==2
        for m=2:M
		u(n,m)=dt.*g(m)+u(n-1,m);
        end
	u(n,1) = u(n,M);       % the left-end point
	u(n,M+1) = u(n,2);     % the right-end poin        
    else
       for m=2:M          % interior nodes     
	       u(n,m)=mu.^2*(u(n-1,m+1)+u(n-1,m-1))+2.*(1-mu.^2).*u(n-1,m)-u(n-2,m);
       end
       	u(n,1) = u(n,M);       % the left-end point
	u(n,M+1) = u(n,2);     % the right-end poin
    end
        figure(1)
	plot(x,u(n,:),'b-');
ylim([-0.5,1.5])
xlabel('x')
ylabel('u')
 t=n*dt;
    title(sprintf('t=%5.2f',t))
hold off
%pause(eps)
drawnow
    if mod(n,40)==0
      udata = [udata; u(n,:)];
tdata = [tdata n*dt];
end
end
tmax=Nt*dt;
figure(2)
waterfall(x,tdata,udata)
%colormap(jet(128)); 
view(5, 40)
axis([xl xr 0 tmax -0.5 3]); 
xlabel('x')
ylabel('t')
zlabel('u')
grid off 


