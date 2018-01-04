%---------------------------------------------------------------
% 
%  Leap-frog scheme to solve the wave equation
%     u_t(x,t) +c u_x(x,t) = 0,        xl < x < xr, 0 < t < tf
%     u(x, 0) = f(x),                xl < x < xr
%  periodic boundary
% 
% The analytic solution is:
%        u(x,t)= f(x-t)=e^(-50(x-t-0.5)^2)
%---------------------------------------------------------------
clear all;                  
close all;

xl=0; xr=1;                 % x domain [xl,xr]
M = 150;                    % M: number of division for x
dx = (xr-xl) / M;           % dx: mesh size
tf = 1.0;                   % final simulation time
Nt = 1000;                  % Nt: number of time steps
dt = tf/Nt;                 
c1 = 50;                     % parameter for the solution
c=1;
mu =c*dt/dx;

% Evaluate the initial conditions
x = xl : dx : xr;              % generate the grid point
f = exp(-c1*(x-0.3).^2);        % dimension f(1:M+1) 
f0= exp(-c1*(x-0.3-dt).^2);
figure(1)
plot(x,f,'b-');
ylim([-0.5,1.5])
xlabel('x')
ylabel('u')
%pause()
% store the solution at all grid points for all time steps
u = zeros(Nt,M+1);   
u(1,:)=f;
udata=f;
tdata=0;
% Find the approximate solution at each time step
for n = 2:Nt
    t = (n-1)*dt;                    % current time
     if n==2               % first time step
        for m=2:M          % interior nodes     
 	       u(n,m) = u(n-1,m) - mu*(u(n-1,m+1)-u(n-1,m-1));
        end     
        u(n,1) = u(n,M);       % the left-end point
        u(n,M+1) = u(n,2);     % the right-end point 
     else 
       for m=2:M          % interior nodes
       u(n,m) = u(n-2,m) - mu*(u(n-1,m+1)-u(n-1,m-1));
       end
       u(n,1) = u(n,M);       % the left-end point
       u(n,M+1) = u(n,2);     % the right-end point 
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
view(10, 60)
axis([xl xr 0 tmax -1 3]); 
xlabel('x')
ylabel('t')
zlabel('u')
grid off 

