% Temperature.m
%
% 05 08 24

clear
disp(['************************************************************' ...
      '************'])
disp('* Fill in the parameters. Pressing Enter gives the default values')
disp(['************************************************************' ...
      '************'])

rho=2300;   % Density
c=800;      % Specific heat 
kappa=1.63; % Heat conduction 
T=3600;     % Solution time

M=input('Give number of gridpoints in space (79 default): ');
if (isempty(M))
  M=79;
end
dt=input('Give the time-step (1 default): ');
if (isempty(dt))
  dt=1;
end
nt=round(T/dt); % Number of steps


dx=0.2/(M+1);                  %space-step
x=[dx:dx:0.2-dx];              %grid points in space


t=0;                           %initiate time
u(1:M,1)=20*ones(M,1);         %initial solution

%%% setup the tridiagonal coefficient-matrix for the right-hand side %%%
c1=dt*kappa/(rho*c*dx^2);  %precompute common factor
B=spdiags(repmat([c1 1-2*c1 c1],M,1),[-1:1], M, M); %form the matrix

for n=1:nt
  u(:)=B*u(:);                  %compute right-hand side
  u(1)=u(1)+c1*20;              %set boundary condition at x=0
  u(M)=u(M)+c1*20-c1*50*t/3600; %set boundary condition at x=0.2 
  t=t+dt;
end

%plot results
clf;
subplot(2,1,1);
plot([0 x 0.2],[20 u' 20-50*t/3600],'bo:')
xlabel('x')
ylabel('Temperature')
title('Temperature at t=3600s')

%load precomputed reference solution
load U.mat;
load X.mat;
hold on;
plot([0 X' 0.2],[20 U' 20-50*t/3600],'r-')
legend('Numerical Solution','Reference Solution')


% Error
subplot(2,1,2)
y=spline(X,U,x)';
plot([0 x 0.2],[0 (y-u)' 0],'go:');axis([0 0.2 min([y-u;0]) max([y-u;0])]);
xlabel('x')
ylabel('Error')
title('Error in the computed solution at t=3600s')
