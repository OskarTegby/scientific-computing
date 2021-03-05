clear all
clf

%
% Space discretization
%
M=24;
dx=2/(M+1);

x=[dx:dx:2-dx];

%
% Time discretization
%
dt=0.0005; 
N=200/dt;

%
% Problem parameters
%
D_D=0.28;
D_E=0.6;
sigma_1=20;
sigma_2=0.0063;
sigma_3=0.04;
sigma_4=0.8;
sigma_1p=0.028;
sigma_4p=0.027;

%
% Load initial state
%
load 'rho_DD.mat'
load 'rho_d.mat'
load 'rho_EE.mat'
load 'rho_e.mat'

i=0;

tic;

for n=1:N
if (mod(n,1000)==0)
  disp(['Timestep nr.',num2str(n),' of ',num2str(N)])
end
  %
  % Approximate second order derivatives with D+D-
  %
  diff1(2:M-1,1)=D_D*(rho_D(1:M-2)-2*rho_D(2:M-1)+rho_D(3:M))/dx^2;
  diff1(1,1)=D_D*(rho_D(2)-rho_D(1))/dx^2;
  diff1(M,1)=D_D*(-rho_D(M)+rho_D(M-1))/dx^2;

  diff3(2:M-1,1)=D_E*(rho_E(1:M-2)-2*rho_E(2:M-1)+rho_E(3:M))/dx^2;
  diff3(1,1)=D_E*(rho_E(2)-rho_E(1))/dx^2;
  diff3(M,1)=D_E*(-rho_E(M)+rho_E(M-1))/dx^2;
  %
  % Evaluate right hand side terms
  %
  term1=sigma_1.*rho_D./(1.0+sigma_1p.*rho_e);
  term2=sigma_2.*rho_e.*rho_d;    
  term3=sigma_3.*rho_E.*rho_D;    
  term4=sigma_4.*rho_e./(1.0+sigma_4p.*rho_D);
  %
  % Time stepping with the Euler forward method
  %
  rho_D=rho_D+dt*(diff1-term1+term2);
  rho_d=rho_d+dt*(term1-term2);
  rho_E=rho_E+dt*(diff3-term3+term4);
  rho_e=rho_e+dt*(term3-term4);
  
  if (mod(n,10000)==0)
    i=i+1;
    rhosave_D(:,i)=rho_D+rho_d;
    rhosave_E(:,i)=rho_E+rho_e;
    t(i)=n*dt;
  end
end

comptime=toc

for k=1:4
  figure(k),clf
  set(gca,'FontSize',16)
end
  
[T X]=meshgrid(t,x);
figure(1)
surf(X,T,rhosave_D)
drawnow
view(2)
shading interp
H=colorbar; set(H,'FontSize',16)
xlabel('x')
ylabel('t')
title('\rho_D+\rho_d')

figure(2)
surf(X,T,rhosave_E)
drawnow
view(2)
shading interp
H=colorbar; set(H,'FontSize',16)
xlabel('x')
ylabel('t')
title('\rho_E+\rho_e')

figure(3)
plot(x,mean(rhosave_D'))
drawnow
xlabel('x')
title('Average of \rho_D+\rho_d')

figure(4)
plot(x,mean(rhosave_E'))
drawnow
xlabel('x')
title('Average of \rho_E+\rho_e')




