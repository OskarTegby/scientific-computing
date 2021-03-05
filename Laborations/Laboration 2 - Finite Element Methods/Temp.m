% Temp.m
clear
disp(['************************************************************' ...
      '************'])
disp('* Fill in the parameters. Pressing Enter gives the default values')
disp(['************************************************************' ...
      '************'])

hprim=input('Give the heat transfer coefficient (1.63 default): ');
if (isempty(hprim))
  hprim=1.63;
end
M=input('Give number of gridpoints in space (5 default): ');
if (isempty(M))
  M=5;
end
T_a=input('Give air temperature (20 default)');
if (isempty(T_a))
  T_a=20;
end

tic;
dx=10/(M+1);
x=[dx:dx:10-dx];

rhs=-hprim*dx^2*T_a*ones(M,1); % Why is T zero?
rhs(1)=rhs(1)-40;
rhs(M)=rhs(M)-200;

A=zeros(M,M);
A(1,1)=-2-hprim*dx^2; % Why do we get this?
A(1,2)=1;
for j=2:M-1
  A(j,j-1)=1;
  A(j,j)=-2-hprim*dx^2;
  A(j,j+1)=1;
end
A(M,M-1)=1;
A(M,M)=-2-hprim*dx^2;

u=A\rhs; % Finding the T_i by solving for the right-hand side using Gaussian elimination.
t=toc;
plot([0 x 10],[40 u' 200],'.-')
xlabel('Length')
ylabel('Temp')
disp(' ');
disp(['Computational Time t = ' num2str(t) ' sec']);
