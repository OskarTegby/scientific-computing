function [solution, tvec]=rusanov(dx,dt,nframes,T,u1_left,u2_left, ...
                                  u3_left,u1_right,u2_right,u3_right); 

Tstop=T-dt;lambda=dt/dx/2;
n=1;np=2;N=round(4/dx)+1;
u1=zeros(3,N);u2=zeros(3,N);u3=zeros(3,N);
f1(1:N)=0;f2(1:N)=0;f3(1:N)=0;
p(1:N)=0;v(1:N)=0;a(1:N)=0;a2(1:N)=0;
h1(1:N)=0;h2(1:N)=0;h3(1:N)=0;

%-- When to save solution
Nstep=round(Tstop/dt);
clear solution tvec
if (nframes>Nstep)
  savestep=1:Nstep;
else
  savestep=round(Nstep/nframes*(1:nframes));
end
savecount=1;


% Initial conditions
t=0;
for i=1:N
    x=dx*(i-1)-2;
    if x>0
       u1(np,i)=u1_right;
       u2(np,i)=u2_right;
       u3(np,i)=u3_right;
    else
       u1(np,i)=u1_left;
       u2(np,i)=u2_left;
       u3(np,i)=u3_left;
    end;
end

% Tidsstegning.
for step=1:Nstep
    t=t+dt;
    temp=n;n=np;np=temp;

    % Hastigheten och trycket.
    v=u2(n,:)./u1(n,:);
    p=0.4*(u3(n,:)-0.5*u2(n,:).^2./u1(n,:));

    % Flodesfunktionen
    f1=u2(n,:);
    f2=u1(n,:).*v.^2+p;
    f3=(u3(n,:)+p).*v;

    % Inre punkter.
    u1(np,2:N-1)=u1(n,2:N-1)-lambda*(f1(3:N)-f1(1:N-2));
    u2(np,2:N-1)=u2(n,2:N-1)-lambda*(f2(3:N)-f2(1:N-2));
    u3(np,2:N-1)=u3(n,2:N-1)-lambda*(f3(3:N)-f3(1:N-2));

    % Numerisk flodesfunktion
    a=abs(v)+sqrt(1.4*p./u1(n,:));
    a2(1:N-1)=(a(1:N-1)+a(2:N))/2;
    h1(1:N-1)=a2(1:N-1).*(u1(n,2:N)-u1(n,1:N-1));
    h2(1:N-1)=a2(1:N-1).*(u2(n,2:N)-u2(n,1:N-1));
    h3(1:N-1)=a2(1:N-1).*(u3(n,2:N)-u3(n,1:N-1));
    u1(np,2:N-1)=u1(np,2:N-1)+lambda*(h1(2:N-1)-h1(1:N-2));
    u2(np,2:N-1)=u2(np,2:N-1)+lambda*(h2(2:N-1)-h2(1:N-2));
    u3(np,2:N-1)=u3(np,2:N-1)+lambda*(h3(2:N-1)-h3(1:N-2));

    % Randvillkor.
    u1(np,1)=u1(np,3);u1(np,2)=u1(np,3);
    u1(np,N)=u1(np,N-2);u1(np,N-1)=u1(np,N-2);
    u2(np,1)=u2(np,3);u2(np,2)=u2(np,3);
    u2(np,N)=u2(np,N-2);u2(np,N-1)=u2(np,N-2);
    u3(np,1)=u3(np,3);u3(np,2)=u3(np,3);
    u3(np,N)=u3(np,N-2);u3(np,N-1)=u3(np,N-2);


    % Spara losning.
    if (step==savestep(savecount)) 
      solution(savecount,1,:)=u1(n,:);
      solution(savecount,2,:)=u2(n,:);
      solution(savecount,3,:)=u3(n,:);  
      tvec(savecount)=t;
      savecount=savecount+1;
    end
    
end

