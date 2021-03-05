function [solution, tvec]=leapfrog(dx,dt,nframes,T,u1_left,u2_left, ...
                                   u3_left,u1_right,u2_right,u3_right); 

Tstop=T-dt;lambda=dt/dx;
nm=1;n=2;np=3;N=round(4/dx)+1;
u1=zeros(3,N);u2=zeros(3,N);u3=zeros(3,N);
f1(1:N)=0;f2(1:N)=0;f3(1:N)=0;
p(1:N)=0;v(1:N)=0;

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
for i=1:N
    x=dx*(i-1)-2;
    if x>0
       u1(n,i)=u1_right;
       u2(n,i)=u2_right;
       u3(n,i)=u3_right;
    else
       u1(n,i)=u1_left;
       u2(n,i)=u2_left;
       u3(n,i)=u3_left;
    end;
end

  
% Time step 1
u1(np,:)=u1(n,:);u2(np,:)=u2(n,:);u3(np,:)=u3(n,:);
if (savestep(savecount)==1) 
  solution(savecount,1,:)=u1(np,:);
  solution(savecount,2,:)=u2(np,:);
  solution(savecount,3,:)=u3(np,:);  
  tvec(savecount)=dt;
  savecount=savecount+1;
end

% Time stepping
t=0;
for step=1:Nstep
    t=t+dt;
    temp=nm;nm=n;n=np;np=temp;

    % Hastigheten och trycket.
    v=u2(n,:)./u1(n,:);
    p=0.4*(u3(n,:)-0.5*u2(n,:).^2./u1(n,:));

    % Flodesfunktionen
    f1=u2(n,:);
    f2=u1(n,:).*v.^2+p;
    f3=(u3(n,:)+p).*v;

    % Inre punkter.
    u1(np,2:N-1)=u1(nm,2:N-1)-lambda*(f1(3:N)-f1(1:N-2));
    u2(np,2:N-1)=u2(nm,2:N-1)-lambda*(f2(3:N)-f2(1:N-2));
    u3(np,2:N-1)=u3(nm,2:N-1)-lambda*(f3(3:N)-f3(1:N-2));
    
        % Artificiell viskositet.
        k=0.05;
    u1(np,3:N-2)=u1(np,3:N-2)-k*(u1(nm,5:N)-4*u1(nm,4:N-1)+...
        6*u1(nm,3:N-2)-4*u1(nm,2:N-3)+u1(nm,1:N-4));
    u2(np,3:N-2)=u2(np,3:N-2)-k*(u2(nm,5:N)-4*u2(nm,4:N-1)+...
        6*u2(nm,3:N-2)-4*u2(nm,2:N-3)+u2(nm,1:N-4));
    u3(np,3:N-2)=u3(np,3:N-2)-k*(u3(nm,5:N)-4*u3(nm,4:N-1)+...
        6*u3(nm,3:N-2)-4*u3(nm,2:N-3)+u3(nm,1:N-4));


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

