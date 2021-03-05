%-- Record Movie

if (size(sol,1)==0)
  disp(['Run solver! ',char(7)]);
  return;
end

figure(1);
scaling = get(auto_ui,'Value');
f1=show(dx,plotvar,1,sol,tvec,u1_left,u2_left,u3_left,u1_right, ...
        u2_right,u3_right);   
if (scaling==0)
  xmin = str2num(get(xmin_ui,'String'));
  xmax = str2num(get(xmax_ui,'String'));
  ymin = str2num(get(ymin_ui,'String'));
  ymax = str2num(get(ymax_ui,'String'));
  if (size(xmin)==0 | size(xmax)==0 | size(ymax)==0 | size(ymin)==0)
    disp(['Give scaling factors!',char(7)]);
    return;
  end
else
  if (plotvar==4) % u2./u1
    xmin=-2; xmax=2; 
    ymin=min(sol(1,2,:)./sol(1,1,:));
    ymax=max(sol(1,2,:)./sol(1,1,:));
    for i=2:size(sol,1)
      ymin=min(ymin,min(sol(i,2,:)./sol(i,1,:)));
      ymax=max(ymax,max(sol(i,2,:)./sol(i,1,:)));
    end
  elseif (plotvar==5) % 0.4*u3-0.5*u2.^2./u1
    xmin=-2; xmax=2; 
    ymin=min(0.4*(sol(1,3,:)-0.5*sol(1,2,:).^2./sol(1,1,:)));
    ymax=max(0.4*(sol(1,3,:)-0.5*sol(1,2,:).^2./sol(1,1,:)));
    for i=2:size(sol,1)
      ymin=min(ymin,min(0.4*(sol(i,3,:)-0.5*sol(i,2,:).^2./sol(i,1,:))));
      ymax=max(ymax,max(0.4*(sol(i,3,:)-0.5*sol(i,2,:).^2./sol(i,1,:))));
    end
  else
    xmin=-2; xmax=2; 
    ymin=min(sol(1,plotvar,:));
    ymax=max(sol(1,plotvar,:));
    for i=2:size(sol,1)
      ymin=min(ymin,min(sol(i,plotvar,:)));
      ymax=max(ymax,max(sol(i,plotvar,:)));
    end
  end
end
maxis=[xmin xmax ymin ymax];      
lines=get(f1,'Children');
axis(maxis)
fmat=moviein(size(sol,1));
for i=1:size(sol,1)
  if (plotvar==1)
    set(lines,'Ydata',sol(i,plotvar,:));
  elseif (plotvar==2)
    set(lines,'Ydata',sol(i,plotvar,:));
  elseif (plotvar==3)
    set(lines,'Ydata',sol(i,plotvar,:));    
  elseif (plotvar==4)
    set(lines,'Ydata',sol(i,2,:)./sol(i,1,:));
  elseif (plotvar==5)  
    set(lines,'Ydata',0.4*(sol(i,3,:)-0.5*sol(i,2,:).^2./sol(i,1,:)));
  end
  fmat(:,i)=getframe(f1);
end
