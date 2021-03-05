function f1=show(dx,var,frame,uu,tvec,u1_left,u2_left,u3_left,u1_right,u2_right,u3_right);

  if (length(uu)==0) 
    
    % Initial conditions
    N=round(4/dx)+1;    
    for i=1:N
      x=dx*(i-1)-2;
      if x>0
        u1(i)=u1_right;
        u2(i)=u2_right;
        u3(i)=u3_right;
      else
        u1(i)=u1_left;
        u2(i)=u2_left;
        u3(i)=u3_left;
      end
    end
    x=-2:dx:2;
    if (var==1)
      f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
      plot(x,u1);title('Density,  t = 0');
    elseif (var==2)
      f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
      plot(x,u2);title('Momentum,  t = 0');
    elseif (var==3)
      f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
      plot(x,u3);title('Energy,  t = 0');
    elseif (var==4)
      f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
      plot(x,u2./u1);title('Velocity,  t = 0');
    elseif (var==5)
      f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
      plot(x,0.4*u3-0.5*u2.^2./u1);
      title('Pressure,  t = 0');
    end

    return;
  end

  t=tvec(frame);
  x=-2:dx:2;
  temp=zeros(1,length(x));
  N=round(4/dx)+1;
  if (var==1)
    f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
    temp(:)=uu(frame,1,:);
    plot(x,temp);title(['Density,  t = ' num2str(t)]);
  elseif (var==2)
    f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
    temp(:)=uu(frame,2,:);
    plot(x,temp);title(['Momentum,  t = ' num2str(t)]);
  elseif (var==3)
    f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
    temp(:)=uu(frame,3,:);
    plot(x,temp);title(['Energy,  t = ' num2str(t)]);
  elseif (var==4)
    f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
    temp(:)=uu(frame,2,:)./uu(frame,1,:);
    plot(x,temp);title(['Velocity,  t = ' num2str(t)]);
  elseif (var==5)
    f1=subplot(2,1,1);set(f1,'Position',[0.13 0.3 0.775 0.6 ]);
    temp(:)=0.4*(uu(frame,3,:)-0.5*uu(frame,2,:).^2./uu(frame,1,:));
    plot(x,temp);
    title(['Pressure,  t = ' num2str(t)]);
  end
  