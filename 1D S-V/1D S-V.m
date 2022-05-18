% function 1D S-V
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.

J=44;
H_left=1.8;
H_right=0.6;
    
interv = 200;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.001;              % hardwired timestep
dx = 0.1;
filestep=1;
nplotstep = 2;           % plot interval
Urecord=zeros(1,interv+2);
Hrecord=zeros(1,interv+2);
Vrecord=zeros(1,interv+2);

% Initialize graphics
[surfplot,top,restart,quit] = initgraphics(interv,dx);

H = zeros(1,interv+2);
for j=1:J
    H(j)=H_left;
end
for j=(J+1):200
    H(j)=H_right;
end
nstep = 0;
n=0;

U = zeros(1,interv+2);   
Hx = zeros(1,interv+1); 
Term1 = zeros(1,interv+2);
Term2 = zeros(1,interv+2);
Ux = zeros(1,interv+1); 

while get(restart,'value')==0 && get(quit,'value')==0
    nstep = nstep + 1;    
    if max(U)>100
        j=1;
    end
    % Reflective boundary conditions
    H(1) = H(2);     
    U(1) = -U(2);  
    H(interv+2) = H(interv+1);  
    U(interv+2) = -U(interv+1);  

    % First half step
    % x direction
    i = 1:interv+1;
    % height
    Hx(i) = (H(i+1)+H(i))/2 - dt/(2*dx)*(U(i+1)-U(i));
    zero_index2=H>0;
    Term1(zero_index2)=U(zero_index2).^2./H(zero_index2);
       
    % x momentum
    Ux(i) = (U(i+1)+U(i))/2 -  ...
            dt/(2*dx)*(Term1(i+1) + g/2*(H(i+1)+ H(i)).*(H(i+1)-H(i))-Term1(i));
   
    % Second half step
    i = 2:interv+1;
    % height
    H(i) = H(i) - (dt/dx)*(Ux(i)-Ux(i-1));
    zero_index= H<=0;
    H(zero_index)=0;
    % x momentum
    zero_index=find(Hx>0);
    Term2(zero_index)=Ux(zero_index).^2./Hx(zero_index);
       
    U(i) = U(i) - (dt/dx)*(Term2(i) + g/2*(Hx(i)+ Hx(i-1)).*(Hx(i)-Hx(i-1)) - Term2(i-1));
    zero_index=find(H<0.01 & abs(U)<0.001);
    U(zero_index)=0;
    H(zero_index)=0;
  
    Velo=abs(U);
    max_velo=max(Velo);
    CFL=dt/dx*1.414*max_velo;
       
    if CFL>1
        CFL=1;
    end
    H(i) = (1-CFL)*H(i)+CFL*(H(i+1)+H(i-1))/2;
    U(i) = (1-CFL)*U(i)+CFL*(U(i+1)+U(i-1))/2;
       
    Urecord(nstep,:)=U;
    Hrecord(nstep,:)=H;
    Vrecord(nstep,:)=Urecord(nstep,:)./Hrecord(nstep,:);

    % Update plot
    if mod(nstep,nplotstep) == 0
        i = 2:interv+1;
        C = abs(U(zero_index));  % Color shows momemtum
        t = nstep*dt;
        tv = norm(C,'fro');
        set(surfplot,'ydata',H(i));
        set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',t,tv))
        drawnow
        if t==100%终止时间条件
            fileDir = 'E:\...\1D S-V data.mat';
            savePath = strcat(fileDir);
            save(savePath,'dt','dx','g','H_left','H_right','Hrecord','J','nstep','Urecord','Vrecord');
            break
        end
    end
    if all(all(isnan(H))), break, end  % Unstable, restart

clear
end
close(gcf)

function [surfplot,top,restart,quit] = initgraphics(n,dx);
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,restart,quit] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

   clf
   shg  %Show most recent graph window
   set(gcf,'numbertitle','off','name','Shallow_water')
   xend=(n-1)*dx;
   x = (0:dx:xend);
   surfplot = plot(x,x); 
   grid off
   axis([0 n*dx 0 2])
   caxis([-1 1])
   shading faceted
   c = (1:64)'/64;
   cyan = [0*c c c];
   colormap(cyan)
   top = title('xxx');
   restart = uicontrol('position',[120 20 80 20],'style','toggle','string','restart');
   quit = uicontrol('position',[20 20 80 20],'style','toggle','string','close');
end