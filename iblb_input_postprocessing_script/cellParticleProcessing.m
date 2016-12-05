clear;
% nf = 400; % fluid nodes 
fluidExist=1;
solidExist=1;
chainExist =1;
nts =50; % time steps
ncells=21;
nForCell=52;
if fluidExist
    fgeom='./fgeom.txt';
    fluid=readGeom(fgeom);
    load fluidRst.txt;
    load fluidForce.txt
    nf = fluid.nf;
    nv=fluid.nv;
    no=fluid.no;
    nb=fluid.nb;

    nbList=fluid.nbList;

    uf=zeros(nf,nts);
    vf=zeros(nf,nts);
    fx_f=zeros(nf,nts);
    fy_f=zeros(nf,nts);
    xf=zeros(nf,1);
    yf=zeros(nf,1);
    usq_f = zeros(nf,nts);
    for i=1:nts
        for j=1:nf
            uf(j,i)=fluidRst((i-1)*nf+j,1);
            vf(j,i)=fluidRst((i-1)*nf+j,2);
            usq(j,i)=sqrt(uf(j,i)^2+vf(j,i)^2);
            fx_f(j,i)=fluidForce((i-1)*nf+j,1);
            fy_f(j,i)=fluidForce((i-1)*nf+j,2);
        end
    end
    if(nb)
        bb=fluid.bb;
        uf(bb+1,:)=0;
        vf(bb+1,:)=0;
    end
% 
    for i=1:nf
        xf(i,1)=nbList(i,10)+1; % convert index from c to matlab
        yf(i,1)=nbList(i,11)+1;
    end

end

% convert it into matrix form
[Ux,Uy]=idx2xy(fluid.lx,fluid.ly,uf,vf,nts);
UU=Ux.^2+Uy.^2;


if solidExist 
    ns =nForCell*ncells;%8;%52;%64;


    load cellRst.txt;
    load cellForce.txt;
%     load cellRef.txt;
    load cellVelocity.txt;
    
    xs=zeros(ns,nts);
    ys=zeros(ns,nts);
    xsR=zeros(ns,nts);
    ysR=zeros(ns,nts);
    fx=zeros(ns,nts);
    fy=zeros(ns,nts);
    us=zeros(ns,nts);
    vs=zeros(ns,nts);
    for i=1:nts  
        for j=1:ns
            xs(j,i)=cellRst((i-1)*ns+j,1)+1;
            ys(j,i)=cellRst((i-1)*ns+j,2)+1;
%             xsR(j,i)=cellRef((i-1)*ns+j,1)+1;
%             ysR(j,i)=cellRef((i-1)*ns+j,2)+1;
            fx(j,i)=cellForce((i-1)*ns+j,1);
            fy(j,i)=cellForce((i-1)*ns+j,2);
            us(j,i)=cellVelocity((i-1)*ns+j,1);
            vs(j,i)=cellVelocity((i-1)*ns+j,2);
        end
    end
end
% % [x,y]=meshgrid(1:n,1:n);
% 

if chainExist 
    np =378;%12*48;%12*24;%52;%64;


    load chainRst.txt;
    load chainForce.txt;
    load chainVelocity.txt;
    
    xp=zeros(np,nts);
    yp=zeros(np,nts);

    fx_p=zeros(np,nts);
    fy_p=zeros(np,nts);
    up=zeros(np,nts);
    vp=zeros(np,nts);
    for i=1:nts  
        for j=1:np
            xp(j,i)=chainRst((i-1)*np+j,1)+1;
            yp(j,i)=chainRst((i-1)*np+j,2)+1;
%             xsR(j,i)=cellRef((i-1)*ns+j,1)+1;
%             ysR(j,i)=cellRef((i-1)*ns+j,2)+1;
            fx_p(j,i)=chainForce((i-1)*np+j,1);
            fy_p(j,i)=chainForce((i-1)*np+j,2);
            up(j,i)=chainVelocity((i-1)*np+j,1);
            vp(j,i)=chainVelocity((i-1)*np+j,2);
        end
%         for k=1:96
%             for m=1:6
%                 j=(k-1)*6+m;
%                 xp(j,i)=chainRst((i-1)*np+j,1)+1;
%                 yp(j,i)=chainRst((i-1)*np+j,2)+1;
%     %             xsR(j,i)=cellRef((i-1)*ns+j,1)+1;
%     %             ysR(j,i)=cellRef((i-1)*ns+j,2)+1;
%                 fx_p(j,i)=chainForce((i-1)*np+j,1);
%                 fy_p(j,i)=chainForce((i-1)*np+j,2);
%     %             us(j,i)=cellVelocity((i-1)*ns+j,1);
%     %             vs(j,i)=cellVelocity((i-1)*ns+j,2);
%             end
%         end
    end
end

% write a movie
% vidObj=VideoWriter('cellFlow.avi');

mov(1:nts) = struct('cdata', [],...
                        'colormap', []);
vidObj=VideoWriter('movie.avi');
vidObj.FrameRate = 5;
open(vidObj);
figure;
for it=1:nts
%     figure;
clf
if fluidExist
%     quiver(xf,yf,uf(:,it),vf(:,it));
    hold on
%     quiver(xf,yf,fx_f(:,it),fy_f(:,it));
    imagesc(UU(:,:,it));
    axis equal
end
if solidExist
    % imagesc(xf,yf,usq(:,it)); 
%     for i=1:ncells %number of cells
%             plot(xs((i-1)*nForCell+1:i*nForCell,it),ys((i-1)*nForCell+1:i*nForCell,it),'or-','MarkerFaceColor',[0,1,0],'MarkerSize',2);
%     end
%     hold on
    plot(xs(:,it),ys(:,it),'or','MarkerFaceColor',[0,1,0],'MarkerSize',2);
%     plot(xs(1,it),ys(1,it),'ko','MarkerFaceColor',[1,0,0]);
%     plot(xsR(:,it),ysR(:,it),'*k-');
%     quiver(xs(:,it),ys(:,it),fx(:,it),fy(:,it),'color','b');
    axis equal
end
if chainExist
    % imagesc(xf,yf,usq(:,it)); 
    hold on
    
    plot(xp(:,it),yp(:,it),'ob','MarkerFaceColor',[1,0,0]);
%     plot(xs(1,it),ys(1,it),'ko','MarkerFaceColor',[1,0,0]);
%     plot(xsR(:,it),ysR(:,it),'*k-');
%     quiver(xp(:,it),yp(:,it),fx_p(:,it),fy_p(:,it),'color','g');
%     quiver(xp(:,it),yp(:,it),up(:,it),vp(:,it),'color','r');
    axis equal
% for i=1:96
% %     for j=1:6
% %         m=(i-1)*6+j;
%         plot(xp((i-1)*6+1:i*6,it),yp((i-1)*6+1:i*6,it),'og-','MarkerFaceColor',[1,0,0],'MarkerSize',2);
% %     end
% end
%     plot(xs(1,it),ys(1,it),'ko','MarkerFaceColor',[1,0,0]);
%     plot(xsR(:,it),ysR(:,it),'*k-');
%     quiver(xp(:,it),yp(:,it),fx_p(:,it),fy_p(:,it),'color','g');
    axis equal
end
% axis([0 40 0 40])
% axis([0 160 0 160])
% axis([60 120 60 100])
FF(it) = getframe;
mov(it)=getframe(gcf);
% writeVideo(vidObj,FF(it));
end
movie(FF,1)
movie2avi(mov, 'movie2.avi', 'compression', 'None','fps',2);
close(vidObj);


% %% center position
% xc=zeros(nts,2);
% for i=1:nts
%     for j=1:ns
%         xc(i,1) = xc(i,1)+xsR(j,i);
%         xc(i,2) = xc(i,2)+ysR(j,i);
%     end
%     xc(i,:) = xc(i,:)/ns;
% end
% vc=zeros(nts-1,2);
% for i=1:nts-1
%     if xc(i+1,1)<xc(i,1)
%         vc(i,:)=xc(i+1,:)-xc(i,:)+401-15;    
%     else    
%         vc(i,:)=xc(i+1,:)-xc(i,:);       
%     end
% end
% dx=1e-4;
% dt=1.667e-3;%2.4e-5;
% rho = 1e3; %kg/m^3
% ts=1e3;%1e4;%
% nu = 0.166667;%0.0024;
% vsim = max(vc(:,1))*dx/(dt*ts)
% 
% g=9.8;
% D=10*dx;
% drho=1e-3*rho;;%1e-4*rho;%0.001;
% drho=drho*ns/3.14/5^2;% ns*drho*g*area = rho'*g*pi*r^2
% vis=1e-3;
% vt = g*D^2*drho/18/vis % eqn(3.6) Numerical simulations of particle sedimentation us-ing the immersed boundary method
% 
% % eqn(3.10)
% 
% % options=optimset('Algorithm','levenberg-marquardt');   % Option to display output
% % [x,fval] = fsolve(@findRe,0.01,options)  % Call solver
% % Re=fsolve(@findRe,0.01)
% Re = vsim*D*rho/vis;
% Cd = 8*pi/Re/log(7.4/Re);
% vt2=sqrt(pi*g*D*drho/2/Cd/rho)
% 
% % immersed fem method, cd = 0.5*rho_w*v^2*area
% Fd = drho*g*pi*D^3/6;
% vt3=sqrt(Fd/(0.5*rho*pi*D^2/4))
% 
% figure;
% maxT=35
% vphy=[0;vc(:,1)]*dx/(dt*ts);
% xt=[0:maxT];
% xt=xt*ts*dt;
% axes('FontSize',14)
% plot(xt,vphy(1:maxT+1),'-','LineWidth',2)
% hold on
% plot([0 maxT]*dt*ts,[vt,vt],'--','LineWidth',2)
% xlabel('\fontsize{14} Time(s)');
% ylabel('\fontsize{14} Velocity(m/s)');
% legend('\fontsize{14} Simulation','\fontsize{14} Stoke prediction')
