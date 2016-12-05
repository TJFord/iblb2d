% distribution analysis
dx=5e-7;
dt=4.16e-8*2.5e5;
lx=fluid.lx
ly=fluid.ly;
msd=zeros(1,nts-1);
for it=1:nts-1
    for ip=1:np
        msd(1,it)=msd(1,it)+(yp(ip,it+1)-yp(ip,1))^2;
    end
    msd(1,it)=msd(1,it)/np;
end
figure;
time=[1:nts-1]*dt;
msd=msd*dx^2;
plot(time, msd(1,:),'b-','LineWidth',2)
% for i=1:np
%     idx = floor(yp(i,1)/binSize)+1;
%     dist(1,idx) = dist(1,idx)+1;
%     idx = floor(yp(i,floor(nts/2))/binSize)+1;
%     dist(2,idx) = dist(2,idx)+1;
%     idx = floor(yp(i,nts)/binSize)+1;
%     dist(3,idx) = dist(3,idx)+1;
% end
% figure;
% plot(dist(1,1:end),'b:','LineWidth',2)
% hold on
% plot(dist(2,1:end),'r-.','LineWidth',2)
% plot(dist(3,1:end),'g-','LineWidth',2)