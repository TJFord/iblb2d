% distribution analysis
lx=fluid.lx
ly=fluid.ly;
binSize=1;
nbin=floor(ly/binSize)+1;
nStep=2:1:nts;
sz=length(nStep);
dist=zeros(sz+1,nbin);

for j=1:sz+1
    for i=1:np
        if j==1
            idx = floor(yp(i,1)/binSize)+1;
            dist(j,idx) = dist(j,idx)+1;
        else
            idx = floor(yp(i,j)/binSize)+1;
            dist(j,idx) = dist(j,idx)+1;
        end

    %     idx = floor(yp(i,1)/binSize)+1;
    %     dist(1,idx) = dist(1,idx)+1;
    %     idx = floor(yp(i,floor(nts/2))/binSize)+1;
    %     dist(2,idx) = dist(2,idx)+1;
    %     idx = floor(yp(i,nts)/binSize)+1;
    %     dist(3,idx) = dist(3,idx)+1;
    end
end
% figure;
for j=1:sz+1
    figure;
    plot(dist(j,1:end),'b-','LineWidth',2);
    % hold on
end
% plot(dist(1,1:end),'b:','LineWidth',2)
% plot(dist(2,1:end),'r-.','LineWidth',2)
% plot(dist(3,1:end),'g-','LineWidth',2)

%% smoothing data
distSm=zeros(sz+1,nbin);
% smooth over space
% for i=1:sz+1
% %     distSm(i,:)=smooth(dist(i,1:end),3);
%     distSm(i,:)=smooth(dist(i,1:end),3);
% end
% smooth over time
for i=1:nbin
%     distSm(i,:)=smooth(dist(i,1:end),3);
    distSm(:,i)=smooth(dist(:,i),3);
end
for j=1:sz+1
    figure;
    plot(distSm(j,1:end),'b-','LineWidth',2);
    axis([0 50 0 30])
% hold on
end
% figure;
% plot(distSm(1,1:end),'b:','LineWidth',2)
% hold on
% plot(distSm(2,1:end),'r-.','LineWidth',2)
% plot(distSm(3,1:end),'g-','LineWidth',2)