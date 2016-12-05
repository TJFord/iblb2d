% worm particle input 
clear
mcell=27e-12*1e-3 % kilo gram; 
dx = 0.5e-6; % m
% ks=1e-6; % N/m
% kb=4e-20; %
% kp=1e-13;
% ks=3e-5; % N/m
% kb=1.2e-18; %
% kp=1e-2;
% ks=1e-6;%1e-8; % N/m
% kb=1e-22;%8e-19; %
% kp=1e-2;
ks=0;%1e-8; % N/m
kb=0;%8e-19; %
kp=0;
kB=1.38e-23;%boltzmann const
d = 100e-9; % particle diameter
T = 300; % temperature in kelvin
vis = 1e-3; %pa*s
D = kB*T/(3*pi*vis*d); % diffusion
% lx=20e-6;
% ly=20e-6;
%ljcut potential parameter
epsilon=0.65e3/6.02e23;
sigma=0.16e-6;
rCut=0.5e-6;
lx=100;%300;%160;%40;
ly=50;%100;%40;

seed=762374;

Dc_p=3e-6;%6e-6; % polymer length % 8um
rho=2e3; % solid density: [kg/m^3];
acc=9.8; % gravity acceleration [m/S^2];

Nb=ceil(Dc_p/dx); % make lattice length and solid mesh to be the same
% Nb=4;
% if mod(Nb,2)~=0
%     Nb=Nb+1;
% end
% dtheta=2*pi/Nb; % 2h/L, angle increment
% index for left and right connected nodes
% m=mcell/Nb; % cell mass for each node
rhoWater=1e3; %kg/m^3;
vol=4/3*pi*(d/2)^3;
m = rhoWater*vol; % mass of water in a cubic lattice 
%% circle cell node position
% for k=0:(Nb-1)
%   theta=k*dtheta;
%   X(k+1,1)=(Dc_p/2)*cos(theta);
%   X(k+1,2)=(Dc_p/2)*sin(theta);  
% end
% %% biconcave shape, see Ref (2)
% for k=0:(Nb-1)
%   theta=k*dtheta;
%   X(k+1,1)=cos(theta);
% %   X(k+1,2)=sin(theta);
%   if sin(theta)>0
%     X(k+1,2)=0.5*sqrt(1-X(k+1,1)^2)*(0.207+2.002*X(k+1,1)^2-1.122*X(k+1,1)^4);
%   else
%     X(k+1,2)=-0.5*sqrt(1-X(k+1,1)^2)*(0.207+2.002*X(k+1,1)^2-1.122*X(k+1,1)^4);  
%   end
% end
% X(:,1)=(PHY.lx/2)+(Dc_p/2)*X(:,1);
% X(:,2)=(PHY.ly/2)+(Dc_p/2)*X(:,2);
% X(:,1)=(Dc_p/2)*X(:,1);
% X(:,2)=(Dc_p/2)*X(:,2);
% Xnew = [X;X(1,:)]/dx;
% Xnew(:,1) = Xnew(:,1) + lx/2;
% Xnew(:,2) = Xnew(:,2) + ly/10;
%% testing 4 points 
% X(1,1)=1;
% X(1,2)=25;
% X(2,1)=98;
% X(2,2)=25;
% X(3,1)=50;
% X(3,2)=1;
% X(4,1)=50;
% X(4,2)=48;
% for k=1:Nb
%     X(k,1)=50;
%     X(k,2)= 1+(k-1)*45;
% end
%% linear worm particle
for k=1:Nb
    X(k,1)=5;%30;
    X(k,2)=k + 5;
end
%% modify node position so that the equilibrium distance are equal
% pz=zeros(Nb,1);
% p=[X pz];
% p=Xnew;
% q = curvspace(p,Nb+1);
% X=q(1:end-1,1:2);
% 
% Xc = X;
X = X - 1;% convert it to c++ code starting with index 0
% X=X*0.9; % could be used to test if cell can restore back to equilibrium
% shape
% calculate geometric properties
% [L_e,beta_e,A_e]=Vec_TN_equilibrium(Xc',Nb);
%% new cell code
pos=X';
bond_list=zeros(3,Nb-1);
bond_list(1:2,:)=[1:Nb-1;2:Nb]; % first node, second node
r=pos(:,bond_list(1,:))-pos(:,bond_list(2,:)); % vector
bond_list(3,:)=sqrt(r(1,:).^2+r(2,:).^2); % length

% angle_list=[i,j,k,angle]', 3 nodes
angle_list=zeros(4,Nb-2);
angle_list(1:3,:)=[1:Nb-2;2:Nb-1;3:Nb]; % node index [left, middle, right]
dr1=pos(:,angle_list(1,:))-pos(:,angle_list(2,:)); % vector 1
dr2=pos(:,angle_list(3,:))-pos(:,angle_list(2,:)); % vector 2
rsq1=dr1(1,:).^2+dr1(2,:).^2;
r1=sqrt(rsq1);
rsq2=dr2(1,:).^2+dr2(2,:).^2;
r2=sqrt(rsq2);
c=dot(dr1,dr2); % dot product of two vectors 
c=c./(r1.*r2); % cosine value
angle_list(4,:)=acos(c); % angle value, different equilibrium value
% angle_list(4,:)=acos(-1.0); % Set all angles to be Pi

% index starting from 0
bond_list(1:2,:) = bond_list(1:2,:)-1;
angle_list(1:3,:)=angle_list(1:3,:)-1;

fid = fopen('MultiWorms.txt','a+');
fprintf(fid,'m');
fprintf(fid,'%20.10g\n',m);
fprintf(fid,'diameter');
fprintf(fid,'%20.10g\n',d);
fprintf(fid,'diffusion');
fprintf(fid,'%20.10g\n',D);
fprintf(fid,'ks');
fprintf(fid,'%20.10g\n',ks);
fprintf(fid,'kb');
fprintf(fid,'%20.10g\n',kb);
fprintf(fid,'kp');
fprintf(fid,'%20.10g\n',kp);
fprintf(fid,'kBT');
fprintf(fid,'%20.10g %10d\n',kB*T,seed);
fprintf(fid,'ljcut');
fprintf(fid,'%20.10g %20.10g %20.10g\n',epsilon,sigma,rCut);
fprintf(fid,'periodic %d %d\n',1,0);
% fprintf(fid,'periodic %d %d\n',1,1);
nCol=4;
nRow=6;%3;
nK=4;%2;
% nCol=1;
% nRow=1;
% nK=1;
ns = nK*nRow*nCol;
fprintf(fid,'ns\n');
fprintf(fid,'%10d\n',ns);
fprintf(fid,'node\n');
fprintf(fid,'%10d\n',ns*Nb);
Xfinal=[];
Xorg=X;
for k=1:nK   
    for j=1:nCol 
        for m = 1:nRow
            Xfinal=[Xfinal X];
            for i=1:Nb   
                fprintf(fid,'%20.15g \t',X(i,1));  
                fprintf(fid,'%20.15g \n',X(i,2)); 
            end 
%             X(:,2) = X(:,2)+14;
            X(:,2) = X(:,2)+7;
        end
        if mod(j,2)
            X(:,2) = Xorg(:,2)+3;
        else
            X(:,2) = Xorg(:,2);
        end
%         X(:,2) = Xorg(:,2);
        X(:,1) = X(:,1)+5;       
    end   
%     X(:,1) = X(:,1)+25;
    X(:,1) = X(:,1)+5;
end
fprintf(fid,'bond\n');
sz=size(bond_list);
fprintf(fid,'%10d\n',ns*sz(2));
for j=1:ns
    for i=1:sz(2)
        fprintf(fid,'%10d \t',bond_list(1,i));  
        fprintf(fid,'%10d \t',bond_list(2,i)); 
        fprintf(fid,'%20.15g\n',bond_list(3,i)); 
    end
    bond_list(1,:) = bond_list(1,:)+Nb;
    bond_list(2,:) = bond_list(2,:)+Nb;
end

fprintf(fid,'angle\n');
sz=size(angle_list);
fprintf(fid,'%10d\n',ns*sz(2));
for j=1:ns
    for i=1:sz(2)
        fprintf(fid,'%10d \t',angle_list(1,i));  
        fprintf(fid,'%10d \t',angle_list(2,i)); 
        fprintf(fid,'%10d \t',angle_list(3,i)); 
        fprintf(fid,'%20.15g\n',angle_list(4,i)); 
    end
    angle_list(1,:) =angle_list(1,:)+Nb;
    angle_list(2,:) =angle_list(2,:)+Nb;
    angle_list(3,:) =angle_list(3,:)+Nb;
end
fclose(fid);