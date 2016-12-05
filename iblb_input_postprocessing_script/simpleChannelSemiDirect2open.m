% a simple square geometry with semi direct indexing
clear
lx= 100;
ly =40;
a=zeros(ly,lx);
a(1,:)=1;
a(end,:)=1;
a(:,1)=1;
a(:,end)=1;

index =zeros(ly,lx);
for i=1:ly % y
    for j=1:lx % x
        index(i,j)= (i-1)*lx+j; 
    end
end

cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];   % unit vectors: x-comp
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]; 
% find adjacent nodes
% [nid, 8 adjacent nodes, coordinates x,y]
nf = lx*ly;%n*n;
% lx=n;
% ly=n;
nt = nf+1;
nlist = zeros(nf,11);
id = 0;
for i=1:ly % y
    for j=1:lx %x
        if index(i,j)>0
            id=id+1;
            nlist(id,1)=id; % first one start from 0
            nlist(id,10:11)=[j,i];
            for k=2:9
                x=j+cx(k);
                y=i+cy(k);   
                if x==0 || y==0 || x==lx+1 || y==ly+1 % out of array bounds
                    nd=nt;
                elseif index(y,x)~=0 % if there is a adjacent node
                    nd=index(y,x); 
                else % if index(ii,jj)==0, out of boundary
                    nd=nt;
                end
                nlist(id,k)=nd;
            end
        end
    end
end

% water parameters , unit length 1 lattice, unit time 1 lb time step
% since lb fluid density = 1, we get unit mass rho_water*dx^3
dx = 1e-6;% m
tau=1;% no units
rho = 1e3;%kg/m^3
dm = rho*dx^3; %kg
vis = 1e-6;%m^2/s
uMax = 1000e-6;% m/s

index = index - 1;
nlist = nlist - 1; % c index starts from 0

% fid = fopen('input.txt','a+');
fid = fopen('inputChannel.txt','a+');
fprintf(fid,'dx');
fprintf(fid,'%20.10f\n',dx);
fprintf(fid,'tau');
fprintf(fid,'%20.10f\n',tau);
fprintf(fid,'rho');
fprintf(fid,'%20.10f\n',rho);
fprintf(fid,'vis');
fprintf(fid,'%20.10f\n',vis);
fprintf(fid,'u');
fprintf(fid,'%20.10f\n',uMax);
fprintf(fid,'lx');
fprintf(fid,'%10d\n',lx);
fprintf(fid,'ly');
fprintf(fid,'%10d\n',ly);
fprintf(fid,'nf');
fprintf(fid,'%10d',nf);
fprintf(fid,'\n');
fprintf(fid,'nbList\n');
for i=1:nf
    for j=1:11
        fprintf(fid,'%10d',nlist(i,j));   
    end
    fprintf(fid,'\n');
end
% IJidx matrix
fprintf(fid,'IJidx\n');
for i=1:ly
    for j=1:lx
        fprintf(fid,'%10d',index(i,j));
    end
    fprintf(fid,'\n');
end
% bounceback node list
inidx=0;
outidx=0;
bbidx=0;
inlet=zeros(1,ly)-1;
outlet=zeros(1,ly)-1;
outlet_neigb=zeros(1,ly)-1;
inlet_neigb=zeros(1,ly)-1;
for i=1:ly
    for j=1:lx
        if j==1 & i>1 & i< ly % inlet node            
            inidx=inidx+1;
            inlet(inidx)=index(i,j);
            inlet_neigb(inidx)=index(i,j+1);
        elseif j==lx & i>1 & i< ly% outlet node         
            outidx=outidx+1;
            outlet(outidx)=index(i,j);
            outlet_neigb(outidx)=index(i,j-1);
        else %other bounce back node
            if a(i,j)==1
                bbidx=bbidx+1;
                bb(bbidx)=index(i,j);
            end
        end
    end
end

inlet = inlet(find(inlet>-1));
sz=length(inlet);
u=zeros(sz,2);
uMax=0.02;
ndir=[-1 0];
L=max(nlist(inlet+1,11));% y
L0=min(nlist(inlet+1,11));
for i=1:sz
    y=nlist(inlet(i)+1,11)-L0;
    u(i,1)=4*uMax/(L*L)*(L*y - y*y);
    u(i,2)=0;
end
nbc = 3; % inlet, outelt and bounce back boudaries
fprintf(fid,'boundaries');
fprintf(fid,'%10d\n',nbc);
Vel.name='velocity';
Vel.norm = [-1,0];
Vel.node = inlet;
Vel.data = u;
% Vel.data = inlet_neigb;
% do it on purpose to test multiple velocity boundareis;

VS(1)=Vel;
% VS(2)=Vel;

outlet = outlet(find(outlet>-1));
outlet_neigb = outlet_neigb(find(outlet_neigb>-1));
Open.name='open';
Open.norm=[1,0];
Open.node=outlet;
Open.data = outlet_neigb;


Open2.name='open';
Open2.norm=[-1,0];
Open2.node=inlet;
Open2.data = inlet_neigb;
OP(1)=Open;
OP(2)=Open2;

bb=bb(find(bb>-1));
BBk.name='bounceback';
BBk.norm=[];
BBk.node=bb;
BBk.data=[];

% write velocity data
nv=0;% 2 means two set of velocities
fprintf(fid,'velocity');
fprintf(fid,'%10d',nv);
% for i=1:nv
%     fprintf(fid,'%10d',length(VS(i).node));
% %     fprintf(fid,'%10d',VS(i).norm);
% end
fprintf(fid,'\n');
for i=1:nv
    sz = length(VS(i).node);
    fprintf(fid,'%10d\n',sz);
    fprintf(fid,'%10d',VS(i).norm(1));
    fprintf(fid,'%10d\n',VS(i).norm(2));
    
    for j=1:sz
        fprintf(fid,'%10d',VS(i).node(j));
        fprintf(fid,'%10f',VS(i).data(j,1));
        fprintf(fid,'%10f',VS(i).data(j,2));
        fprintf(fid,'\n');   
    end
    fprintf(fid,'\n');   
end

%% write open boundary
% sz = length(Open.node);
no = 2;% 1;
fprintf(fid,'open');
fprintf(fid,'%10d',no);
fprintf(fid,'\n');
for i=1:no
    sz = length(OP(i).node);
    fprintf(fid,'%10d\n',sz);
    fprintf(fid,'%10d',OP(i).norm(1));
    fprintf(fid,'%10d',OP(i).norm(2));
    fprintf(fid,'\n');
    for j=1:sz
        fprintf(fid,'%10d',OP(i).node(j));
        fprintf(fid,'%10d',OP(i).data(j));
        fprintf(fid,'\n');   
    end
    fprintf(fid,'\n');   
end


%% write bounceback boundary
sz = length(BBk.node);
fprintf(fid,'bounceback');
fprintf(fid,'%10d',sz);
fprintf(fid,'\n');
for j=1:sz
    fprintf(fid,'%10d',BBk.node(j));
end
fprintf(fid,'\n');  
fclose(fid);