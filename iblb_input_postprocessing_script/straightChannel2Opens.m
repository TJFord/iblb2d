% gimp, draw 2d pic, doesn't work for many obstacles
% image/mode/gray/ windows bmp save
clear
% a=imread('shape2.bmp');
% a=imread('p20.bmp');
% a=imread('cylinder.bmp');
% a=a/255;
% a=round(a);
clear
% a=imread('shape2.bmp');
% a=imread('p20.bmp');
% a=imread('cylinder3.bmp');
% a=imread('cylinder2.bmp');
% a=imread('vonKarman.bmp');
% a=imread('channel.bmp');
a=imread('40x100.bmp');
% a=imread('cylinder.bmp');
% BW1 = edge(a,'sobel');
% imshow(BW1)
% BW = edge(a,'canny');
% figure, imshow(BW2)
% BW = edge(a,'log');
% figure, imshow(BW3)
BW = edge(a,'roberts');
% figure, imshow(BW4)
a = BW;

[Rows,Cols]=size(a);

[x,y]=meshgrid(1:Rows,1:Cols);

xmin=Rows;
xmax=0;
ymin=Cols;
ymax=0;

b=zeros(Rows,Cols);
wall=0;
%% find xmin, xmax, ymin, ymax
for i=1:Rows
    for j=1:Cols
        if a(i,j)==1            
            if xmin>i, xmin = i;end
            if xmax<i, xmax = i;end
            if ymin>j, ymin = j;end
            if ymax<j, ymax = j;end
             b(i,j)=a(i,j); % save boundary flag 1
        end                    
    end
end

for i=1:Rows   
        index=(find(a(i,:)==1));
        if (~isempty(index))
            if length(index) == (index(end)-index(1)+1)
                wall = 1; %straight wall
            else
                wall = 0; % nonstraight wall
            end 
    
            if wall % if there is a continual wall, fill before and after with -1
                if (index(1)-1)>0 % region before wall
                    b(i,1:index(1)-1)=-1;
                end
                if (Cols-index(end))>0 % region after wall
                    b(i,index(end)+1:Cols)=-1;
                end

            else
                for j=1:Cols
                    if a(i,j)==0
                        c = 1; % assuming it is outside
                        pre = find(a(i,1:j)==1);
                        pst = find(a(i,j:end)==1);
                        if isempty(pre) | isempty(pst)
                            c = 1; % no previous bc or no successing bc
                        else
                            sz1=length(pre);
                            sz2=length(pst);
                            n1=find(diff(pre)~=1); % check if boundary is continuous, countinous means only one boundary
                            n2=find(diff(pst)~=1); %
                            if isempty(n1) | isempty(n2) % if continuous means it is inside, for porous media, it doesn't work
                                c=0;
%                             elseif ~mod(n1,2) | ~mod (n2,2)
%                                 c=0;
                            end
                        end; 
                        
%                         pre = find(a(i))
%                         c=1;% outside indicator
%                             
%                             id1=find(a(i,k:end)==1);
%                             if ~isempty(id1)
%                                 sz = length(id1);                    
%                                 if sz == 1, c = 0;end % inside
%                                 if sz == id1(end)-id1(1)+1, c=0;end
%                             end
                        if c, b(i,j)=-1;end
                    end
                end
            end
        else
            b(i,:)=-1;
        end
  
end
%% write out simple geometry, for checking purpose
fid = fopen('geom.txt','a+');
%fprintf(fid,'Time step %d\n',time);
fprintf(fid,'%10d',xmax-xmin+1,ymax-ymin+1);
fprintf(fid,'\n');
for i=xmin:xmax
    for j=ymin:ymax
        fprintf(fid,'%2d ',b(i,j));        
    end 
    fprintf(fid,'\n');
end
fclose(fid);


%% prepare input for c++ code
c= b(xmin:xmax,ymin:ymax);
Rows = xmax-xmin+1;
Cols = ymax-ymin+1;
nf = Rows*Cols-length(find(c==-1)); % total active fluid nodes, including boundary
nt = nf + 1; % extra one to fill missing adjacent nodes
index = zeros(Rows,Cols); %[id, i, j]
id = 0;
for i=1:Rows
    for j=1:Cols
        if c(i,j) ~= -1
            id=id+1;
            index(i,j) = id;
        end
    end
end
cx = [  0,   1,  0, -1,  0,    1,  -1,  -1,   1];   % unit vectors: x-comp
cy = [  0,   0,  1,  0, -1,    1,   1,  -1,  -1]; 
% find adjacent nodes
% [nid, 8 adjacent nodes, coordinates x,y]
nlist = zeros(nf,11);
id = 0;
for i=1:Rows% rows->y
    for j=1:Cols % cols->x
        if index(i,j)>0
            id=id+1;
            nlist(id,1)=id;
            nlist(id,10:11)=[j,i];%[i,j];
            for k=2:9
%                 ii=i+cx(k);
%                 jj=j+cy(k);   
                ii=i+cy(k); % y
                jj=j+cx(k); % x
                if ii==0 | jj==0 | ii==Rows+1 | jj==Cols+1 % out of array bounds
                    nd=nt;
                elseif index(ii,jj)~=0 % if there is a adjacent node
                    nd=index(ii,jj); 
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
dx = 1e-2;% m
tau=0.5072;% no units
rho = 1e3;%kg/m^3
dm = rho*dx^3; %kg
vis = 1e-6;% kinetic viscosity m^2/s
uMax = 5e-4;% m/s
index = index -1;
nlist = nlist -1;

fid = fopen('channel.txt','a+');
fprintf(fid,'scheme bgk\n');
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
%% pay attention to this 
fprintf(fid,'%10d\n',Cols);
fprintf(fid,'ly');
fprintf(fid,'%10d\n',Rows);
fprintf(fid,'nf');
fprintf(fid,'%10d\n',nf);
fprintf(fid,'nt');
fprintf(fid,'%10d\n',nt);
fprintf(fid,'nbList\n');
for i=1:nf
    for j=1:11
        fprintf(fid,'%10d',nlist(i,j));   
    end
    fprintf(fid,'\n');
end
% IJidx matrix
fprintf(fid,'IJidx\n');
for i=1:Rows
    for j=1:Cols
        fprintf(fid,'%10d',index(i,j));
    end
    fprintf(fid,'\n');
end
% bounceback node list
inidx=0;
outidx=0;
bbidx=0;
inlet=zeros(1,Rows)-1;
inlet_neigb=zeros(1,Rows)-1;
outlet=zeros(1,Rows)-1;
outlet_neigb=zeros(1,Rows)-1;
for i=1:Rows
    for j=1:Cols
        if j==1 & i >1 & i < Rows % inlet node
            if c(i,j)==1
                inidx=inidx+1;
                inlet(inidx)=index(i,j);
%                 inlet_neigb(inidx)=index(i,j+1);
                inlet_neigb(inidx)=index(i,Cols);
            end
        elseif j==Cols & i >1 & i < Rows % outlet node
            if c(i,j)==1
                outidx=outidx+1;
                outlet(outidx)=index(i,j);
%                 outlet_neigb(outidx)=index(i,j-1);% left neighbor node id for right open bc 
                outlet_neigb(outidx)=index(i,1);% left neighbor node id for right open bc 
            end
        else %other bounce back node
            if c(i,j)==1
                bbidx=bbidx+1;
                bb(bbidx)=index(i,j);
            end
        end
    end
end


inlet = inlet(find(inlet>-1));

sz=length(inlet);
% parabolic velocity for inlet
u=zeros(sz,2);
%% can be set by physical inputs
% uMax = 0.02;% this can be specified by physical value
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
% do it on purpose to test multiple velocity boundareis;

VS(1)=Vel;
VS(2)=Vel;


outlet = outlet(find(outlet>-1));
outlet_neigb = outlet_neigb(find(outlet_neigb>-1));
Open.name='open';
Open.norm=[1,0];
Open.node=outlet;
Open.data = outlet_neigb;

inlet = inlet(find(inlet>-1));
inlet_neigb = inlet_neigb(find(inlet_neigb>-1));

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

%bcS(3,1)=struct('name',{},'norm',{},'node',{},'data',{});
% bcS(1)=Vel;
% bcS(2)=Open;
% bcS(3)=BBk;

% write boundary conditions
% fprintf(fid,'boundaries');
% fprintf(fid,'%10d\n',nbc);
%% write velocity data
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
        fprintf(fid,'%20.10f',VS(i).data(j,1));
        fprintf(fid,'%20.10f',VS(i).data(j,2));
        fprintf(fid,'\n');   
    end
    fprintf(fid,'\n');   
end
%% write open boundary
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
% for i=1:nv
%     fprintf(fid,'%s',bcS(i).name);
%     fprintf(fid,'%10d',length(bcS(i).node));
%     if(~isempty(bcS(i).norm))
%         fprintf(fid,'%10d',bcS(i).norm);
%     end
%     fprintf(fid,'\n');
%     for j=1:length(bcS(i).node)
%         fprintf(fid,'%10d',bcS(i).node(j));
%         if(~isempty(bcS(i).data))
%             fprintf(fid,'%10f',bcS(i).data(j,1));
%             fprintf(fid,'%10f',bcS(i).data(j,2));
%         end
%         fprintf(fid,'\n');
%     end    
% end

% for i=1:nbc
%     fprintf(fid,'%s',bcS(i).name);
%     fprintf(fid,'%10d',length(bcS(i).node));
%     if(~isempty(bcS(i).norm))
%         fprintf(fid,'%10d',bcS(i).norm);
%     end
%     fprintf(fid,'\n');
%     for j=1:length(bcS(i).node)
%         fprintf(fid,'%10d',bcS(i).node(j));
%         if(~isempty(bcS(i).data))
%             fprintf(fid,'%10f',bcS(i).data(j,1));
%             fprintf(fid,'%10f',bcS(i).data(j,2));
%         end
%         fprintf(fid,'\n');
%     end    
% end

fclose(fid);