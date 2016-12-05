function [f]=readGeom(input)
file = input;
fid = fopen(file);
while  ~feof(fid)
    next=fgetl(fid);
    if strcmp(next,'lx')
        lx = str2num(fgetl(fid));
    elseif strcmp(next,'ly')
        ly =str2num( fgetl(fid));
    elseif strcmp(next,'nf')   
        nf =str2num( fgetl(fid));
    elseif strcmp(next,'nbList')
        nbList = zeros(nf,11);
        for i=1:nf
            nbList(i,:)=str2num(fgetl(fid));           
        end
    elseif strcmp(next,'boundaries')
        nbc =str2num( fgetl(fid));
    elseif strcmp(next,'velocity')
        nv =str2num( fgetl(fid));
        if(nv>0)
            v_size = str2num(fgetl(fid));
            v_norm = str2num(fgetl(fid));
            vel=zeros(v_size,3);
            for i=1:v_size
                vel(i,:)=str2num(fgetl(fid));
            end
        end
    elseif strcmp(next,'open')  
        no = str2num(fgetl(fid));
        if(no>0)
            o_size =str2num( fgetl(fid));
            o_norm = str2num(fgetl(fid));
            open=zeros(o_size,2);
            for i=1:o_size
                open(i,:)=str2num(fgetl(fid));
            end
        end
    elseif strcmp(next,'bounceback') 
        nb =str2num(fgetl(fid));
        if (nb>0)
            bb_size = str2num(fgetl(fid));
            bb=zeros(1,bb_size);
            for i=1:bb_size
                bb(i)=str2num(fgetl(fid));
            end
        end
    end
    
end
fclose(fid);

    f.lx=lx;
    f.ly=ly;
    f.nf=nf;
    f.nbList = nbList;
    f.nbc=nbc; 
    f.nv=nv;
    if(nv)     
        f.v_norm=v_norm;
        f.vel=vel;
    end
    f.no=no;
    if(no)
        f.o_norm=o_norm;
        f.open=open;
    end
    f.nb=nb;
    if(nb) 
        f.bb=bb;
    end
