function [Ux,Uy]=idx2xy(lx,ly,uf,vf,nts)
Ux=zeros(ly,lx,nts);
Uy=zeros(ly,lx,nts);
for k=1:nts
    for i=1:ly
        for j=1:lx
            Ux(i,j,k)=uf((i-1)*lx+j,k);
            Uy(i,j,k)=vf((i-1)*lx+j,k);
        end
    end
end