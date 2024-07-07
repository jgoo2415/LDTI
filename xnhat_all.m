function [xnhat, dnhat, dvdd, dddx, trH]=xnhat_all(n,nx,ny,nz,nt,DD_array,x0,t,delta,steps,h,h1,h2,v0)

xnhat=cell(nt,1); % each cell, a 3*(steps+1) matrix
dnhat=cell(nt,1); % each cell, a 6*(steps+1) matrix
dvdd=cell(nt,1); % each cell, a 3*6*(steps+1) tensor
dddx=cell(nt,1); % each cell, a 6*3*(steps+1) tensor
trH=cell(nt,1); % each cell, a 6*1*(steps+1) tensor

for q0=1:nt
    xnhat{q0,1}=zeros(3, steps+1);
    xnhat{q0,1}(:,1)=x0; % regardless of time t, set initial values 
    
    dnhat{q0,1}=zeros(6, steps+1);
    dvdd{q0,1}=zeros(3,6,steps+1);
    dddx{q0,1}=zeros(6,3,steps+1);
    trH{q0,1}=zeros(6,1,steps+1);
end

for q0=1:nt
    v2=v0;
    for k=1:steps
        [mat, dmdu, dm2du2]=dnhat_all(DD_array,n,nx,ny,nz,nt,[xnhat{q0,1}(:,k); t(q0)],h,h1,h2); 
        dnhat{q0,1}(:,k)=mat;
    
        [v1,lambda]=vhat(mat); 
        if v2'*v1<0
            v1=-v1;
        end
        
        dvdd{q0,1}(:,:,k)=dvdm(mat,v1,lambda); % gv is (3,6)
        dddx{q0,1}(:,:,k)=dmdu(:,1:3);
        trH{q0,1}(:,:,k)=dm2du2; % trH is (6,1) 

        xnhat{q0,1}(:,k+1)=xnhat{q0,1}(:,k)+delta*v1;
        v2=v1;
    end
    
   [mat, dmdu, dm2du2]=dnhat_all(DD_array,n,nx,ny,nz,nt,[xnhat{q0,1}(:,steps+1); t(q0)],h,h1,h2);  
   dnhat{q0,1}(:,steps+1)=mat;
   
   [v1,lambda]=vhat(mat); 
   if v2'*v1<0
       v1=-v1;
   end
   dvdd{q0,1}(:,:,steps+1)=dvdm(mat,v1,lambda); % gv is (3,6)
   dddx{q0,1}(:,:,steps+1)=dmdu(:,1:3);
   trH{q0,1}(:,:,steps+1)=dm2du2; % trH is (6,1) 
end
end


