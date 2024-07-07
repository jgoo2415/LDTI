function [y]=init_template1b(nx,ny,nz,B)
n=nx*ny*nz;
mat=zeros(n,6);

for i=1:nx
    for j=1:ny
        for k=1:nz
            epsilon=0.05; 
            xx=(i/nx)*(1/0.5);
            yy=(j/ny)*(1/0.5);
            
            v1=sqrt(xx^2+yy^2);
            v2=(i/nx)^2/((0.5-epsilon)^2)+(j/ny)^2/((0.5-epsilon)^2);
            v3=(i/nx)^2/((0.5+epsilon)^2)+(j/ny)^2/((0.5+epsilon)^2);
  
          if (abs(k/nz-0.5)<epsilon)&&(v2>1)&&(v3<1)
            myv=[yy/v1 xx/v1 0; -xx/v1 yy/v1 0; 0 0 1];  
            mym=myv*diag([10,2,1])*myv';
            mat(nx*ny*(k-1)+nx*(j-1)+i,1)=mym(1,1);
            mat(nx*ny*(k-1)+nx*(j-1)+i,2)=mym(1,2);
            mat(nx*ny*(k-1)+nx*(j-1)+i,3)=mym(1,3);
            mat(nx*ny*(k-1)+nx*(j-1)+i,4)=mym(2,2);
            mat(nx*ny*(k-1)+nx*(j-1)+i,5)=mym(2,3);
            mat(nx*ny*(k-1)+nx*(j-1)+i,6)=mym(3,3);
 
          else
            mat(nx*ny*(k-1)+nx*(j-1)+i,1)=1;
            mat(nx*ny*(k-1)+nx*(j-1)+i,4)=1;
            mat(nx*ny*(k-1)+nx*(j-1)+i,6)=1;  
          end
        end
    end
end

y=zeros(n,48);
sigma=0.5*eye(48)+0.5*ones(48);
sigma2=sigma^(1/2);

for i=1:n
    z=randn(48,1);
    y(i,:)=B*mat(i,:)'+sigma2*z;
end

end

