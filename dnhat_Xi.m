function [mat]=dnhat_Xi(DD_array,n,nx,ny,nz,nt,h,parform)

mat=cell(nt,1);  % dnhat of all U_i, i=1,..., n
% each cell of mat is (6,nx*ny*nz);

for q0=1:nt
    mat{q0}=zeros(6, nx*ny*nz);
end

%u is i0/nx, j0/ny, k0/nz, q0/nt
cut1=16*h^2;
cut2=1/(2*h^2);
factor=1/(((2*pi)^2)*n*(h^4));

parpool('local',parform);
parfor (q0=1:nt, parform)
    for i0=1:nx
        for j0=1:ny
            for k0=1:nz

                u=[i0/nx, j0/ny, k0/nz, q0/nt];
                ind0=nx*ny*(k0-1)+nx*(j0-1)+i0;
                for q=1:nt
                    for i=1:nx
                        for j=1:ny
                            for k=1:nz                            
                                inner=((u(1)-i/nx)^2+(u(2)-j/ny)^2+(u(3)-k/nz)^2+(u(4)-q/nt)^2);
                                if inner<=cut1
                                    kernel=exp(-cut2*inner);
                                    ind=nx*ny*(k-1)+nx*(j-1)+i;
                                    mat{q0}(:,ind0)=mat{q0}(:,ind0)+kernel*DD_array{q}(:,ind);                                     
                                end
                            end
                        end
                    end
                end

           
            end
        end
    end
 mat{q0}=mat{q0}*factor; % nh^(d+1)  
end
delete(gcp('nocreate'))
end