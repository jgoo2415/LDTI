function gam0gam0t=gamgamt_all_prep(n,nx,ny,nz,nt,dhat_array,xnhat,t,steps,h,y_array,B)

gam0gam0t=cell(nt,1); % save only a and b
for ind=1:nt
    gam0gam0t{ind}=zeros([6,6, steps+1]);        
end
   
BBinv=inv(B'*B);
factor=1/(((2*pi)^2)*n*(h^4));


for q0=1:nt
    for s=1:(steps+1)
        u=xnhat{q0}(:,s);
        sig=zeros(48,48);
        for q=1:nt
            for i=1:nx
                for j=1:ny
                    for k=1:nz
                        inner=((u(1)-i/nx)^2+(u(2)-j/ny)^2+(u(3)-k/nz)^2+(t(q0)-q/nt)^2);
                        if inner<=16*h^2
                            kernel=exp(-0.5*(inner/h/h));
                   
                            yy=(y_array{q}(nx*ny*(k-1)+nx*(j-1)+i, :)'-B*dhat_array{q}(:,nx*ny*(k-1)+nx*(j-1)+i));
                            sig=sig+kernel*(yy*yy');
                        end
                    end
                end
            end
        end
        gam0gam0t{q0}(:,:,s)=(BBinv)*B'*sig*B*(BBinv)'*factor;
    end    
end


end 

