function [mu, chat]=xnhatconf(dnhat, dvdd, dddx, gamgam, trH, delta, steps, beta, nt)

mu=cell(nt,1);
chat=cell(nt,1);

for q0=1:nt
    mu{q0,1}=zeros(3, steps+1);
    chat{q0,1}=zeros(3,3, steps+1, steps+1);
end

factor1=0.5*sqrt(beta);
factor2=1/(8*pi*sqrt(pi));

for q0=1:nt
    mat=dnhat{q0,1}(:,2);
    sig=gamgam{q0,1}(:,:,2);
    covpart=dvdd{q0,1}(:,:,2)*(mat*mat'+sig)*dvdd{q0,1}(:,:,2)';
    mu2=dvdd{q0,1}(:,:,2)*trH{q0,1}(:,:,2);
    
    mu{q0,1}(:,2)=delta*factor1*mu2;
    chat{q0,1}(:,:,2,2)=delta*factor2*covpart;
    
    for p=2:steps  
        mu1=dvdd{q0,1}(:,:,p)*dddx{q0,1}(:,:,p)*mu{q0,1}(:,p); % mu1 is (3,1)
	    mu2=dvdd{q0,1}(:,:,p)*trH{q0,1}(:,:,p);	% mu2 is (3,1)
    
        mat=dnhat{q0,1}(:,p);
        sig=gamgam{q0,1}(:,:,p);
        covpart=dvdd{q0,1}(:,:,p)*(mat*mat'+sig)*dvdd{q0,1}(:,:,p)'; % covpart is (3,3)
        covpart2=dvdd{q0,1}(:,:,p)*dddx{q0,1}(:,:,p); % covpart2 is (3,3)
        
        mu{q0,1}(:, p+1)=mu{q0,1}(:,p)+delta*mu1+delta*factor1*mu2;
        chat{q0,1}(:,:,p+1,p+1)=chat{q0,1}(:,:,p,p)+delta*factor2*covpart+delta*covpart2*chat{q0,1}(:,:,p,p)+...
        delta*chat{q0,1}(:,:,p,p)*covpart2';   
    end
end


for q0=1:nt
for p=2:steps
    for k=1:p
        covpart1=dvdd{q0,1}(:,:,p)*dddx{q0,1}(:,:,p); 
        chat{q0,1}(:,:,p+1,k)=chat{q0,1}(:,:,p,k)+delta*covpart1*chat{q0,1}(:,:,p,k);
        chat{q0,1}(:,:,k,p+1)=chat{q0,1}(:,:,p+1,k);
    end
end
end



end