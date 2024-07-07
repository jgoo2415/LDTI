function covhat0=cov0(cov2,a,b,steps,wa,wb)
covhat0=zeros(steps+1,steps+1);

for p=2:(steps+1)
    for k=2:(steps+1)
        covhat0(p,k)=wb'*cov2{b}(:,:,p,k)*wb+wa'*cov2{a}(:,:,p,k)*wa;        
    end
end
end