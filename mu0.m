function muhat0=mu0(muhat,a,b,t,steps,wa,wb,wt)

muhat0=zeros(steps+1,1); 
intmuhat=zeros(1, steps+1); % integral of muhat(s,t) between a and b

if b>a+6
    for k=1:(steps+1)
        for i=1:2:(b-a-1)
            intmuhat(:,k)=intmuhat(:,k)+4*wt(:,a+i)'*muhat{a+i}(:,k);
        end
        for j=2:2:(b-a-2)
            intmuhat(:,k)=intmuhat(:,k)+2*wt(:,a+j)'*muhat{a+j}(:,k);
        end
    end
end


if b>a+6
    H=t(a+1)-t(a);
    for k=1:(steps+1)
        intmuhat(:,k)=intmuhat(:,k)+wt(:,a)'*muhat{a}(:,k)+wt(:,b)'*muhat{b}(:,k);
        intmuhat(:,k)=intmuhat(:,k)*H/3;
    end
end

for k=1:(steps+1)
    muhat0(k,1)=wb'*muhat{b}(:,k)-wa'*muhat{a}(:,k)-intmuhat(:,k);
end
end