function wnhat=W0(xnhat,a,b,t,steps,wa,wb,wt,n,h)

wnhat=zeros(steps+1,1); 
intxnhat=zeros(1, steps+1); % integral of Xnhat(s,t) between a and b

if b>a+6
    for k=1:(steps+1)
        for i=1:2:(b-a-1)
            intxnhat(:,k)=intxnhat(:,k)+4*wt(:,a+i)'*xnhat{a+i}(:,k);
        end
        for j=2:2:(b-a-2)
            intxnhat(:,k)=intxnhat(:,k)+2*wt(:,a+j)'*xnhat{a+j}(:,k);
        end
    end
end


if b>a+6
    H=t(a+1)-t(a);
    for k=1:(steps+1)
        intxnhat(:,k)=intxnhat(:,k)+wt(:,a)'*xnhat{a}(:,k)+wt(:,b)'*xnhat{b}(:,k);
        intxnhat(:,k)=intxnhat(:,k)*H/3;
    end
end

for k=1:(steps+1)
    wnhat(k,1)=wb'*xnhat{b}(:,k)-wa'*xnhat{a}(:,k)-intxnhat(:,k);
    wnhat(k,1)=sqrt(n*h^3)*wnhat(k,1);
end
end