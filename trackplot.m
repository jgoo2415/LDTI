function y=trackplot(curve, mu, chat, alpha, n, h, steps, col)

r=sqrt(chi2inv(1-alpha,3));
phi=(0:0.05:1)*pi-pi/2;
theta=(0:0.05:1)*2*pi;
len=length(phi);

z=zeros(len^2,3);
for i=1:len
	for j=1:len
		z(len*(i-1)+j,1)=r*cos(phi(i))*cos(theta(j));
		z(len*(i-1)+j,2)=r*cos(phi(i))*sin(theta(j));
		z(len*(i-1)+j,3)=r*sin(phi(i));
	end
end
for i=2:steps
    ellipse=repmat((curve(:,i)-mu(:,i)/(sqrt(n*h^3))),1,len^2)-real(sqrtm(chat(:,:,i,i)))*z'./(sqrt(n*h^3));
 	plot3(ellipse(1,1:end),ellipse(2,1:end),ellipse(3,1:end), col);
    hold on
end
y=0;
plot3(curve(1,1:steps),curve(2,1:steps),curve(3,1:steps),'k','linewidth',1);
end


