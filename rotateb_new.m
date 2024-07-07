function [B]=rotateb_new(dir)
D=zeros(dir,3);
N=1;
a=4*pi/dir;
d=sqrt(a);
Mv=round(pi/d);
dv=pi/Mv;
dp=a/dv;
for m=0:(Mv-1)
 V=pi*(m+0.5)/Mv;
 Mp=round(2*pi*sin(V)/dp);
 for n=0:(Mp-1)
  P=2*pi*n/Mp;
  D(N,:)=[sin(V)*cos(P) sin(V)*sin(P) cos(V)];  
  N=N+1;
 end;
end;
B=zeros(dir,6);
for i=1:dir
 tmpmat=D(i,:)'*D(i,:);
 B(i,1) = tmpmat(1,1);
 B(i,2) = 2*tmpmat(1,2);
 B(i,3) = 2*tmpmat(1,3);
 B(i,4) = tmpmat(2,2);
 B(i,5) = 2*tmpmat(2,3);
 B(i,6) = tmpmat(3,3);
end;
B=-B; 
