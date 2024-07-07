%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Update: 07/07/2024
% Author: Juna Goo (junagoo@boisestate.edu)
% Reference: Goo, J., Sakhanenko, L. and Zhu, D. C. (2022) 
% A chi-square type test for time-invariant fiber pathways of the brain. 
% Statistical Inference for Stochastic Processes, 25, 449–469.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Example Code for time-dependent curves 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theorem 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(2024) % Control random number generator

nx=40; ny=40; nz=40; nt=19; n=nx*ny*nz*nt; mt=9;
t=(1:1:nt)/nt; % equally spaced time points

x0=[0.5*cos(0.3);0.5*sin(0.3);0.5]; % initial value
delta=0.015; steps=30; % Euler's method 

h=0.0167;  % bandwidth
beta=n*h^7; % beta
h1=(n/beta)^(-1/9);
h2=(n/beta)^(-1/10);


dir=48;% the number of gradient directions
B=rotateb_new(dir); % create a 48*6 B-matrix

y_array=cell(nt,1); % each cell, a (nx*ny*nz)*48 matrix; 
DD_array=cell(nt,1); % each cell, a 6*(nx*ny*nz) matrix; 

for ind=1:mt
    y_array{ind,1}=init_template1b(nx,ny,nz,B);  % generate the response vector (under the null hypothesis)
    DD_array{ind,1}=dtilda(y_array{ind,1},B);  % ols of D  
end

yvar=0.45;
for ind=(mt+1):nt
    y_array{ind,1}=init_template2b(nx,ny,nz,B,yvar);  % generate the response vector (under the alternative hypothesis)
    DD_array{ind,1}=dtilda(y_array{ind,1},B);  % ols of D  
end

v0=[1;1;1]; % an initial eigenvector (This can be modified depending on where the fiber trajectory needs to be traced.)

[xnhat, dnhat, dvdd, dddx, trH]=xnhat_all(n,nx,ny,nz,nt,DD_array,x0,t,delta,steps,h,h1,h2,v0); % Outputs are estimated integral curve, diffusion tensor kernel estimator, dv/dD, dD/du, and d^2D/du2. 


for ind=1:mt
    plot3(xnhat{ind}(1,:),xnhat{ind}(2,:),xnhat{ind}(3,:),'b') % estimated curves over the first half of time points
    hold on
end
for ind=(mt+1):nt
    plot3(xnhat{ind}(1,:),xnhat{ind}(2,:),xnhat{ind}(3,:),'r') % estimated curves over the second half of time points
    hold on
end

view([0 90]) % To project it onto a 2D-plane

parform=4; % number of workers for parallel computing
dnhat_array=dnhat_Xi(DD_array(:,1),n,nx,ny,nz,nt,h,parform); % To get dnhat(Xi) for all i=1,...,n. 
gam0gam0t=gamgamt_all_prep(n,nx,ny,nz,nt,dnhat_array,xnhat,t,steps,h,y_array(:,1),B); % sigma nhat 

[muhat,chat]=xnhatconf(dnhat,dvdd,dddx,gam0gam0t(:,1),trH,delta,steps,beta,nt); % estimated mean function and covariance function

ind=1; % at the 1st time point 
alpha=0.05; 
y=trackplot(xnhat{ind,1},muhat{ind,1},chat{ind,1},alpha,n,h,steps,'c'); % 95% confindence ellipsoids
hold on
ind=mt+1;
y=trackplot(xnhat{ind,1},muhat{ind,1},chat{ind,1},alpha,n,h,steps,'m'); % 95% confindence ellipsoids

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theorem 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chi_alt_svd=zeros(1,3); % To save the test statistic value under H0

a=1; b=19; 
wa=cell(1,3); wb=cell(1,3); wt=cell(1,3); % weight function examples

wa{1}=[t(a);t(a);t(a)]; wb{1}=[t(b);t(b);t(b)]; wt{1}=ones(3,nt);  % linear 
wa{2}=[exp(t(a));exp(t(a));exp(t(a))]; wb{2}=[exp(t(b));exp(t(b));exp(t(b))]; wt{2}=[exp(t);exp(t);exp(t)]; % exponential
wa{3}=ones(3,1); wb{3}=ones(3,1); wt{3}=zeros(3,nt); % constant

m=2; % the number of singular values used

for rep=1:3
    what=W0(xnhat,a,b,t,steps,wa{rep},wb{rep},wt{rep},n,h); 
    covhat=cov0(chat,a,b,steps,wa{rep},wb{rep}); 
    meanhat=mu0(muhat,a,b,t,steps,wa{rep},wb{rep},wt{rep});
    diff=what(2:steps+1,1)-meanhat(2:steps+1,1);
    P=covhat(2:steps+1,2:steps+1);
    [U S V]=svds(P,m); % truncated svd
    chi_alt_svd(1,rep)=diff'*pinv(U*S*V')*diff; % the observed value of the test statistic 
end 

chi_alt_svd 
chi2inv(0.95,m) % critical value to compare


