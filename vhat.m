%v is (3,1)-vec, mat is (6,1)-vec for 3,3 symm matrix; 11, 12, 13, 22, 23, 33;
function [v,lambda]=vhat(m)
[vec,d]=eig([m(1,1) m(2,1) m(3,1); m(2,1) m(4,1) m(5,1); m(3,1) m(5,1) m(6,1)]);
[lam, Ind]=sort([d(1,1) d(2,2) d(3,3)]);
v_eig=vec(:,Ind);
lambda=lam(3);
v=v_eig(:,3)./norm(v_eig(:,3));
end