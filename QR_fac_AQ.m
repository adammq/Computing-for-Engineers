function [Q, R] = QR_fac_AQ(A)
%QR factorization
[r,c]=size(A);
if r~=c
    fprintf('input must be square matrix');
    Q = 'No Solution';
    R = 'No Solution';
    return
end
Q=eye(r);
R=A;

for i = 1:r-1
    %-------------------------c
    c=zeros(r,1);
    c(i:r,1)=R(i:r,i);
    %-------------------------e
    e=zeros(r,1);
    e(i,1)=1*sign(c(i,1));
    %-------------------------v
    v=c+norm(c)*e;
    %-------------------------H
    H=eye(r)-(2/(v'*v))*(v*v');
    %-------------------------Q
    Q=Q*H;
    %-------------------------R
    R=H*R;
end
end