function lambda = QR_eig_AQ(A)
%find eigenvalues using QR factorization
for i=1:40
    [Q,R]=QR_fac_AQ(A);
    A=R*Q;
end
lambda=diag(A);