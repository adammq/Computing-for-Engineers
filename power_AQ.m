function [lambda, V] = power_AQ(A,X0,tol,max_iter)
%Input - A is an nxn matrix
%      - X0 is the nx1 starting vector
%      - tol is the tolerance in norm of the difference between two
%        successive calculated eigenvectors
%      - max_iter is the maximum number of iterations
%
%Output - lambda is the dominant eigenvalue
%       - V is the eigenvector corresponding to the dominant eigenvalue

for i = 1:max_iter
    Y=A*X0;
    [~,index]=max(abs(Y));
    l=Y(index);
    X=Y/l;
    
    if norm(X-X0) < tol
        fprintf("Tolerance met after %d iterations\n",i);
        lambda=l;
        V=X;
        return
    end
    
    X0=X;
end

fprintf("Tolerance not met.\n");
lambda="No Solution";
V="No Solution";
end