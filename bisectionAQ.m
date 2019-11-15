function x_root = bisectionAQ(x_est,tol,max_iter)
%NEWTONAQ makes use of Newton's method to approximate where deflection is
%maximum along a beam subjected to a triangular distributed load.
%%
%assigned values to variables
E=70e9;   %Pa
I=52.9e-6; %m^(-4)
w=20e3;   %N/m
L=4;      %m

%define deflection function and its derivative
F=@(x) w*x/(360*E*I*L)*(7*L^4-10*L^2*x^2+3*x^4); %deflection function
f=@(x) w/(360*E*I*L)*(7*L^4-30*L^2*x^2+15*x^4);  %derivative function

%define x interval and corresponding outputs
x_bound=[0,L];
f_bound=[f(x_bound(1)),f(x_bound(2))];

%define x vector
x=zeros(1,max_iter);

%initialize x(1)
x(1)=x_est;

%%
for i=1:max_iter
    if i==1
        continue;
    end
    
    %making use of some math logic to speedily tighten the interval - same
    %sign means on same side of interval
    if f(x(i-1))*f_bound(1)>0
        x_bound(1)=x(i-1);
        f_bound(1)=f(x_bound(1));
    else
        x_bound(2)=x(i-1);
        f_bound(2)=f(x_bound(2));
    end
    
    %generating next approximation
    x(i)=(x_bound(1)+x_bound(2))/2;
    
    %calculating estimated relative error
    est_rel_err(i-1)=abs(x(i)-x(i-1))/x(i);
    
    %checking to see if our tolerence has been met
    if est_rel_err(i-1) <= tol
        
        %store approximated root, print relevant metrics
        x_root = x(i);
        fprintf("Tolerance reached in %i iterations.\n",i);
        fprintf("Root approximated as x = %fm, with deflection of %fm.\n",x(i),F(x(i)));
        fprintf("Generating convergence plot ...\n");
        
        %generate Convergence Plot
        semilogy(est_rel_err);
        title("Convergence Plot");
        xlabel("Iterations");
        ylabel("Estimated Relative Error");
        
        return;
    end
end

%dealing with the exceptions ...
if max_iter==1 || (max_iter>=1 & est_rel_err(end)>tol)
    fprintf("Tolerance not met - need more iterations.\n");
    x_root = NaN;
end

end %function bisectionAQ