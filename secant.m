function [x2,iter] = secant(func,x1,x2,tol)

%Evaluate function at x1
f1 = func(x1);
%Main Loop
dx = inf;
iter = 0;
while abs(dx) >= tol
    f2 = func(x2);
    dx = f2*(x2-x1)/(f2-f1);
    %Adjust points
    x1 =x2;
    f1 =f2;

    %Calculate New x2
    x2 = x2-dx;
    %Counter Iterations
    iter = iter+1;
end



end