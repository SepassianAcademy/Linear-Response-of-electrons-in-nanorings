%findingk and q unbound
function [qu_Ln,ku_Ln] =  findingkunbound(L,V_0,rho,r_max, N_max,Le_max)

% 
%  V_0 = 8;
%  L = 1;
% rho = 0.9;
% r_max = 3;
% Le_max = 10;

%N_max = 1000;
tolerance1 = 1e-10;
tolerance2 = 1e-6;
tolerance3 = 1e-10;
s_max = sqrt(2*V_0)+Le_max;


W_4 = @(k,q,x) besselj(L+1, q*x)*bessely(L, k*x)*q-bessely(L+1,k*x)*besselj(L, q*x)*k;
W_5 = @(k,q,x) q*bessely(L,k*x)*bessely(L+1,q*x)-bessely(L+1,k*x)*k*bessely(L,q*x);
W_6 = @(k,q,x) besselj(L+1,k*x)*besselj(L,q*x)*k-besselj(L+1,q*x)*besselj(L,k*x)*q;
W_7 = @(k,q,x) besselj(L, k*x)*bessely(L+1,q*x)*q-bessely(L,q*x)*besselj(L+1,k*x)*k;


f = @(k,q) ((W_5(k,q,1)*besselj(L,k*rho)-W_7(k,q,1)*bessely(L,k*rho))*besselj(L,q*r_max)-bessely(L,q*r_max)*(W_4(k,q,1)*besselj(L,k*rho) + ...
    bessely(L,k*rho)*W_6(k,q,1)))*q*besselj(L+1,q*rho)-besselj(L,q*rho)*k*((W_5(k,q,1)*besselj(L,q*r_max)-bessely(L,q*r_max)*W_4(k,q,1))*besselj(L+1,k*rho) + ...
    (-besselj(L,q*r_max)*W_7(k,q,1)-bessely(L,q*r_max)*W_6(k,q,1))*bessely(L+1, k*rho));

q = @(k) sqrt(k^2-2*V_0);

% Define the anonymous function for root finding
root_func = @(k) f(k, q(k));

% Create a vector of values for k
k = linspace(sqrt(2*V_0)+10^(-6), s_max, N_max*10);  % Reduced range for visualization


roots = [];
for i = (1:length(k)-1)
      [root,iter] = secant(root_func,k(i),k(i+1), tolerance1);
      roots = [roots; real(root)]; % Append to roots1
end

roots = sort(roots);
%Remove NaN values from roots
roots = roots(~isnan(roots));
filtered_roots = [];

% Filter the roots based on the given condition
for j = 1:length(roots)
    if (roots(j) < s_max) && (roots(j) > sqrt(2 * V_0)) && (abs(root_func(roots(j)))<tolerance3)
        filtered_roots = [filtered_roots, roots(j)];
    end
end

% Initialize filtered_roots2 with the first element of filtered_roots
if ~isempty(filtered_roots)
    filtered_roots2 = filtered_roots(1);

    % Filter out roots that are within the specified tolerance
    for i = 2:length(filtered_roots)
        if all(abs(filtered_roots(i) - filtered_roots2) > tolerance2)
            filtered_roots2 = [filtered_roots2; filtered_roots(i)];
        end
    end
else
    filtered_roots2 = [];
end




ku_Ln = filtered_roots2;
qu_Ln = (ku_Ln.^2-2*V_0).^(1/2);


