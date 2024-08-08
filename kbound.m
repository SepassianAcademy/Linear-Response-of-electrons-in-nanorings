%Findingk bound
function [qb_Ln,kb_Ln] = kbound(L,V_0,rho,N_max)
% V_0 = 2;
% L = 0;
% rho = 0.9;% r_max = 10;
% N_max = 1000;
tolerance2 = 10^(-8);

W =  @(k,x)  besselj(L, k*x)*bessely(L, k) - bessely(L, k*x)*besselj(L, k);
W1 = @(k,x) bessely(L + 1, k)*besselj(L, k*x) - bessely(L, k*x)*besselj(L + 1, k);
W2 = @(k,q,x) bessely(L + 1, k*x)*k* besselk(L, q) - q*bessely(L, k)*besselk(L + 1, q*x);

W3 = @(k,q,x) k*besselk(L, q)*besselj(L + 1, k*x) - q*besselj(L, k)*besselk(L + 1, q*x);

q = @(k) sqrt(2*V_0-k^2);

func = @(k)(-q(k)*W(k, rho)*besselk(L + 1, q(k)) + k*besselk(L, q(k))*W1(k, rho))*q(k)*besseli(L + 1, q(k)*rho) +....W(k, rho)
    besseli(L, q(k)*rho)*k*(W2(k, q(k), 1)*besselj(L + 1, k*rho) - bessely(L + 1, k*rho)*W3(k, q(k), 1));


initial_guesses1 = linspace(0.001, sqrt(2*V_0), N_max); % Initial guesses for roots1


roots1=[];
for i = 1:length(initial_guesses1)
    try
        root = fzero(func, initial_guesses1(i));
        if ~isnan(root)
             roots1 = [roots1;root]; % add to roots1 the found root
        end
    catch
        % Do nothing if fzero fails
    end
end

% Initialize filtered_roots2 with the first element of filtered_roots
filtered_roots2 = roots1(1);

    % Filter out roots that are within the specified tolerance
for i = 2:length(roots1)
        if all(abs(roots1(i) - filtered_roots2) > tolerance2)
            filtered_roots2 = [filtered_roots2; roots1(i)];
        end
end



kb_Ln = filtered_roots2;
qb_Ln = (2*V_0-kb_Ln.^2).^(1/2);

%end
