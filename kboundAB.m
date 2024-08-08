%Findingk bound
function [kb_Ln] = kboundAB(L,rho,N_max,Le)
% V_0 = 2;
% L = 0;
% rho = 0.9;% r_max = 10;
% N_max = 1000;
tolerance2 = 10^(-8);

W =  @(k,x)  besselj(L, k*x)*bessely(L, k) - bessely(L, k*x)*besselj(L, k);


func = @(k)(W(k, rho));


initial_guesses1 = linspace(0.001, Le, N_max); % Initial guesses for roots1


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

%end
