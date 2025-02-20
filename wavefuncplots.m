
n = 1; L = 0; V_0m = [4,20,100]; omegam = [1.17,1.1313,0.97]; 
rhom =[0.1,0.5,0.9]; rho = 0.5; rp1 = linspace(0,rho,500*rho); 
rp2 =linspace(rho,1,500*(1-rho)); rp3 = linspace(1,3,1000); 
matp1 = zeros(3,length(rp1)); mat01 = zeros(3,length(rp1)); 
matp2 = zeros(3,length(rp2)); mat02 = zeros(3,length(rp2)); 
matp3 =zeros(3,length(rp3)); mat03 = zeros(3,length(rp2)); 

for l=1:length(V_0m)
rho = rhom(2); V_0 = V_0m(l); omega = omegam(l); 
[k, q] = findkandq(L,V_0, rho); 
[N, B0, C0, D0] = constantsunperturbed(L, k, q, rho); 
k_omega = sqrt(k^2 + 2*(vpa(omega, 50) + 0.05 * 1i)); 
s = sqrt(2*V_0 -k_omega^2); 
[C_11, C_21, C_22, C_3] = constantsperturbed(N, B0, C0, D0,k, k_omega, q, s, rho); 
phi_1 = @(r) besseli(1, s*r)*C_11 +N*(2*r*(q^2-s^2)*besseli(0, q*r)-4*q*besseli(1, q*r))/((q^2-s^2)^2);
phi_2 = @(r) C_21*besselj(1,k_omega*r)+C_22*bessely(1,k_omega*r)+B0*2*(r*(-k^2+k_omega^2)*besselj(0, k*r) + 2*k*besselj(1,k*r))/((k-k_omega)^2*(k+k_omega)^2)+...
C0*2*(r*(-k^2+k_omega^2)*bessely(0, k*r)+2*k*bessely(1,k*r))/((k-k_omega)^2*(k+k_omega)^2);

phi_3 = @(r) C_3*besselk(1,s*r)+2*D0*besselk(0, q*r)*r/(q^2-s^2) +...
4*D0*besselk(1, q*r)*q/(q^2-s^2)^2; phi_01 = @(r) N*besseli(0,q*r);
phi_02 = @(r) B0*besselj(0,k*r)+C0*bessely(0,k*r); phi_03 = @(r)D0*besselk(0,q*r);


for i=1:length(rp1)
    matp1(l,i) = abs(phi_1(rp1(i)))^2; mat01(l,i) =abs(phi_01(rp1(i)))^2;
end

for i=1:length(rp2)
    matp2(l,i) = abs(phi_2(rp2(i)))^2; mat02(l,i) = abs(phi_02(rp2(i)))^2;
end

for i=1:length(rp3)
    matp3(l,i) = abs(phi_3(rp3(i)))^2; mat03(l,i) =abs(phi_03(rp3(i)))^2;
end

end

%Define your colormap (e.g., lines)

colors = lines(length(V_0m));

% Initialize the figure and hold on to it
figure;
hold on;

% Loop through each set of data for l = 1:length(V_0m)
%     % Plot each set of data with the same color for the same l plot(rp1,
%     matp1(l, :), 'Color', colors(l, :)) %plot(rp1, mat01(l, :), '--',
%     'Color', colors(l, :)) plot(rp2, matp2(l, :), 'Color', colors(l, :))
%     %plot(rp2, mat02(l, :), '--', 'Color', colors(l, :)) plot(rp3,
%     matp3(l, :), 'Color', colors(l, :)) %plot(rp3, mat03(l, :), '--',
%     'Color', colors(l, :))
% end


% % Add dummy plots for the legend
% h(1) = plot(NaN, NaN, '-', 'Color', colors(1, :));
% h(2) = plot(NaN, NaN, '--', 'Color', colors(1, :));
% h(3) = plot(NaN, NaN, '-', 'Color', colors(2, :));
% %h(4) = plot(NaN, NaN, '--', 'Color', colors(2, :));
% %h(5) = plot(NaN, NaN, '-', 'Color', colors(3, :));
% %h(6) = plot(NaN, NaN, '--', 'Color', colors(3, :));
% Assuming the variables `rp1`, `rp2`, `rp3`, `matp1`, `matp2`, `matp3`, and `V_0m` are already defined

% Assuming the variables `rp1`, `rp2`, `rp3`, `matp1`, `matp2`, `matp3`, and `V_0m` are already defined
% Assuming the variables `rp1`, `rp2`, `rp3`, `matp1`, `matp2`, `matp3`, and `V_0m` are already defined

% Colors for the plots
colors = {'b', 'r', 'black'};

% Line thickness
lineThickness = 1.0;

% Create the figure and hold on
figure;
hold on;

% Create plot handles for the legend
h1 = plot(NaN, NaN, 'b', 'LineWidth', lineThickness);
h2 = plot(NaN, NaN, 'r', 'LineWidth', lineThickness);
h3 = plot(NaN, NaN, 'black', 'LineWidth', lineThickness);

for l = 1:length(V_0m)
    % Plot each set of data with the same color for the same l
    plot(rp1, matp1(l, :), 'Color', colors{l}, 'LineWidth', lineThickness);
    plot(rp2, matp2(l, :), 'Color', colors{l}, 'LineWidth', lineThickness);
    plot(rp3, matp3(l, :), 'Color', colors{l}, 'LineWidth', lineThickness);
end

% Labeling the axes with LaTeX interpreter
xlabel('$r/R$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$|\phi|^2$', 'Interpreter', 'latex', 'FontSize', 22);

% Adding the title with LaTeX interpreter
title('Radial wavefunctions, $\rho = 0.5$', 'Interpreter', 'latex', 'FontSize', 20);

% Adding the legend
legend([h1, h2, h3], ...
       {'$V_0 = 4$,\,$\omega = 1.17$', ...
        '$V_0 = 20$,\,$\omega = 1.13$', ...
        '$V_0 = 100$,\,$\omega = 0.97$'}, ...
       'Interpreter', 'latex', 'FontSize', 14);

% Ensure the hold is released
hold off;

