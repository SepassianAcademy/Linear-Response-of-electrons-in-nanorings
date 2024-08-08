N=1000;
phi = linspace(-2,2,N);
rhom = [0.1,0.2,0.3,0.4,0.5,0.6,0.7];
%V_0 = [20,40,400,1000,10000];
V_0 = 20;
N_max = 10000;
L =0;
wavenumbers = zeros(length(rhom),length(phi),2);

for i=1:length(rhom)
rho = rhom(i);
  for l = 1:length(phi)
    current_L = L+phi(l);
    [kb_Ln] = kboundGS(current_L,V_0,rho,N_max);
    wavenumbers(i,l,1) = kb_Ln(1);
  end

end

hold on
for i=1:length(rhom)
    plot(phi,wavenumbers(i,:,1))
end

xlabel('Normalized flux $\Phi/\Phi_0$', 'Interpreter', 'latex','FontSize',16);
ylabel('Wavenumber $k_{LN}R$', 'Interpreter', 'latex','FontSize',16);
lgd = legend('$\rho = 0.1$','$\rho = 0.2$', '$\rho = 0.3$','$\rho = 0.4$','$\rho = 0.5$','$\rho = 0.6$', '$\rho = 0.7$','Interpreter', 'latex', 'Location', 'southeast','FontSize',10);
%lgd = legend('$V_0 = 20$', '$V_0 = 40$','$V_0 = 400$', '$V_0 = 1000$','$V_0 = 10000$','Interpreter', 'latex', 'Location', 'southeast','FontSize',10);
title('Bound-state wavenumbers, $V_0=20$.', 'Interpreter', 'latex','FontSize',20);
% Specify the full path for the SVG file
fullFilePath = '~/Documents/AAU/8.semester/P8/Results/resultsABring/wavenumberrho0102030507V_020phiN1000.svg';

% Save the plot as an SVG file
saveas(gcf, fullFilePath, 'svg');
file_path = '~/Documents/AAU/8.semester/P8/Results/resultsABring/wavenumberrho0102030507V_020phiN1000.csv';
writematrix(polvalabs, file_path);

