
colors = lines(length(rhom));
hold on;
for i = 1:length(rhom)
    % Get the color for the current iteration
    currentColor = colors(i,:);
    
    % Plot with the current color
    plot(omegal(:), alphasol(i,:), 'Color', currentColor);
end

% Second plot loop
for i = 1:length(rhom)
    % Get the color for the current iteration
    currentColor = colors(i,:);
    
    % Plot with the current color
    plot(omegaval(:), polvalabs(i,:), 'o', 'MarkerSize', 3, 'Color', currentColor);
end
hold off;


legend('$\rho =0.2$, Analytical','$\rho =0.3$, Analytical', '$\rho =0.5$, Analytical','$\rho =0.7$, Analytical','$\rho =0.9$, Analytical','$\rho =0.2$, Numerical','$\rho =0.3$, Numerical', '$\rho =0.5$, Numerical','$\rho =0.7$, Numerical','$\rho =0.9$, Numerical', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('$\mathrm{Frequency} \ \omega R^2$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Polarizability} \ |\alpha(\Phi,\omega)|/R^4$', 'Interpreter', 'latex', 'FontSize', 18);
title('Dynamic Polarizability, $\Phi=0$', 'Interpreter', 'latex','FontSize',22);
fullFilePath = '~/Documents/AAU/8.semester/P8/Results/resultsABring/dypolnumanalrho0203050709PHi0.svg';
 
% Save the plot as an SVG file
saveas(gcf, fullFilePath, 'svg');fullFilePath = '~/Documents/AAU/8.semester/P8/Results/resultsABring/dynamicpolrho08090990999PHI0.svg';
 % Save the plot as an SVG file
saveas(gcf, fullFilePath, 'svg');