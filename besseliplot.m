% Define the range of r values
r = linspace(0, 1, 1000); % Adjust the range and number of points as needed

% Compute the function values
y = 1.1116*besseli(0, 1.6664*r);

% Plot the function
figure;
plot(r, y, 'LineWidth', 2);
xlabel('$r$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$1.6154 \cdot I_0(1.6664 \cdot r)$', 'Interpreter', 'latex', 'FontSize', 18);
title('Plot of $1.6154 \cdot I_0(1.6664 \cdot r)$', 'Interpreter', 'latex', 'FontSize', 22);
grid on;
