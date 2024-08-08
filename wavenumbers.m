%Plotting the wavevector with boundary solutions
%l the l in paper
%phi is phi/phi_0 from the paper
numroots = 2;
L = @(l1,phi) l1+phi;
wl = @(L,k,x) bessely(L, k)*besselj(L, k*x)-bessely(L, k*x)*besselj(L,k);
l = [-1,0,1];


philin = linspace(0,0.5,1000);
kvalues = zeros(size(l,2), size(philin, 2), numroots);


%for i=1:(2*abs(l)+1)
for i=1:3
    %l1 = l(i)+i-1;
    for j=1:size(philin,2)
        for z=1:numroots
        current_L = L(l(i),philin(j));
        func = @(k) wl(current_L,k,0.2);
        kvalues(i,j,z) = fzero(func,3); %initial guess will either be 5 where z is 1 or 10 where z=2, this ensures that we get the correct root.
        end
    end
end
%  
figure;
hold on


% 
for i=1:3
    %i = 3;
     yplotvalues = kvalues(i, :,1);
     plot(philin, yplotvalues, 'red');
     yplotvalue1 = kvalues(i, :,2);
     plot(philin, yplotvalue1, 'blue');
end
% xlabel('Normalized flux $\Phi/\Phi_0$', 'Interpreter', 'latex');
% ylabel('Wavenumber $k_{LN}R$', 'Interpreter', 'latex');
% lgd = legend('$n = 1$', '$n = 2$','Interpreter', 'latex', 'Location', 'northeast');
% title('Bound-state wavenumbers as a function of normalized flux, $\rho=R_{inner}/R_{outer}=0.2$.', 'Interpreter', 'latex');
% % Specify the full path for the SVG file
% fullFilePath = '~/Documents/AAU/8.semester/P8/Results/resultsABring/wavenumberrho02.svg';
% 
% % Save the plot as an SVG file
% saveas(gcf, fullFilePath, 'svg');



saveas(gcf, fullFilePath);
