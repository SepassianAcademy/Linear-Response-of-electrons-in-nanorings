%Polarizationplots
L = @(l,phi) l+phi;
wl = @(L,k,x) bessely(L, k)*besselj(L, k*x)-bessely(L, k*x)*besselj(L,k);
 wlplus = @(L,k_Ln,x) 1*(bessely(L, k_Ln)*besselj(L+1, k_Ln*x)-bessely(L+1, k_Ln*x)*besselj(L,k_Ln));
 wlminus = @(L,k_Ln,x) -1*(bessely(L, k_Ln)*besselj(L-1, k_Ln*x)-bessely(L-1, k_Ln*x)*besselj(L,k_Ln));

%remember k_omega depends on l, n and omega. 
%remember that k_Ln depends on l and n.  

beta = @(L, k_Ln, rho) besselj(L, k_Ln)/(besselj(L,k_Ln*rho));

alphaplus1  = @(L,omega,komega,k_Ln,rho) k_Ln^2*komega*(beta(L,k_Ln,rho)^2*wlminus(L+1,komega,rho)+rho*wlplus(L,komega,rho)-4*beta(L,k_Ln,rho)/(pi*komega));
alphaplusd1 = @(L,omega,komega,k_Ln,rho) 4*omega^4*(beta(L,k_Ln,rho)^2-1)*rho*wl(L+1,komega,rho);
alphaplus2  = @(L,omega,komega,k_Ln,rho) k_Ln^2*komega*(beta(L,k_Ln,rho)^2*wlplus(L-1,komega,rho)+rho*wlminus(L,komega,rho)-4*beta(L,k_Ln,rho)/(pi*komega));
alphaplusd2 = @(L,omega,komega,k_Ln,rho) 4*omega^4*(beta(L,k_Ln,rho)^2-1)*rho*wl(L-1,komega,rho); 

alphaplus  = @(L,omega,komega,k_Ln,rho) -1/(2*omega^2)+ alphaplus1(L,omega,komega,k_Ln,rho)/(alphaplusd1(L,omega,komega,k_Ln,rho)) + alphaplus2(L,omega,komega,k_Ln,rho)/(alphaplusd2(L,omega,komega,k_Ln,rho));

omegal = linspace(0,3, 5000);
cphi = [0, 0.16, 0.32, 0.48];
alphasol = zeros(size(cphi,2), size(omegal,2));
l = 0;
rootnum=1;
rhom = [0.2,0.3,0.5,0.7,0.9];
%L = l+phi/phl = 0;
ftestval = zeros(1,size(omegal,2));
ftestval2 = zeros(1,size(omegal,2));

for i=1:size(cphi,2)
    for j=1:size(omegal,2)
        current_L = 0;
        func = @(k_Ln) wl(current_L,k_Ln,rho);
        k_Ln = fzero(func,3);
        komega = ((k_Ln)^2 + 2*(omegal(j)+0.05i))^(1/2);
        komega2 = ((k_Ln)^2 + 2*(-omegal(j)-0.05i))^(1/2);
        alphasol(i,j) = abs(alphaplus(current_L, (omegal(j)+0.05i), komega, k_Ln, rho)+alphaplus(current_L, (-omegal(j)-0.05i), komega2, k_Ln, rho));
        
    end
end
hold on
for i=1:size(cphi, 2)
    %Get the color for the current iteration

    %Plot with the current color
    plot(omegal(:), alphasol(i,:))
end
legend('$\Phi/\Phi_0 =0$','$\Phi/\Phi_0 =0.16$', '$\Phi/\Phi_0 =0.32$','$\Phi/\Phi_0 =0.48$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$\mathrm{Frequency} \ \omega R^2$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Polarizability} \ |\alpha(\Phi,\omega)|/R^4$', 'Interpreter', 'latex', 'FontSize', 18);
title('Dynamic Polarizability, $\rho=R_{inner}/R_{outer}=0.9$.', 'Interpreter', 'latex');
% fullFilePath = '~/Documents/AAU/8.semester/P8/Results/resultsABring/dynamicpolrho09test.svg';
% 
% % Save the plot as an SVG file
% saveas(gcf, fullFilePath, 'svg');
% 


saveas(gcf, fullFilePath);
