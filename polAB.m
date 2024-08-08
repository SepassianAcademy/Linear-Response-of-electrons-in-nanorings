omega_min = 0 + 10^(-5); % No change, confirmed correct
omega_max = 4;
N = 500;
G = 0.05;
rhom = [0.8,0.9,0.99];
phi = 0;
% rhom = [0.8, 0.9, 0.99, 0.999];
omegaval = linspace(omega_min, omega_max, N);
polval = zeros(length(rhom), length(omegaval));
Le = 1000;
N_max = 1000;

for l = 1:length(rhom)
    rho = rhom(l);
    L0 = 0; % Removed unnecessary 0+
    [kb_0n] = kboundAB(L0, rho, N_max, Le);
    k_01 = kb_0n(1);
    E_01 = k_01^2 / 2;
    kb_L1p = kboundAB(L0+1, rho, N_max, Le);
    kb_L1m = kboundAB(L0-1, rho, N_max, Le);
    polfunc = @(E_n, E_m, g_nm, omega) 2*g_nm /((E_m -E_n)^2-(omega + 1i*G)^2);
    for j = 1:length(omegaval)
        for k = 1:length(kb_L1p)
            k_Lm1m = kb_L1m(k);
            k_Lp1m = kb_L1p(k);
            [matep, matem] = matelement(L0, rho, k_01, k_Lm1m, k_Lp1m);
            g_nmp = abs(matep)^2 * (k_Lp1m^2-k_01^2) / 2;
            g_nmm = abs(matem)^2 * (k_Lm1m^2-k_01^2) / 2;
            E_mm = k_Lm1m^2/2;
            E_mp = k_Lp1m^2/2;
            polval(l, j) = polval(l, j) + polfunc(E_01, E_mm, g_nmm, omegaval(j))+polfunc(E_01, E_mp, g_nmp, omegaval(j));
        end
    end
end
hold on
hold on % Moved to before the plotting loop
polvalabs = abs(polval);
for l = 1:length(rhom)
    plot(omegaval, polvalabs(l, :))
end
%legend('$\rho =0.2$','$\rho =0.3$', '$\rho =0.5$','$\rho =0.7$','$\rho =0.9$', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$\mathrm{Frequency} \ \omega R^2$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Polarizability} \ |\alpha(\Phi,\omega)|/R^4$', 'Interpreter', 'latex', 'FontSize', 18);
title('Dynamic Polarizability, $\Phi=0$.', 'Interpreter', 'latex','FontSize',22);



% file_path = '~/Documents/AAU/8.semester/P8/Results/resultsABring/absalpharho0.2phi00.160.320.48N1000.csv';
% writematrix(polvalabs, file_path);

