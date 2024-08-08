
omega_min = 0+10^(-5);
omega_max = 3;
N = 500;
G = 0.05;
V_0 = 5000;
L0=0;
L1 = 1;
%rho = 0.9;
rhom = [0.8,0.9,0.99,0.999];
r_max = 10;
Le_max = 10;
N_max = 20; %Number of points for root finding in bound and unbound cases.

polval = zeros(length(rhom),N);
omegaval = linspace(omega_min,omega_max,N);
for s=1:length(rhom)
[qb_0n,kb_0n] = kbound(L0,V_0,rhom(s),N_max);
k_01 = kb_0n(1);
E_01 = k_01^2/2;
[qb_L1n, kb_L1n] = kbound(L1,V_0,rhom(s),N_max);
%[qu_L1n, ku_L1n] = findingkunbound(L1,V_0,rhom(s),r_max, N_max,Le_max);
froots1m = [kb_L1n];

polfunc = @(E_n, E_m, g_nm,omega) 2*g_nm/((E_m-E_n)^2-(omega+1i*G)^2);


  for j=1:length(omegaval)
     for k=1:length(froots1m)
        k_1m = froots1m(k);
         if k_1m^2 < 2*V_0
           q_1m = sqrt(2*V_0-k_1m^2); 
           q_01 = sqrt(2*V_0-k_01^2);
           [N0b, B0b, C0b, D0b] = constantsbound(0,k_01,q_01,rhom(s),r_max);
           [N1b, B1b, C1b, D1b] = constantsbound(1,k_1m,q_1m,rhom(s),r_max);
           matelb = matrixelementbound(k_01,k_1m,q_01,q_1m,r_max,rhom(s),N0b,B0b,C0b,D0b,N1b,B1b,C1b,D1b);
           g_nm = abs(matelb)^2*(k_1m^2-k_01^2);%1/2 not need because g_nm = 2*dipolematrixelement*energydifference in atomic units
           E_m = k_1m^2/2;
           polval(s,j) = polval(s,j) + polfunc(E_01,E_m,g_nm,omegaval(j));
         else 
           q_1m = sqrt(k_1m^2-2*V_0); 
           q_01 = sqrt(k_01^2-2*V_0);
           [N0b, B0b, C0b, D0b] = constantsbound(0,k_01,q_01,rhom(s),r_max);
           [N1ub, B1ub, C1ub, D1ub,E1ub] = constantsunbound(1,k_1m,q_1m,rhom(s),r_max);
           matelu = matrixelementunbound(k_01,k_1m,q_01,q_1m,r_max,rhom(s),N0b,B0b,C0b,D0b,N1ub,B1ub,C1ub,D1ub,E1ub);
           g_nm = abs(matelu)^2*(k_1m^2-k_01^2); %1/2 not need because g_nm = 2*dipolematrixelement*energydifference in atomic units
           E_m = k_1m^2/2;
           polval(s,j) = polval(s,j)+polfunc(E_01,E_m,g_nm,omegaval(j));
        end
     end
  end

end

polvalabs = abs(polval);

file_path = '~/Documents/AAU/8.semester/P8/Results/resultsABringfinitepotentialnumerical/Data/absalpharho0.80.90.990.999V_05000N500.csv';
% file_path2 = '~/Documents/AAU/8.semester/P8/Results/resultsABringfinitepotentialanalytical/Data/absalpharho0.20.30.50.70.9V_08.csv';
% Read the CSV file into a table
data_table = readtable(file_path2);

% Convert the table to a matrix
data_matrix = table2array(data_table);
hold on
plot(omegaval(:),polvalabs(1:length(rhom),:))
omegaval2 = linspace(omega_min,omega_max,1000);
%plot(omegaval2,data_matrix(1:length(rhom),:))

writematrix(polvalabs, file_path);

