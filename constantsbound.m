function [N, B_0, C_0, D_0] = constantsbound(L,k,q,rho,r_max)
% L=0;
% k=4;
% V_0 = 10;
% rho = 0.8;
% r_max= 10;
% q = sqrt(-k^2+2*V_0);
 % Define relevant functions
W_0 = @(x) bessely(L,k)*besselj(L,k*x)-bessely(L,k*x)*besselj(L,k); 
W_8 = @(x) bessely(L+1,k*x)*besseli(L,q*x)*k+besseli(L+1,q*x)*bessely(L,k*x)*q;
W_9 = @(x) besselj(L+1,k*x)*besseli(L,q*x)*k+besseli(L+1,q*x)*besselj(L,k*x)*q;
W_10 = @(x) besselj(L+1,k*x)*bessely(L,k)-bessely(L+1,k*x)*besselj(L,k);

% Define the constants
B_1 = -pi*rho*W_8(rho)/2;
C_1 = pi*rho*W_9(rho)/2;
D_1 = pi*rho/(2*besselk(L,q))*(q*W_0(rho)*besseli(L+1,q*rho)+besseli(L,q*rho)*k*W_10(rho));
 % Relevant integrals
 integ1 = @(a,Z_Lm1, Z_L,Z_L1,Zp_Lm1, Zp_L,Zp_L1,r) 1/2 * r^2 * (Z_L*Zp_L - 1/2 * Z_Lm1*Zp_L1 - 1/2 * Z_L1*Zp_Lm1);



 % Calculating the normalization constant
 T_11 = @(r) integ1(q, besseli(L-1,r*q), besseli(L,r*q), besseli(L+1,r*q), besseli(L-1,r*q), besseli(L,r*q), besseli(L+1, r*q), r);
 S_21 = @(r) integ1(k, besselj(L-1, r*k), besselj(L, r*k), besselj(L+1, r*k), besselj(L-1, r*k), besselj(L, r*k), besselj(L+1, r*k), r);
 S_22 = @(r) integ1(k, besselj(L-1, r*k), besselj(L, r*k), besselj(L+1, r*k), bessely(L-1, r*k), bessely(L, r*k), bessely(L+1, r*k), r);
 S_23 = @(r) integ1(k, bessely(L-1, r*k), bessely(L, r*k), bessely(L+1, r*k), bessely(L-1, r*k), bessely(L, r*k), bessely(L+1, r*k), r);
 S_31 = @(r) integ1(q, besselk(L-1, r*q), besselk(L, r*q), besselk(L+1, r*q), besselk(L-1, r*q), besselk(L, r*q), besselk(L+1, r*q), r);

 g_1 = @(r) T_11(r);
 g_2 = @(r) B_1^2*S_21(r) + 2 * B_1*C_1*S_22(r)+C_1^2*S_23(r);
 g_3 = @(r) D_1^2*S_31(r);

 N = (g_3(r_max)-g_3(1)+g_2(1)-g_2(rho)+g_1(rho)-g_1(0))^(-1/2);
 B_0 = N*B_1;
 C_0 = N*C_1;
 D_0 = N*D_1;

end
