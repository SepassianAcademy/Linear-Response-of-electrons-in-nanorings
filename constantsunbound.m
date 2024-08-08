function [N,B,C,D,E] = constantsunbound(L,k,q,rho,r_max)
% L=10;
% k=20;
% V_0 = 30;
% q = sqrt(k^2-2*V_0);
% rho = 0.8;
% r_max= 10;
% 
% q = sqrt(k^2-2*V_0);

%Relevant functions
W_0 = @(x) bessely(L,k)*besselj(L,k*x)-bessely(L,k*x)*besselj(L,k);
W_1 = @(x) bessely(L+1,k)*besselj(L,k*x)-bessely(L,k*x)*besselj(L+1, k);
W_4 = @(x) besselj(L + 1, q*x)*bessely(L, k*x)*q - bessely(L + 1, k*x)*besselj(L, q*x)*k;
W_5 = @(x) q*bessely(L,k*x)*bessely(L+1,q*x)-bessely(L+1,k*x)*k*bessely(L,q*x);
W_6 = @(x) besselj(L+1,k*x)*besselj(L,q*x)*k-besselj(L+1,q*x)*besselj(L,k*x)*q;
W_7 = @(x) besselj(L, k*x)*bessely(L+1,q*x)*q-bessely(L,q*x)*besselj(L+1,k*x)*k;


%Relevant integrals 
integ1 = @(a,Z_Lm1, Z_L,Z_L1,Zp_Lm1, Zp_L,Zp_L1,r) 1/2*r^2*(Z_L*Zp_L-1/2*Z_Lm1*Zp_L1-1/2*Z_L1*Zp_Lm1);

B1 = pi*rho*W_4(rho)/2;
C1 = pi*rho*W_6(rho)/2;

D1 = -1/4*pi^2*rho*(-(q*W_0(rho)*bessely(L+1,q)-k*bessely(L, q)*W_1(rho))*q*besselj(L+1,q*rho) ...
    + k*besselj(L,q*rho)*(W_5(1)*besselj(L+1,k*rho)-bessely(L+1,k*rho)*W_7(1)));

E1 = 1/4*pi^2*rho*(-(q*W_0(rho)*besselj(L+1,q)-k*besselj(L,q)*W_1(rho))*q*besselj(L+1,q*rho) ...
    + k*besselj(L,q*rho)*(W_4(1)*besselj(L+1,k*rho)+bessely(L+1,k*rho)*W_6(1)));


% Calculating the normalization constant
T_11 = @(r) integ1(q, besselj(L-1,r*q), besselj(L,r*q),besselj(L+1,r*q),besselj(L-1,r*q),besselj(L,r*q),besselj(L+1, r*q),r);
S_21 = @(r) integ1(k,besselj(L-1, r*k), besselj(L, r*k),besselj(L+1, r*k), besselj(L-1, r*k), besselj(L, r*k), besselj(L+1, r*k), r);
S_22 = @(r) integ1(k, besselj(L-1, r*k),besselj(L,r*k), besselj(L+1,r*k), bessely(L-1, r*k), bessely(L, r*k), bessely(L+1, r*k), r);
S_23 = @(r) integ1(k, bessely(L-1, r*k),bessely(L,r*k), bessely(L+1,r*k), bessely(L-1, r*k), bessely(L, r*k), bessely(L+1, r*k), r);
S_31 = @(r) integ1(q,besselj(L-1, r*q), besselj(L, r*q),besselj(L+1, r*q), besselj(L-1, r*q), besselj(L, r*q), besselj(L+1, r*q), r);
S_32 = @(r) integ1(q, besselj(L-1, r*q),besselj(L,r*q), besselj(L+1,r*q), bessely(L-1, r*q), bessely(L, r*q), bessely(L+1, r*q), r);
S_33 = @(r) integ1(q, bessely(L-1, r*q),bessely(L,r*q), bessely(L+1,r*q), bessely(L-1, r*q), bessely(L, r*q), bessely(L+1, r*q), r);

g_1 = @(r) T_11(r);
g_2 =  @(r) B1^2*S_21(r)+2*B1*C1*S_22(r) + C1^2*S_23(r);
g_3 = @(r) D1^2*S_31(r)+2*D1*E1*S_32(r)+E1^2*S_33(r);

N = (g_3(r_max)-g_3(1)+g_2(1)-g_2(rho)+g_1(rho)-g_1(0))^(-1/2);

B = N*B1;
C = N*C1;
D = N*D1;
E = N*E1;
