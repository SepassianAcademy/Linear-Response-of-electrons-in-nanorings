function mateu = matrixelementunbound(k_01,k_1m, q_01,q_1m,r_max,rho,N0,B0,C0,D0,N1,B1,C1,D1,E1)
% k_01 =2.0;
% k_1m = 30.0;
% V_0 = 8;
% q_01 = sqrt(2*V_0-k_01^2);
% q_1m = sqrt(k_1m^2-2*V_0);
% r_max = 10;
% rho = 0.9;
% N0 = 1;
% B0 = 2;
% C0 = 3;
% D0 = 4;
% N1 = 5;
% B1 = 6;
% C1 = 7;
% D1 = 8;
% E1 = 9;




integ = @(a, b, L, Z_L, Z_L1, Zp_L, Zp_L1, r) r^2*(b*Z_L*Zp_L+a*Z_L1*Zp_L1)/(a^2-b^2)+2*r*(((L+1)*b^2-L*a^2)*Z_L*Zp_L1-b*a*Z_L1*Zp_L)/(a^2-b^2)^2;


integKJ = @(a,b,x)-(a*(x*(a^2+b^2)*besselj(1,b*x)+2*b*besselj(0,b*x))*besselk(1, a*x)+ ...
          besselk(0,a*x)*(-2*besselj(1,b*x)*b+x*besselj(0,b*x)*(a^2+b^2))*b)*x/(a^2+b^2)^2;

integKY = @(a,b,x)-(a*(x*(a^2+b^2)*bessely(1,b*x)+2*b*bessely(0,b*x))*besselk(1, a*x)+ ...
          besselk(0,a*x)*(-2*bessely(1,b*x)*b+x*bessely(0,b*x)*(a^2+b^2))*b)*x/(a^2+b^2)^2;

f_1 = @(r) N0*N1*integ(1i*q_01, q_1m, 0, besselj(0,1i*q_01*r),besselj(1, 1i*q_01*r), besselj(0, q_1m*r), besselj(1,q_1m*r), r);
f_21 = @(r) B0*B1*integ(k_01, k_1m, 0, besselj(0,k_01*r),besselj(1,k_01*r),besselj(0, k_1m*r), besselj(1, k_1m*r), r);
f_22 = @(r) B0*C1*integ(k_01, k_1m, 0, besselj(0, k_01*r),besselj(1, k_01*r), bessely(0, k_1m*r), bessely(1, k_1m*r), r);
f_23 = @(r) C0*B1*integ(k_01, k_1m, 0, bessely(0, k_01*r),bessely(1, k_01*r), besselj(0, k_1m*r), besselj(1, k_1m*r), r);
f_24 = @(r) C0*C1*integ(k_01, k_1m, 0, bessely(0, k_01*r), bessely(1, k_01*r), bessely(0, k_1m*r), bessely(1, k_1m*r), r);

f_2 = @(r) f_21(r)+f_22(r)+f_23(r)+f_24(r);




f_3 = @(r) D0*D1*integKJ(q_01,q_1m,r)+D0*E1*integKY(q_01,q_1m,r);
mateu = 1/2*(f_3(r_max)-f_3(1)+f_2(1)-f_2(rho)+f_1(rho)-f_1(0));

