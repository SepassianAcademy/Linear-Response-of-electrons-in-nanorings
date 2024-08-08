function [matep,matem] = matelement(L,rho,k_Ln,k_Lm1m, k_Lp1m)



matem = (2*k_Ln*k_Lm1m*(rho*besselj(L,k_Ln*rho)*besselj(L-1,k_Lm1m*rho)-besselj(L,k_Ln)*besselj(L-1,k_Lm1m)))/(rho*(k_Ln^2-k_Lm1m^2)^2*...
        sqrt(besselj(L,k_Ln*rho)^2-besselj(0,k_Ln)^2)*sqrt(besselj(L-1,k_Lm1m*rho)^2-besselj(L-1,k_Lm1m)^2));


matep = (2*k_Ln*k_Lp1m*(rho*besselj(L,k_Ln*rho)*besselj(L+1,k_Lp1m*rho)-besselj(L,k_Ln)*besselj(L+1,k_Lp1m)))/(rho*(k_Ln^2-k_Lp1m^2)^2*...
        sqrt(besselj(L,k_Ln*rho)^2-besselj(0,k_Ln)^2)*sqrt(besselj(L+1,k_Lp1m*rho)^2-besselj(L+1,k_Lp1m)^2));

