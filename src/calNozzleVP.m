function [Ve,As_At,Pe]=calNozzleVP(Pc,lambda_e,Aac,Pa,P1_P0,P2,P3,K,K_sub,K_plus,At,Ae)
% function to calculate one-dimensional isentropic flow of nozzle 
% The ejection velocity, nozzle normal shock wave position and outlet pressure are calculated. 
% Ve: nozzle exit velocity
% As_At: the normal shock wave inside the nozzle ( the combustion chamber pressure will only appear in P2 and P3 )
% pe: outlet airflow pressure 
% Pc, lambda_e, Aac, Pa, P1_P0, P2, P3, K, K_sub, K_plus, At, Ae are all calculation parameters of airflow freezing flow.
%
if (Pc <= Pa)
    lambda=0;
    As_At=1;
    Pe=Pa;
elseif (Pc <= P3)
    lambda=fzero(@(x) (1-K_sub*x*x/K_plus)^(K/K_sub)-Pa/Pc,0.5);
    if isnan(lambda)
        lambda=0;
    end
    As_At=1;
    Pe=Pa;
elseif ((P3 < Pc)&&(Pc <= P2))
    lambda=fzero(@(x) x*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*x*x)^(1/K_sub)/(1-K_sub*x*x/K_plus)^(K/K_sub)-Pc/Pa*At/Ae,0.5);
    sigma=At/Ae/(lambda*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*lambda*lambda)^(1/K_sub));
    Ma_s1=fzero(@(x)  (2*K*x*x/K_plus-K_sub/K_plus)^(-1/K_sub)*(K_plus*x*x/(K_sub*x*x+2))^(K/K_sub)-sigma,2);
    As_At=Ma_s1*(2/K_plus*(1+K_sub/2*Ma_s1*Ma_s1))^(-K_plus/2/K_sub);
    Pe=Pa;
else
    lambda=lambda_e;
    As_At=0;
    Pe=Pc*P1_P0;
end
Ve=lambda*Aac;
end