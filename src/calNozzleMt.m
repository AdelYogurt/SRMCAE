function mt=calNozzleMt(Pc,P3,Pa,At,Ae,NLC,K,C_K_0,C_K_1,C_K_2,R,T)
% function of calculate the one-dimensional isentropic flow
% calculate mass ejection rate of gas
%
% input: Pc is the pressure of combustion chamber. P3, Pa, At, Ae, NLC, K,
% C_K_0, C_K_1, C_K_2, R, T, are all the calculation parameters of frozen
% flow When the pressure of the combustion chamber is less than P3, the
% throat of the nozzle is subsonic flow, and the rest of the throat is
% supersonic. When the pressure is less than atmospheric pressure, no gas
% is ejected.
%
if (P3 < Pc)
    mt=NLC*C_K_0*Pc*At/sqrt(R*T); % mt
elseif (Pc <= Pa)
    mt=0;
else
    mt=NLC*Pc*Ae*sqrt(2*K/(K-1)*...
        ((Pa/Pc)^C_K_1-(Pa/Pc)^C_K_2))/sqrt(R*T); % mt
end
end