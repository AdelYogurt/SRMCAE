function dy=equationInteriorBallistics...
    (t,y,...
    V10,Dt,De,ARt,ARe,NE,T,R,K,C_K_0,C_K_1,C_K_2,P3,...
    Rou,omega_V,burn_rate_data,burn_area_data,Pa)
% Interior ballistic differential equations
%
% y(1) is the combustion chamber pressure
% y(2) is the total mass of fuel that has been ejected
% y(3:i) is the total mass of fuel that has been burned
% y(i:j) is the thickness of meat that has been burned
%
% dy(1) is the pressure increment of combustion chamber per unit time
% dy(2) is the fuel ejection mass per unit time
% dy(3:i) is the fuel combustion mass per unit time
% dy(i:j) is the combustion distance per unit time
%
Pc=y(1);me=y(2);mb=y(3);e=y(4);

% interpolation calculation of burning surface area
sigma=interpolationBrunArea(burn_area_data,e);
V1=V10-omega_V+mb/Rou;
dy(4)=interpolationBrunRate(Pc,burn_rate_data); % e£¬burning distance per unit time
dy(3)=sigma*dy(4)*Rou; % mb, propellant burned mass per unit time

At=(Dt+ARt*t*2).^2*pi/4;
Ae=(De+ARe*t*2).^2*pi/4;
dy(2)=calNozzleMt(Pc,P3,Pa,At,Ae,NE,K,C_K_0,C_K_1,C_K_2,R,T); % mt, gas ejection mass per unit time
dy(1)=R*T*(dy(3)-dy(2))/V1;
dy=[dy(1);dy(2);dy(3);dy(4)];
end