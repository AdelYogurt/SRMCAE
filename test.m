clc;
clear;
close all hidden;

user_data=UserData();
user_data=dataInputProcess(user_data);

omega_V = user_data.grain_volume_initial;
Pa = user_data.perssure_atmosphere;
T = user_data.gas_temperature;
M = user_data.gas_molecular_weight;
K = user_data.gas_specific_heat_ratio;
Rou = user_data.powder_density;
V10 = user_data.combustion_chamber_volume;
Dt = user_data.Dt;
De = user_data.De;
ARt = user_data.ARt;
ARe = user_data.ARe;
NLC = user_data.nozzle_loss_coefficient;

R = 8314.472 / M;
At = Dt * Dt * pi / 4;% ÔÝÊ±²»¿¼ÂÇÉÕÊ´
Ae = De * De * pi / 4;% ÔÝÊ±²»¿¼ÂÇÉÕÊ´
K_plus = K + 1;
K_sub = K - 1;
C_K_0 = sqrt(K) * (2 / (K_plus)) ^ ((K_plus) / 2 /(K_sub));
C_K_1 = 2 / K;
C_K_2 = (K_plus) / K;
KFC = sqrt(K / R * ((2 / K_plus) ^ (K_plus / K_sub)));
Aac = sqrt(2 * K * R * T / K_plus);

lambda_e=fzero(@(x) x*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*x*x)^(1/K_sub)-At/Ae,2);
P1_P0=(1-K_sub*lambda_e*lambda_e/K_plus)^(K/K_sub);
P1=Pa/P1_P0;
Mae=sqrt(2*lambda_e*lambda_e/K_plus/(1-K_sub*lambda_e*lambda_e/K_plus));
P2_P0=P1_P0*(2*K*Mae*Mae/K_plus-K_sub/K_plus);
P2=Pa/P2_P0;
lambda_e_sub=fzero(@(x) x*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*x*x)^(1/K_sub)-At/Ae,0.5);
P3_P0=(1-K_sub*lambda_e_sub*lambda_e_sub/K_plus)^(K/K_sub);
P3=Pa/P3_P0;

Pc_list=(linspace(0.10,10,11))*1e6;
Sb=21000*1e-6;
mt_list=zeros(1,11);
mb_list=zeros(1,11);
Pc=fzero(@(Pc) Sb*burnRateInterpolation(Pc,user_data.burning_rate_data)*Rou-nozzleGetMt(Pc,P3,Pa,At,Ae,NLC,K,C_K_0,C_K_1,C_K_2,R,T),10e6);
for Pc_index=1:size(Pc_list,2)
    Pc=Pc_list(Pc_index);
    mt_list(Pc_index)=nozzleGetMt(Pc,P3,Pa,At,Ae,NLC,K,C_K_0,C_K_1,C_K_2,R,T);
    mb_list(Pc_index)=Sb*burnRateInterpolation(Pc,user_data.burning_rate_data)*Rou;
end
plot(Pc_list,mt_list,Pc_list,mb_list);


