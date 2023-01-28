clc;
clear;
close all hidden;

user_data=UserData();
% 初始化计算
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
At = Dt * Dt * pi / 4;% 暂时不考虑烧蚀
Ae = De * De * pi / 4;% 暂时不考虑烧蚀
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

max_Pc=10;

Pc_cal_list=(Pa*1e-6:0.02:max_Pc)'*1e6;
thrust_cal_list=zeros(size(Pc_cal_list,1),1);

for Pc_i =1:size(Pc_cal_list,1)
    Pc=Pc_cal_list(Pc_i);
    me = nozzleGetMt(Pc,P3,Pa,At,Ae,NLC,K,C_K_0,C_K_1,C_K_2,R,T);
    [Ve,As_At,Pe]=nozzleGetVP(Pc,lambda_e,Aac,Pa,P1_P0,P2,P3,K,K_sub,K_plus,At,Ae);
    F = me * Ve + (Pe - Pa) * Ae;
    thrust_cal_list(Pc_i)=F;
end

load("thrust_data");
thrust_data=thrust_data_120_1;
time=thrust_data(:,1);

result_Pc=thrust_data(:,1);
for result_i =1:size(result_Pc,1)
    thrust=thrust_data(result_i,2);
    result_Pc(result_i,1)=findNearPc(thrust,thrust_cal_list,Pc_cal_list);
end

% 平衡压强法计算初始燃面
result_BA=thrust_data(:,1);
for result_i =1:size(result_BA,1)
    Pc=result_Pc(result_i,1);
    me = nozzleGetMt(Pc,P3,Pa,At,Ae,NLC,K,C_K_0,C_K_1,C_K_2,R,T);
    result_BA(result_i,1)=me/burnRateInterpolation(Pc,user_data.burning_rate_data)/Rou;
end

% plot(time*1e3,result_Pc*1e-6);
% xlabel('\fontsize{8}\bft  (ms)');
% ylabel('\fontsize{8}\bfPc  (MPa)');
% title('\fontsize{8}\bft-Pc曲线'); 

plot(time*1e3,result_BA*1e6);
xlabel('\fontsize{8}\bft  (ms)');
ylabel('\fontsize{8}\bfBA  (mm^2)');
title('\fontsize{8}\bft-BA曲线'); 

grid on;

function Pc=findNearPc(thrust,thrust_cal_list,Pc_cal_list)
error_min=1e6;
for Pc_index=1:size(Pc_cal_list,1)
    error=abs(thrust-thrust_cal_list(Pc_index,1));
    if (error<error_min)
        error_min=error;
        Pc=Pc_cal_list(Pc_index,1);
    end
end
end