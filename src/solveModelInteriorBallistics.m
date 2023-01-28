function solveModelInteriorBallistics(user_model)
% to solve the interior ballistics, calculate the combustion chamber
% pressure, total gas production, total gas ejection, meat thickness,
% combustion area, mass flow rate of gas generation, mass flow rate of gas
% ejection, velocity of nozzle airflow, thrust and time curve of normal
% shock wave position in nozzle. 
% All data are recorded in user_model after
% calculation.
%
if nargin < 1
    error('solveModelInternalBallistics: lack user_model input')
end

% get parameter from user_model
Pa = user_model.perssure_atmosphere;
% combustion chamber parameter
LC = user_model.LC;
DC = user_model.DC;
Dt = user_model.Dt;
De = user_model.De;
NE = user_model.NE;
user_model.chamber_volume=LC*DC*DC/4*pi;
V10 = user_model.chamber_volume;
% corrosion rate
ARt = user_model.ARt;
ARe = user_model.ARe;
% grain geometry data
burn_area_data = user_model.burn_area_data;
user_model.grain_volume_initial=trapz(burn_area_data(:,1),burn_area_data(:,2));
omega_V = user_model.grain_volume_initial;
% propellant gas data
T = user_model.GT;
M = user_model.GMW;
K = user_model.GK;
Rou = user_model.PD;
% propellant burn rate data
burn_rate_data = user_model.burn_rate_data;

if (V10 < omega_V)
   error('solveModelInternalBallistics: combustion chamber volume less than grain volume');
end

% to judge when to stop solve
e_max=max(burn_area_data(:,1));

% Pre-processing, calculate the parameters required to solve the equation
R = 8314.472 / M;
% erosive is not considered for the time being
At = Dt * Dt * pi / 4;
Ae = De * De * pi / 4;
K_plus = K + 1;
K_sub = K - 1;
C_K_0 = sqrt(K) * (2 / (K_plus)) ^ ((K_plus) / 2 /(K_sub));
C_K_1 = 2 / K;
C_K_2 = (K_plus) / K;
KFC = sqrt(K / R * ((2 / K_plus) ^ (K_plus / K_sub)));
Aac = sqrt(2 * K * R * T / K_plus);

% The velocity coefficient when the outlet is supersonic
lambda_e=fzero(@(x) x*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*x*x)^(1/K_sub)-At/Ae,2); 
P1_P0=(1-K_sub*lambda_e*lambda_e/K_plus)^(K/K_sub);
P1=Pa/P1_P0;
% The outlet Mach number of the normal shock wave at the nozzle exit
% increases to atmospheric pressure after passing through the normal shock
% wave.
Mae=sqrt(2*lambda_e*lambda_e/K_plus/(1-K_sub*lambda_e*lambda_e/K_plus));
P2_P0=P1_P0*(2*K*Mae*Mae/K_plus-K_sub/K_plus);
P2=Pa/P2_P0;
% The velocity coefficient when the outlet is subsonic
lambda_e_sub=fzero(@(x) x*(K_plus/2)^(1/K_sub)*(1-K_sub/K_plus*x*x)^(1/K_sub)-At/Ae,0.5);
P3_P0=(1-K_sub*lambda_e_sub*lambda_e_sub/K_plus)^(K/K_sub);
P3=Pa/P3_P0;

% The time curves of combustion chamber pressure, total gas production,
% total gas ejection and meat thickness are calculated by solving the
% interior ballistic equation.
solve_time=5;
solve_result=ode45(@(t,y) equationInteriorBallistics...
    (t,y,...
    V10,Dt,De,ARt,ARe,NE,T,R,K,C_K_0,C_K_1,C_K_2,P3,...
    Rou,omega_V,burn_rate_data,burn_area_data,Pa),[0,solve_time],[Pa;0;0;0]);
while ((max(solve_result.y(4,:)) < e_max) && (solve_time < 1000))
    % if motor no burn out, extend solve time
    solve_time=solve_time+5;
    solve_result=odextend(solve_result, @(t,y) equationInternalBallistics...
        (t,y,...
        V10,Dt,De,ARt,ARe,NE,T,R,K,C_K_0,C_K_1,C_K_2,P3,...
        Rou,omega_V,burn_rate_data,burn_area_data,Pa),solve_time);
end
t_list=solve_result.x;
Pc_list=solve_result.y(1,:);
mt1_total_list=solve_result.y(2,:);
mb_total=solve_result.y(3,:);
e_list=solve_result.y(4,:);

% locate when engine engine burn out
for perssure_out_index=1:length(t_list)
    if Pc_list(perssure_out_index) > Pa*1.05
        break;
    end
end
for perssure_out_index=perssure_out_index:length(t_list)
    if Pc_list(perssure_out_index) < Pa*1.05
        break;
    end
end
t_list(perssure_out_index:end)=[];
Pc_list(perssure_out_index:end)=[];
mb_total(perssure_out_index:end)=[];
mt1_total_list(perssure_out_index:end)=[];
e_list(perssure_out_index:end)=[];

% The time curves of combustion area, mass flow rate of gas generation,
% mass flow rate of gas velocity, thrust and normal shock wave position in
% nozzle are calculated by post-processing.
user_model.t_list = t_list;
user_model.Pc_list = Pc_list;
user_model.mb_total_list = mb_total;
user_model.mt_total_list = mt1_total_list;
user_model.e_list = e_list;
user_model.BA_list = zeros(1,length(t_list));
user_model.mb_list = zeros(1,length(t_list));
user_model.mt_list = zeros(1,length(t_list));
user_model.Ve_list = zeros(1,length(t_list));
user_model.F_list = zeros(1,length(t_list));
user_model.As_At_list = zeros(1,length(t_list));
for y_index = 1:length(t_list)
    % Interpolation calculation of burning surface area curve
    user_model.BA_list(y_index) = interpolationBrunArea(burn_area_data,user_model.e_list(y_index));
    % Calculate gas generation flow per unit time
    user_model.mb_list(y_index) = interpolationBrunRate(user_model.Pc_list(y_index),burn_rate_data);
    % Calculate the gas ejection flow per unit time
    user_model.mt_list(y_index) = calNozzleMt(user_model.Pc_list(y_index),P3,Pa,At,Ae,NE,K,C_K_0,C_K_1,C_K_2,R,T);
    % Calculate the outlet velocity, shock wave in the nozzle, outlet pressure
    [Ve,As_At,Pe]=calNozzleVP(Pc_list(y_index),lambda_e,Aac,Pa,P1_P0,P2,P3,K,K_sub,K_plus,At,Ae);
    user_model.Ve_list(y_index) = Ve;
    user_model.F_list(y_index) = user_model.Ve_list(y_index) * user_model.mt_list(y_index)+ (Pe - Pa) * Ae;
    user_model.As_At_list(y_index) =As_At;
end

end