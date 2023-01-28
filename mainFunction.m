clc;
clear;
close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\input']);

cfg_filename='36Motor.cfg';

% read data from cfg file
user_model=preModelCFG(cfg_filename);

C_BV=0.5; % ȼ��˥��ϵ��������ʵ���������KNSBΪ0.5
burn_data_number=size(user_model.burn_rate_data,1);
user_model.burn_rate_data=user_model.burn_rate_data.*[ones(burn_data_number,1)*C_BV,ones(burn_data_number,1),ones(burn_data_number,1)];

% user_model.burn_area_data=erosiveCombustion(user_model.burn_area_data,user_model.LC/user_model.DC); % ������ʴȼ�ն�ȼ��������ߵĸı�

grain_mass=user_model.grain_volume_initial*user_model.PD;
disp("grain_volume: "+user_model.grain_volume_initial);
disp("grain_mass: "+grain_mass);

solveModelInteriorBallistics(user_model); % �������

displayModel(user_model); % ���ݿ��ӻ�
