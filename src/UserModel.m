classdef UserModel < handle
    properties
        % input data from cfg file
        input_cfg_filename;
        input_model;

        % environmenet parameter
        perssure_atmosphere = 101325;
        
        % combustion chamber parameter
        LC = 0; % characteristic length of chamber
        DC = 0; % characteristic inner diameter of chamber
        Dt = 0; % nozzle throat diameter
        De = 0; % nozzle outlet diameter
        NE = 0.85; % nozzle efficiency
        chamber_volume = 0; % volume of combustion chamber

        % corrosion rate
        ARt = 0;
        ARe = 0;
        
        % grain geometry data
        burn_area_data;
        grain_volume_initial = 0;

        % propellant gas data
        GT = 0; % average temperature of propellant gas
        GMW = 0; % propellant gas average molecular weight
        GK = 0; % specific heat ratio of propellant gas
        PD = 0; % propellant density
        
        % propellant burn rate data
        burn_rate_data;
        
        % result
        t_list; % ʱ��
        Pc_list;% ȼ����ѹǿ
        mb_total_list;% ��ȼ��������
        mt_total_list;% ��ȼ�������
        e_list;% ���
        BA_list;%ȼ�����
        mb_list;% ȼ����������
        mt_list;% ȼ���������
        Ve_list;% ȼ�������ٶ�
        F_list;% ����
        As_At_list;% �����������λ��
    end
    
end