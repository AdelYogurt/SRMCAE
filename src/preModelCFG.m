function user_model=preModelCFG(cfg_filename)
% read model definition from cfg file
%
if nargin < 1
    error('readModelCFG: need input cfg file');
end

if length(cfg_filename) > 4
    if ~strcmpi(cfg_filename((end-3):end),'.cfg')
        cfg_filename=[cfg_filename,'.cfg'];
    end
else
    cfg_filename=[cfg_filename,'.cfg'];
end

% check input file
if exist(cfg_filename,'file') ~= 2
    error('readModelCFG: cfg file do not exist')
end

% initialize model
input_model=struct();

cfg_file=fopen(cfg_filename,'r');

% read cfg define line by line
while (~feof(cfg_file))
    string_read=regexprep(fgetl(cfg_file),'\s',''); % read char list and deblank
    if ~isempty(string_read) && string_read(1) ~= '%'
        string_list=strsplit(string_read,{'=','{','(',')','}',';','''','|',','});
        parameter=string_list{1};
        value=string_list(2:end);
        if isempty(value{end})
            value=value(1:end-1);
        end

        if isempty(value)
            error('readModelCFG: definition lack value');
        end

        % add parameter, if is number, convert char to number
        digital_value=zeros(1,length(value));
        digital_flag=1;
        for value_index=1:length(value)
            digital_value(value_index)=str2double(value{value_index});
            if isnan(digital_value(value_index))
                digital_flag=0;
                break;
            end
        end

        % if is string
        if digital_flag
            input_model.(parameter)=digital_value;
        else % if is number
            if length(value) == 1
                input_model.(parameter)=value{1};
            else
                input_model.(parameter)=value;
            end
        end
    end
end

fclose(cfg_file);
clear('cfg_file');

% process input value
user_model=UserModel();
user_model.input_cfg_filename=cfg_filename;
user_model.input_model=input_model;

% chamber
user_model.LC=readCFGParameter(input_model,'LEN_CHAMBER');
user_model.DC=readCFGParameter(input_model,'DIA_CHAMBER');
user_model.Dt=readCFGParameter(input_model,'DIA_NOZZLE_THROAT');
user_model.De=readCFGParameter(input_model,'DIA_NOZZLE_OUTLET');
if isfield(input_model,'NOZZLE_EFFICIENCY')
    user_model.NE=input_model.NOZZLE_EFFICIENCY;
else
    user_model.NE=1;
end

% propellant gas
user_model.GT=readCFGParameter(input_model,'GAS_TEMPERATURE');
user_model.GMW=readCFGParameter(input_model,'GAS_MOLECULAR');
user_model.GK=readCFGParameter(input_model,'GAS_SPECIFIC_RATIO');
user_model.PD=readCFGParameter(input_model,'PROPELLANT_DENSITY');

% propellant burn rate data
burn_rate_data=readCFGParameter(input_model,'PROPELLANT_BURN_RATE');
if mod(length(burn_rate_data),3) == 0
    user_model.burn_rate_data=reshape(burn_rate_data,3,length(burn_rate_data)/3)';
else
    error('preModelCFG: PROPELLANT_BURN_RATE data format incorrect');
end

% load grain data
switch input_model.GRAIN_FORMAT
    case 'xlsx'
        burn_area_data=readtable(input_model.GRAIN_FILE);
        burn_area_data=table2array(burn_area_data);
    case 'csv'
        burn_area_data=importdata(input_model.GRAIN_FILE);
    case 'txt'
        burn_area_data=importdata(input_model.GRAIN_FILE);
end
user_model.burn_area_data=checkInputGrain(burn_area_data);

% calculate geometry volume
user_model.grain_volume_initial=trapz(burn_area_data(:,1),burn_area_data(:,2));

end

function value=readCFGParameter(input_model,parameter)
if isfield(input_model,parameter)
    value=input_model.(parameter);
else
    error(['preModelCFG: cfg parameter lack ',parameter]);
end
end

function burn_area_data=checkInputGrain(burn_area_data)
% check burn area data input from grain file if corespondent to program
% format, if not, throw out error
%
if isempty(burn_area_data)
    error('inputCheck: burn_area_data is empty');
end
[rank,colume]=size(burn_area_data);
if colume~=2
    error('inputCheck: burn_area_data must be n x 2 matrix');
end

% first colume of first rank must be zero
if burn_area_data(1,1) ~= 0
    burn_area_data=[0,burn_area_data(1,2);burn_area_data];
end
end