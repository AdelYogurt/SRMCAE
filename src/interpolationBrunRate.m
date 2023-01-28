function burning_rate=interpolationBrunRate(Pc,burning_rate_data)
% The burning rate calculation function uses the Vienna formula
%
% input:
% Pc: combustion chamber pressure, burning_rate_data: burning rate data
%
% burning rate data: the first is the burning rate coefficient, the second is the
% pressure exponent, and the third is the maximum pressure range.
% 
Pc_narrow=Pc*1e-6;
index=1;
while ((index < size(burning_rate_data,1)) && (Pc_narrow > burning_rate_data(index,3)))
    index=index+1;
end
if Pc < 0
    burning_rate=0;
else
    burning_rate=burning_rate_data(index,1)*Pc_narrow^burning_rate_data(index,2)+0;
end
end