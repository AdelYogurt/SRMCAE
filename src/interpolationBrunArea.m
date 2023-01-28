function burn_area=interpolationBrunArea(burn_area_data,e)
% linear interpolation for calculate burn area
%
% input:
% burn_area_data (lenght x 2 matrix), e (thickness of burning grain)
%
e_index=size(burn_area_data,1); % 从最后一行开始搜索，查找比e小所在的行
while ((e_index > 1) && (burn_area_data(e_index,1) > e))
    e_index=e_index-1;
end
if ((e_index == size(burn_area_data,1)) || ...
        (e_index == 0))
    burn_area=0;
else
    % linear interpolation
    burn_area=burn_area_data(e_index,2)+...
        (burn_area_data(e_index+1,2)-burn_area_data(e_index,2))*...
        (e-burn_area_data(e_index,1))/...
        (burn_area_data(e_index+1,1)-burn_area_data(e_index,1));
end
end