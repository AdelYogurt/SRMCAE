function displayModel(user_model)
% Data visualization, curve drawing
% 
figure(1);
subplot(2,2,1);  
plot(user_model.t_list*1e3,user_model.Pc_list*1e-6,'linewidth',2); 
grid on;  
xlabel('\fontsize{8}\bft  (ms)'); 
ylabel('\fontsize{8}\bfPc  (MPa)'); 
title('\fontsize{8}\bft-Pc曲线'); 
subplot(2,2,2)  
plot(user_model.t_list*1e3,user_model.Ve_list,'linewidth',2); 
grid on; 
xlabel('\fontsize{8}\bft  (ms)');
ylabel('\fontsize{8}\bfve  (m/s)');
title('\fontsize{8}\bft-ve曲线'); 
subplot(2,2,3)  
plot(user_model.t_list*1e3,user_model.F_list,'linewidth',2); 
grid on;  
xlabel('\fontsize{8}\bft  (ms)'); 
ylabel('\fontsize{8}\bfF  (N)');
title('\fontsize{8}\bft-F曲线'); 
subplot(2,2,4) 
plot(user_model.t_list*1e3,user_model.BA_list*1e6,'linewidth',2);
grid on; 
xlabel('\fontsize{8}\bft  (ms)');
ylabel('\fontsize{8}\bfBA  (mm^2)');
title('\fontsize{8}\bft-BA曲线')


tspan = length(user_model.t_list)/20; 
tspan = 1:ceil(tspan):length(user_model.t_list);
tspan(end) = length(user_model.t_list);  
disp('        t(ms)     Pc(MPa)'); 
format short g;  
result=[(user_model.t_list(tspan)*1000)',(user_model.Pc_list(tspan)*1e-6)'];
disp(result)
end