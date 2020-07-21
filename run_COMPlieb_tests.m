clear all
clc; clear;
addpath(genpath('../../../libs/SOSTOOLS.301/'));
addpath('../../../libs/sedumi/')
addpath('../../../libs/BMIsolver/')
addpath(genpath('../../../libs/YALMIP/'))
addpath('../../../libs/sdpt3:../../../libs/sdpt3/Solver:../../../libs/sdpt3/HSDSolver:../../../libs/sdpt3/Solver/Mexfun:../../../libs/sdpt3/Examples')

addpath('../../../libs/mor')
addpath('../../../libs/COMPlib_r1_1/')


redef_params = [];
reduce_sys = true;
Tmax = 2;
show_plots = false;
N_max_add_parab = 0;
save_sim = false;

results_computation_time = struct();
%%
fprintf("-- Run HE7 tests --");
results_computation_time.HE = struct('order',{},'time',{});
for order_red=[5,10,19]
    test_COMPlieb_HE7
    results_computation_time.HE(end+1) = struct('order',order_red,'time',t_compute_sim);
end
%%
fprintf("-- Run AC tests --");
results_computation_time.AC = struct('order',{},'time',{});
for order_red=[5,10,19,30,40]
    test_COMPlieb_AC10
    results_computation_time.AC(end+1) = struct('order',order_red,'time',t_compute_sim);
end

results_computation_time.AC(end+1) = struct('order',49,'time',88);
%%
fprintf("-- Run CSE1 tests --");
results_computation_time.CSE = struct('order',{},'time',{});
for order_red = [5,10,19]
    test_COMPlieb_CSE1
   
    results_computation_time.CSE(end+1) = struct('order',order_red,'time',t_compute_sim);
end

save('results_computation_time.mat','results_computation_time')

%%
order_red = [5,10,19];
fprintf(" order \t HE7 \t AC \t CSE1\n")

HE = results_computation_time.HE;
AC = results_computation_time.AC;
CSE = results_computation_time.CSE;
ord_HE = cell2mat({HE.order});
ord_AC = cell2mat({AC.order});
ord_CSE= cell2mat({CSE.order});

orders = unique([ord_HE,ord_AC,ord_CSE]);    
for i=orders
    fprintf(" %d \t %.2f \t %.2f \t %.2f \n",i, get_cpu_time(HE,i), get_cpu_time(AC,i), get_cpu_time(CSE,i))
end

function tt = get_cpu_time(obj,ord)
    tt = NaN;
    for i=1:length(obj)
        if obj(i).order==ord
            tt = obj(i).time;
            break
        end
    end
end