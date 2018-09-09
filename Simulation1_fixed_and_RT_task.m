
%% Simulation1_fixed_and_RT_task
%
% This script shows psychophysical kernel under fixed duration and reaction
% time task with and without non-decision time (corresponds to Fig. 3 in
% the paper).
%
% The original simulation is computationally intensive (10^6 trials). As a
% short cut, this script runs a smaller number of trials (10^4) but applys
% a boxcar smoothing (50 ms) to reduce the noise in kernel. For the actual
% production of the result, no smoothing should be applied.
% 
%
%
% Copyright, Kiani Lab, NYU

%% Initialize
clear;

smoothing_wid = 50; %ms
xrange = [0 1200];
yrange = 0:1:3;



%% Run simulation

% Fixed duration task
p.termination_rule = {'Fixed', 1000}; % 1000ms
p.B = [-30 30]; % Bound
p.non_dec_time = 0;
p.error_no_reach = false;

fprintf('Running fixed duration task simulation...\n');
fix_sim = DDM_Kernel_Simulation(p);

% RT task without non-decision time
clear p;
p.termination_rule = {'RT', NaN};
p.t_max = 5000;
p.B = [-30 30]; % Bound
p.non_dec_time = 0;
p.cut_off_RT = 1000;

fprintf('Running RT task without non-decision time simulation...\n');
RT_sim = DDM_Kernel_Simulation(p);

% RT task with non-decision time
clear p;
p.termination_rule = {'RT', NaN};
p.t_max = 5000;
p.B = [-30 30]; % Bound
p.non_dec_time = 300;
p.non_dec_time_sd = 100;
p.cut_off_RT = 1000;

fprintf('Running RT task with non-decision time simulation...\n');
RT_ndec_sim = DDM_Kernel_Simulation(p);



%% Show figure
figure('color', 'w', 'position', [100 100 400 600]);

subplot(3,2,1);
show_kernel(fix_sim, 'stim', smoothing_wid, xrange, yrange);
title('Fixed duration task');

subplot(3,2,3);
show_kernel(RT_sim, 'stim', smoothing_wid, xrange, yrange);
title('RT task without non-decision time');

subplot(3,2,4);
show_kernel(RT_sim, 'resp', smoothing_wid, xrange, yrange);

subplot(3,2,5);
show_kernel(RT_ndec_sim, 'stim', smoothing_wid, xrange, yrange);
title('RT task with non-decision time');

subplot(3,2,6);
show_kernel(RT_ndec_sim, 'resp', smoothing_wid, xrange, yrange);







