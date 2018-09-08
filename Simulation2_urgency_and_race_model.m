
%% Simulation2_urgency_and_race_model
%
% This script shows psychophysical kernel under DDM with urgency and under race model
% (corresponds to Fig. 6 in the paper).
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

% DDM with urgency
p.t_max = 10000;
B = 80;
urg_half = 700;
urg_asymp = 80;
b = (0:p.t_max-1)'*urg_asymp./((0:p.t_max-1)'+urg_half);
p.B = [-B+b(:), B-b(:)]; % Bound
p.non_dec_time = 0;
p.cut_off_RT = 1000;

fprintf('Running DDM with urgency simulation...\n');
urg_sim = DDM_Kernel_Simulation(p);

% Race: input correlation
clear p;
p.t_max = 10000;
p.B = [-30 30];
p.rB = [Inf -Inf];
p.interaction = 0;
p.rho = -0.2;
p.stim_noise = 1;

p.non_dec_time = 0;
p.cut_off_RT = 1000;

fprintf('Running Race model with input correlation 0.2...\n');
input_sim = RACE_Kernel_Simulation(p);


% Race: reflective bound
clear p;
p.t_max = 10000;
p.B = [-30 30];
p.rB = [10, -10];
p.rho = -1;
p.interaction = 0;
p.non_dec_time = 0;
p.cut_off_RT = 1000;

fprintf('Running Race model with reflective bound -10...\n');
reflect_sim = RACE_Kernel_Simulation(p);


% Race: leak and mutual inhibition
clear p;
p.t_max = 10000;
p.stim_noise = 0;
p.B = [-60 60];
p.B0 = 30;
p.rB = [0 -0];
p.rho = -1;
p.interaction = 0.003;
p.decay = 0.003;
p.booster = p.B0 * (p.interaction + p.decay);
p.non_dec_time = 0;
p.cut_off_RT = 1000;

fprintf('Running Race model with leak/mutual inhibition = 1...\n');
lm_sim = RACE_Kernel_Simulation(p);


%% Show figure
figure('color', 'w', 'position', [100 100 400 700]);

subplot(4,2,1);
show_kernel(urg_sim, 'stim', smoothing_wid, xrange, yrange);
title('Urgency');

subplot(4,2,2);
show_kernel(urg_sim, 'resp', smoothing_wid, xrange, yrange);


subplot(4,2,3);
show_kernel(input_sim, 'stim', smoothing_wid, xrange, yrange);
title('Race: Input correlation 0.2');

subplot(4,2,4);
show_kernel(input_sim, 'resp', smoothing_wid, xrange, yrange);


subplot(4,2,5);
show_kernel(reflect_sim, 'stim', smoothing_wid, xrange, yrange);
title('Race: Reflective bound -10');

subplot(4,2,6);
show_kernel(reflect_sim, 'resp', smoothing_wid, xrange, yrange);


subplot(4,2,7);
show_kernel(lm_sim, 'stim', smoothing_wid, xrange, yrange);
title('Race: Leak/Inhibition = 1');

subplot(4,2,8);
show_kernel(lm_sim, 'resp', smoothing_wid, xrange, yrange);





