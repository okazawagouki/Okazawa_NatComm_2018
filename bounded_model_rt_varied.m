function bounded_model_rt_varied(ID)
B = 10:5:50;

addpath ../functions/;


ID = str2double(ID);
B = B(ID);

folder = 'bounded_model_rt_varied/';
mkdir(folder);


% p.iters = 1000;
p.iters = 10000 * 100;
p.split_trials = 10000;
p.t_max = 10000;

p.t_frame = 1;
p.non_dec_time = 0;
p.subtract_time = 0;
p.non_dec_time_sd = 0;
p.error_no_reach = false;

p.cut_off_decision = false; % make decision when reach the time limit.
p.cut_off_RT = p.t_max; % fix at 1s

p.B = [-B, B]; % bounded
sim = DDM_Kernel_Simulation(p);

save(sprintf('%s/sim%02d.mat', folder, ID), 'sim', 'B');







