function sim = DDM_Kernel_Simulation(p)
% function sim = DDM_Kernel_Simulation(p)
%   run a simulation of psychophysical reverse correlation under
%   the assumption of drift diffusion model.
%
%   Input:
%       p -> a structure defining model parameters.
%            see below for the available parameters and their
%            default values.
%   Output:
%       sim -> a structure that contains simulation output.
%       sim.kernel.stim{i}: kernel aligned to stimulus onset when the model choice is i (1 or 2)
%       sim.kernel.resp{i}: (RT task only) kernel aligned to response when the model choice is i (1 or 2).
%       sim.bound_crossed_percent: percent the model crossed the bound
%                       within the simulated time (should be close to 100% to obtain an accurate result).
%       sim.cut_off_RT: the duration of kernel (default: cut off at median reaction time).
%       sim.param : parameters used to run the simulation (i.e. input p).
%
%
% Copyright (2018), Kiani Lab, NYU


%% default parameters
def.iters = 10000;   % number of simulated trials
def.t_max = 5000;    % number of simulated time steps in each trial (ms)
def.dt = 1;         % time unit, in ms

def.coh = 0;    % stimulus strength
def.k = 1;      % sensitivity parameter (drift rate = coh * k)
def.k0 = 0;         % base drift of the accumulators, equivalent of **bias**
def.B = [-30 30];   % lower and upper bounds (1 x 2) or (t_max x 2)
def.B0 = 0;         % base line of the accumulators (correspond to bias)
def.B0_SD = 0;      % fluctuation of baseline of the accumulators
def.sigma = 1;      % standard deviation of fluctuation of stimuli
def.stim_noise = 0; % noise added to the stimulus fluctuation (sensory noise)
def.dec_noise = 0;  % noise added to decision variable (decision noise)
def.w = 1;          % weight (1 x 1) or (t_max x 1)
def.non_dec_time = 0;    % average non-decision time
def.non_dec_time_sd = 0; % SD of non-decision time
def.termination_rule = {'RT', NaN}; % rule RT task:{'RT',NaN}, fixed:{'Fixed', stim_duration}

def.seed = NaN; % If specified, use this value to initialize matlab random
def.cut_off_RT = nan; % explicitly determine cut off RT (otherwise, median RT will be the cut off time)

%% setup
if nargin < 1
    p = def;
else
    fs = fieldnames(def);
    for n=1:length(fs)
        if ~isfield(p, fs{n})
            p.(fs{n}) = def.(fs{n});
        end
    end
end

sttime = tic;

    %fixed bound height over time. change it to simulate collapsing bounds
if size(p.B,1) == 1
    p.bound_height = repmat(p.B,[p.t_max,1]);
else
    p.bound_height = p.B;
end

    %fixed weights over time. change it to simulate dynamic weighting of sensory evidence 
if length(p.w) == 1
    p.weight = ones(1, p.t_max) * p.w;
else
    p.weight = p.w;
end

sim = struct();

if ~isnan(p.seed)
    %fix the random seed if you want to always obtain the identical result
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',p.seed));
end

%% run simulation
    % create stimuli
k = p.coh * p.dt;
s = p.sigma * sqrt(p.dt);

  % stimulus fluctuation
  % Sf .. stimulus fluctuation subtracted mean coherence level.
  % E .. sensory evidence used in DDM.
  % for the kernel analysis, Sf will be used.
Sf = normrnd(k, s, [p.t_max, p.iters]); % Sb = base
E = Sf' * p.k + p.k0;
Sf = (Sf - k)';

    %apply temporal weights
V = (E + randn(size(E)) * p.stim_noise) .* repmat(p.weight, [p.iters 1]);
    %add the initial value
V(:,1) = V(:,1) + p.B0 + randn(size(V,1),1) * p.B0_SD;

    %add decision noise
V = V + randn(size(V)) * p.dec_noise;


    %calculate the cumulative evidence
V = cumsum(V, 2);

    %find the choice and bound crossing time for each trial 
choice = nan(p.iters, 1);
bound_crossing = zeros(p.iters, 1);
RT = nan(p.iters, 1); % reaction time including non-decision time
dec_time = nan(p.iters, 1); % the exact timing of decision (bound crossing)

    %find bound crossing
for tr=1:p.iters
    t = find(V(tr,:)>=p.bound_height(:,2)' | V(tr,:)<=p.bound_height(:,1)', 1);
    if ~isempty(t)
        bound_crossing(tr) = 1;
        dec_time(tr) = t;
        if V(tr,t) >= p.bound_height(t,2)
            choice(tr) = 1;
        else
            choice(tr) = 2;
        end
    end
end

if strcmp(p.termination_rule{1}, 'Fixed') % in case of fixed duration task
    nI = find(bound_crossing==0);
    L = V(nI, p.termination_rule{2}) >= 0;
    choice(nI(L)) = 1;
    choice(nI(~L)) = 2;
    dec_time(nI) = p.termination_rule{2};
end

    % determine RT
if strcmp(p.termination_rule{1}, 'RT') % if RT task, RT is decision time + non decision time
    RT = dec_time + p.non_dec_time + round(randn(size(RT)) * p.non_dec_time_sd);
    RT(RT<=0) = 0; % clip at 0

    idx = RT > p.t_max;
    RT(idx) = NaN; % remove trials whose decision time is larger than stimulus duration or smaller than 0.
    bound_crossing(idx) = 0;
elseif strcmp(p.termination_rule{1}, 'Fixed') % if fixed duration task, RT is stimulus duration + non decision time
    RT = p.termination_rule{2} + p.non_dec_time + round(randn(size(RT)) * p.non_dec_time_sd);
    RT(RT<=0) = 0; % clip at 0
end


    %remove sensory fluctuation after bound crossing
for tr = find(bound_crossing)'
    t = RT(tr)+1;
    if t >= 1 && t <= size(Sf,2)
        Sf(tr, t:end) = NaN;
    end
end

%% kernel analysis
    %save the results
sim.choice = nansum(choice==2)/sum(~isnan(choice));
sim.medianRT = nanmedian(RT);
sim.sdRT = nanstd(RT);

if isnan(p.cut_off_RT)
    cut_off_RT = floor(nanmedian(RT));
    if isnan(cut_off_RT) || cut_off_RT > p.t_max 
        cut_off_RT = p.t_max;
    end
else
    cut_off_RT = p.cut_off_RT;
end


    %psychophysical kernel aligned to stimulus onset
sim.kernel.stim{1} = nanmean(Sf(choice==1,1:cut_off_RT),1);
sim.kernel.stim{2} = nanmean(Sf(choice==2,1:cut_off_RT),1);


if strcmp(p.termination_rule{1}, 'RT') % if RT task
        %average psychophysical kernel aligned to choice 
    S_choice = nan(p.iters,cut_off_RT);
    for tr = 1 : p.iters
        if bound_crossing(tr)==1
            st = max(0,RT(tr) - cut_off_RT)+1;
            en = min(p.t_max, RT(tr));
            valid_portion = st:en;
            S_choice(tr,end-length(valid_portion)+1:end) = Sf(tr,valid_portion);
        end
    end
    sim.kernel.resp{1} = nanmean(S_choice(choice==1,:),1);
    sim.kernel.resp{2} = nanmean(S_choice(choice==2,:),1);
else
    sim.kernel.resp{1} = [];
    sim.kernel.resp{2} = [];
end

sim.bound_crossed_percent = sum(bound_crossing)/length(bound_crossing) * 100;
sim.cut_off_RT = cut_off_RT;
sim.param = p;

fprintf('\tpercent cross bound: %1.1f%%, median RT = %1.3f (time spent %1.1fsec)\n', sim.bound_crossed_percent, nanmedian(RT), toc(sttime));
if sim.bound_crossed_percent < 99
    warning('Less than 99% of trials crossed the bound. The result could be inaccurate.');
end



