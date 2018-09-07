function sim = RACE_Kernel_Simulation(p)
% p is a struct to define parameter. See below for the default parameters.

%% default parameters
def.iters = 10000;   % number of simulated trials
def.t_max = 5000;    % number of simulated time steps in each trial
def.split_trials = 0; % split trials when the memory consumption is large (if larger than 0, this corresponds to number of trials to split).
def.parfor = false;

def.dt = 1;         % time unit, in ms
def.t_frame = 1;   % duration of one noise frame

def.termination_rule = {'RT', NaN}; % rule RT task:{'RT',NaN}, fixed:{'Fixed', stim_duration}

def.coh = 0;    % coherence (if not scalar, coh is randomly chosen for eath trial)
def.k = 1;      % sensitivity parameter (drift rate = coh * k)
def.k0 = 0;         % base drift of the accumulators, equivalent of **bias**
def.B = [-30 30];   % decision bound (1 x 2) or (t_max x 2)
def.rB = [Inf -Inf]; % reflection bound (1 x 2) or (t_max x 2)
def.B0 = 0;         % base line of the accumulators (correspond to bias)

def.rho = -1;   % correlation between the two decision variables (-1 = diffusion model)
def.interaction = 0; % 0= no interaction. 1=very strong interaction.
def.decay = 0; % 0= no decay. 1=very strong decay.
def.booster = 0; % constant input to balance the inhibition and excitation.

def.sigma = 1;      % standard deviation of fluctuation of stimuli
def.signal_compress = NaN;  % non-linearity in signal. put function handle to define comporessing function.
def.stim_noise = 0; % noise added to the stimulus fluctuation (sensory noise)
def.dec_noise = 0;  % noise added to decision variable (decision noise)
def.w = 1;          % weight (1 x 1) or (t_max x 1)
def.non_dec_time = 0;    % average non-decision time
def.non_dec_time_sd = 0; % SD of non-decision time
def.non_dec_time_dist = 'normal'; % shape of non-decision time distribution. 'normal' or 'log normal'.
def.subtract_time = 0;   % Time subtracted from RT during kernel analysis to compensate for non-decision time
def.include_dec_frame = 1;  % whether to include the coherence exactly at the time of crossing bound.
                            % 1.. include the time when the bound crossed.
                            % 0.. do not include the time when the bound crossed.
                            % >1 .. include the time after the bound cross (t-1 from the bound cross)
                            % <0 .. include only the time before the bound cross (t-1 from the bound cross)
    %what to do with decision variable and momentary evidence vector when bound crossing happens 
def.v_past_bound =  NaN; %can be 'bound' or NaN. choosing them depends on whether dealing with a RT or fix duration task. 
                        %Choose NaN if you want a trial that reaches the bound not to contribute to average traces after bound crossing. 
                        %Choose 'bound' if you want it to contribute but be fixed at the bound value
def.s_past_bound =  NaN; %can be 'none' or NaN. choosing them depends on whether dealing with a RT or fix duration task.  
                        %Choose NaN if you want a trial that reaches the bound not to contribute to average traces after bound crossing. 
                        %Choose 'none' if you want it to contribute even past bound crossing. 
def.cut_off_decision = false;   % if true, make a decision when the process reaches to t_max even when it does not reach the bounds.
def.get_raw_data = false;   %whether to get raw S, DV, and other computed values.
def.seed = NaN;
def.cut_off_RT = nan; % explicitly determine cut off RT
def.error_no_reach = true; % end the program with error when less than 95% of trials reach the bound.

%% setup
if nargin < 1
    p = def;
else
    p = safeStructAssign(def, p);
end
    % split trials to reduce memory consumption
if p.split_trials
    disp('split_trials is on. Trials will be split.');
    sim = RACE_Kernel_Simulation_split_trial(p);
    return;
end
sttime = tic;

    %fixed bound height over time. change it to simulate collapsing bounds
if size(p.B,1) == 1
    p.bound_height = repmat(p.B,[p.t_max,1]);
else
    p.bound_height = p.B;
    p = rmfield(p, 'B');
end
p.bound_height(:,1) = -p.bound_height(:,1); % flip bound sign (because this is a race model)
rB = p.rB;
rB(:,1) = -rB(:,1); % flip bound sign (because this is a race model)

    %fixed weights over time. change it to simulate dynamic weighting of sensory evidence 
if length(p.w) == 1
    p.weight = ones(1, p.t_max) * p.w;
else
    p.weight = p.w;
    p = rmfield(p, 'w');
end

sim = struct();

if ~isnan(p.seed)
    %fix the random seed to reduce variability in the model exploration pahse. 
    %For final simulations, you should use random seed
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',p.seed));
end

%% run simulation
    % create stimuli
if isscalar(p.coh)
    coh = p.coh;
else
    idx = ceil(rand(1, p.iters) * length(p.coh));
    p.coh = p.coh(:)';
    coh = ones(ceil(p.t_max/p.t_frame),1) * p.coh(idx);
end
k = coh * p.dt;
s = p.sigma * sqrt(p.dt);

  % stimulus fluctuation
  % Sb .. stimulus fluctuation subtracted mean coherence level. Temporal resolution = t_frame ms.
  % S .. real stimulus values (sensory evidence) used in DDM. Temporal resolution = 1ms.
  % for the kernel analysis, Sb will be used.
Sb = normrnd(k, s, [ceil(p.t_max/p.t_frame), p.iters]); % Sb = base
S = Sb(:) * ones(1, p.t_frame); % S = S of 1ms resolution
S = reshape(S', [ceil(p.t_max/p.t_frame)*p.t_frame, p.iters])';
S = S * p.k + p.k0;
Sb = (Sb - k)';

fprintf('Sensory matrix generated (%1.1fsec)\n', toc(sttime));
    %non-linear compression of signal
if isa(p.signal_compress, 'function_handle')
    S = p.signal_compress(S);
end

    %apply temporal weights
if p.stim_noise == 0
    warning('Stim noise is 0. The two race processes will be complitely anti-correlated.');
    noise_rho = -1;
else
    noise_rho = (p.rho * (1 + p.stim_noise^2) + 1)/(p.stim_noise^2); % see Race_model_correlation_math.docx
end
Noise = mvnrnd(zeros(length(S(:)),2),[1 noise_rho; noise_rho 1] * (p.stim_noise^2));
V1 = (S) .* repmat(p.weight, [p.iters 1]) + reshape(Noise(:,1), size(S));
V2 = (-S) .* repmat(p.weight, [p.iters 1]) + reshape(Noise(:,2), size(S));
clear Noise;
if ~p.get_raw_data
    clear S;
end
    %add the initial value
V1(:,1) = V1(:,1) + p.B0;
V2(:,1) = V2(:,1) + p.B0;

    %add decision noise
V1 = V1 + randn(size(V1)) * p.dec_noise;
V2 = V2 + randn(size(V2)) * p.dec_noise;


    %calculate the cumulative evidence
V1c = V1;
V2c = V2;
for t=2:size(V1,2)
    V1c(:,t) = V1c(:,t-1) + V1(:,t) - V2c(:,t-1) * p.interaction - V1c(:,t-1) * p.decay + p.booster;
    V2c(:,t) = V2c(:,t-1) + V2(:,t) - V1c(:,t-1) * p.interaction - V2c(:,t-1) * p.decay + p.booster;
    V1c(V1c(:,t) < rB(1), t) = rB(1);
    V2c(V2c(:,t) < rB(2), t) = rB(2);
end

fprintf('Decision variable generated (%1.1fsec)\n', toc(sttime));

clear V1 V2;
    %find the choice and bound crossing time for each trial 
choice = nan(p.iters, 1);
bound_crossing = zeros(p.iters, 1);
RT = nan(p.iters, 1); % reaction time including non-decision time
true_dec_time = nan(p.iters, 1); % the exact timing of decision (bound crossing)
est_dec_time = nan(p.iters, 1); % decision time estimated in kernel analysis (RT - p.subtract_time)

    %first find trials in which bound crossing took place
[I, J] = find(V1c>=repmat(p.bound_height(:,1)',[p.iters 1]) | V2c>=repmat(p.bound_height(:,2)',[p.iters 1]));
    %now find the first bound crossing time on these trials, ultimately we'd like to be
    %have an index to the trials with bound crossing. Also, J must be be the first time
    %of bound crossing on those trials 
if getMatlabVersion<8.3
    [I,ind] = sort(I, 1, 'ascend');
    I = flipud(I);
    J = flipud(J(ind));
    [I,ind] = unique(I);
    J = J(ind);
else
    [I,ind] = sort(I, 1, 'ascend');
    J = J(ind);
    [I,ind] = unique(I);
    J = J(ind);
end

    % determine choice & true decision time
bound_crossing(I) = 1;
L = V1c(sub2ind(size(V1c),I,J)) >= p.bound_height(J,1);
choice(I(L)) = 1;
L = V2c(sub2ind(size(V2c),I,J)) >= p.bound_height(J,2);
choice(I(L)) = 2;
true_dec_time(I) = J;

if strcmp(p.termination_rule{1}, 'Fixed') % in case of fixed duration task
    nI = find(bound_crossing==0);
    L = V1c(nI, p.termination_rule{2}) >= V2c(nI, p.termination_rule{2});
    choice(nI(L)) = 1;
    choice(nI(~L)) = 2;
    true_dec_time(nI) = p.termination_rule{2};
end
if p.cut_off_decision
    nI = find(bound_crossing==0);
    L = V1c(nI, end) >= V2c(nI, end);
    choice(nI(L)) = 1;
    choice(nI(~L)) = 2;
    true_dec_time(nI) = p.t_max;
end

    % determine RT  
if strcmp(p.termination_rule{1}, 'RT') % if RT task, RT is decision time + non decision time
    if strcmp(p.non_dec_time_dist, 'normal')
        RT = true_dec_time + p.non_dec_time + round(randn(size(RT)) * p.non_dec_time_sd);
    elseif strcmp(p.non_dec_time_dist, 'log normal')
        ev = @(x)abs(std(random('logn', ones(10000,1) * log(p.non_dec_time), ones(10000,1) * x)) - p.non_dec_time_sd);
        logsd = fminsearch(ev, 0.3);
        RT = true_dec_time + round(random('logn', ones(size(RT)) * log(p.non_dec_time), ones(size(RT)) * logsd));
    end
    RT(RT<=0) = 0; % clip at 0

    if ~p.cut_off_decision
        idx = RT >= p.t_max;
        RT(idx) = NaN; % remove trials whose decision time is larger than stimulus duration or smaller than 0.
        bound_crossing(idx) = 0;
    end
elseif strcmp(p.termination_rule{1}, 'Fixed') % if fixed duration task, RT is stimulus duration + non decision time
    if strcmp(p.non_dec_time_dist, 'normal')
        RT = p.termination_rule{2} + p.non_dec_time + round(randn(size(RT)) * p.non_dec_time_sd);
    elseif strcmp(p.non_dec_time_dist, 'log normal')
        ev = @(x)abs(std(random('logn', ones(10000,1) * log(p.non_dec_time), ones(10000,1) * x)) - p.non_dec_time_sd);
        logsd = fminsearch(ev, 0.3);
        RT = p.termination_rule{2} + round(random('logn', ones(size(RT)) * log(p.non_dec_time), ones(size(RT)) * logsd));
    end
end

est_dec_time(I) = RT(I) - p.subtract_time;

fprintf('RT computed (%1.1fsec)\n', toc(sttime));


    %determine what happens to decision variable and evidence traces after bound crossing 
for trial = find(bound_crossing)'
    if any(isnan(p.v_past_bound))
        %if you want to cut the path off after bound crossing 
        t = est_dec_time(trial)+p.include_dec_frame;
        if t < 1, t = 1; end
        if t <= size(V1c,2)
            V1c(trial, t:end) = NaN;
            V2c(trial, t:end) = NaN;
        end
    end
    if any(isnan(p.s_past_bound))
        t = ceil(est_dec_time(trial)/p.t_frame)+p.include_dec_frame;
        if t < 1, t = 1; end
        if t <= size(Sb,2)
            Sb(trial, t:end) = NaN;
        end
    end
end

%% kernel analysis
    %save the results
sim.choice = nansum(choice==2)/sum(~isnan(choice));
sim.medianRT = nanmedian(RT);
sim.sdRT = nanstd(RT);
    %average decision variable (path) aligned to stimulus onset
bt_median = floor(nanmedian(est_dec_time));
if isnan(bt_median) % no trial crossed the bound
    bt_median = p.t_max;
end
if isnan(p.cut_off_RT)
    bt_mediant = round(bt_median/p.t_frame);
else
    bt_mediant = p.cut_off_RT;
    bt_median = p.cut_off_RT * p.t_frame;
end


sim.path.motion{1} = nanmean(V1c(choice==1,1:bt_median),1);
sim.path.motion{2} = nanmean(V2c(choice==2,1:bt_median),1);
    %psychophysical kernel aligned to stimulus onset
sim.kernel.motion{1} = nanmean(Sb(choice==1,1:bt_mediant),1);
sim.kernel.motion{2} = nanmean(Sb(choice==2,1:bt_mediant),1);
    %average decision variable and psychophysical kernel aligned to choice 
V1_choice = nan(p.iters,bt_median);
V2_choice = nan(p.iters,bt_median);
S_choice = nan(p.iters,bt_mediant);
for trial = 1 : p.iters
    if bound_crossing(trial)==1
        valid_portion = max(0,est_dec_time(trial)-bt_median)+1:est_dec_time(trial);
        V1_choice(trial,end-length(valid_portion)+1:end) = V1c(trial,valid_portion);
        V2_choice(trial,end-length(valid_portion)+1:end) = V2c(trial,valid_portion);

        st = max(0,round(est_dec_time(trial)/p.t_frame) - bt_mediant)+1;
        en = min(floor(p.t_max/p.t_frame), round(est_dec_time(trial)/p.t_frame));
        valid_portion2 = st:en;
        S_choice(trial,end-length(valid_portion2)+1:end) = Sb(trial,valid_portion2);
    end
end
sim.bound_crossed_percent = sum(bound_crossing)/length(bound_crossing) * 100;
sim.path.choice{1} = nanmean(V1_choice(choice==1,:),1);
sim.path.choice{2} = nanmean(V2_choice(choice==2,:),1);
sim.kernel.choice{1} = nanmean(S_choice(choice==1,:),1);
sim.kernel.choice{2} = nanmean(S_choice(choice==2,:),1);
sim.t_frame = p.t_frame;
sim.cut_off_RT = bt_mediant;
sim.param = p;

if p.get_raw_data
    sim.Sb = Sb;
    sim.S = S;
    sim.V1c = V1c;
    sim.V2c = V2c;
    sim.bound_crossing = bound_crossing;
    sim.response = choice;
    sim.RT = RT;
    sim.est_dec_time = est_dec_time;
    sim.true_dec_time = true_dec_time;
    sim.coh = coh(1,:)'; % coherence level
end

fprintf('percent cross bound: %1.1f%%, median RT = %1.3f (time spent %1.1fsec)\n', sim.bound_crossed_percent, nanmedian(RT), toc(sttime));
if sim.bound_crossed_percent < 95
    if p.error_no_reach
        error('Bound crossed percent is less than 95%. You should use longer t_max');
    else
        warning('Bound crossed percent is less than 95%. You should use longer t_max');
    end
end


function sim = RACE_Kernel_Simulation_split_trial(p)
iters = p.iters;
max_tr = p.split_trials;
if mod(iters, max_tr)~=0
    error('iters cannot be split into split_trials');
end
N = iters/max_tr;

seed = p.seed;
if ~isnan(seed)
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed));
end
p.seed = NaN;

% compute N simulations
orig_iters = iters;
orig_split = p.split_trials;

p.split_trials = 0;
p.iters = max_tr;
sim_all = cell(N,1);
sim_all{1} = RACE_Kernel_Simulation(p);
p.cut_off_RT = sim_all{1}.cut_off_RT; % match RT across simulations
if p.parfor
    parfor n=2:N
        sim_all{n} = RACE_Kernel_Simulation(p);
    end
else
    for n=2:N
        sim_all{n} = RACE_Kernel_Simulation(p);
    end
end

% average N simulations
sim = sim_all{1};
for n=2:N
    sim.path.motion{1} = sim.path.motion{1} + sim_all{n}.path.motion{1};
    sim.path.motion{2} = sim.path.motion{2} + sim_all{n}.path.motion{2};
    sim.path.choice{1} = sim.path.choice{1} + sim_all{n}.path.choice{1};
    sim.path.choice{2} = sim.path.choice{2} + sim_all{n}.path.choice{2};
    sim.kernel.motion{1} = sim.kernel.motion{1} + sim_all{n}.kernel.motion{1};
    sim.kernel.motion{2} = sim.kernel.motion{2} + sim_all{n}.kernel.motion{2};
    sim.kernel.choice{1} = sim.kernel.choice{1} + sim_all{n}.kernel.choice{1};
    sim.kernel.choice{2} = sim.kernel.choice{2} + sim_all{n}.kernel.choice{2};
end
sim.path.motion{1} = sim.path.motion{1}/N;
sim.path.motion{2} = sim.path.motion{2}/N;
sim.path.choice{1} = sim.path.choice{1}/N;
sim.path.choice{2} = sim.path.choice{2}/N;
sim.kernel.motion{1} = sim.kernel.motion{1}/N;
sim.kernel.motion{2} = sim.kernel.motion{2}/N;
sim.kernel.choice{1} = sim.kernel.choice{1}/N;
sim.kernel.choice{2} = sim.kernel.choice{2}/N;

sim.param.iters = orig_iters;
sim.param.split_trials = orig_split;




