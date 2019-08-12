% depending on randomness, coherencycpt finds the peak frequency of
% spike-field coherence. it is also quite jaggedy. be aware of peaks at
% 1/n frequency where n is small. regardless of noise condition, spectrum
% is quite flat between +/- 10 Hz of peak frequency

clear;

%% set up 
% chronux parameters
params.tapers = [4 7];
params.fpass = [5 100];
params.pad = 1;
params.Fs = 1000;
params.trialave = 1;

% data parameters
Fs = params.Fs;
f = 50;
nTime = 401;
nTrial = 10;
sigma = 0.1;
spikeSigma = 0.005;

%% create lfp data, sine wave with perfect frequency
t = 1:nTime;
y = nan(nTime, nTrial); 
r = rand(1, nTrial)*(Fs/f); % phase shift 
for i = 1:nTrial
    y(:,i) = sin((t'+r(i))/(Fs/f)*2*pi) + rand(nTime, 1)*2*sigma-sigma; % add gaussian noise
end

%% create spikes, one per cycle with jitter
s = struct();
for i = 1:nTrial
    s(i).times = (r(i) + 0:Fs/f:(nTime-r(i)))'/Fs; % based on phase shift
    s(i).times = s(i).times + rand(size(s(i).times))*spikeSigma;
end

%% plot one trial
figure;
hold on;

% plot LFP
plot(y(:,1));

% plot spikes
a = s(1).times;
for i = 1:numel(a)
    plot([a a]*Fs, [-1 1], 'LineWidth', 2);
end

%% plot spike-field coherence, mean across trials
[C,~,~,~,~,fAxis] = coherencycpt(y, s, params);

figure;
plot(fAxis, C);

%% %%%%%%%%%%%%%%%%%%%%%
% if the oscillation has perfect frequency but noisy amplitude, the
% coherence is more random, with a flat zone around the peak
% frequency, but still can have other high values, even if the spikes are
% perfectly oscillatory

%% set up 
% chronux parameters
params.tapers = [4 7];
params.fpass = [5 100];
params.pad = 1;
params.Fs = 1000;
params.trialave = 1;

% data parameters
Fs = params.Fs;
f = 50;
nTime = 401;
nTrial = 10;
sigma = 1; %%%
spikeSigma = 0; %%%

%% create lfp data, sine wave with perfect frequency, noisy amplitude
t = 1:nTime;
y = nan(nTime, nTrial); 
r = rand(1, nTrial)*(Fs/f); % phase shift 
for i = 1:nTrial
    y(:,i) = sin((t'+r(i))/(Fs/f)*2*pi) + rand(nTime, 1)*2*sigma-sigma; % add gaussian noise
end

%% create spikes, one per cycle with jitter
s = struct();
for i = 1:nTrial
    s(i).times = (r(i) + 0:Fs/f:(nTime-r(i)))'/Fs; % based on phase shift
    s(i).times = s(i).times + rand(size(s(i).times))*spikeSigma;
end

%% plot one trial
figure;
hold on;

% plot LFP
plot(y(:,1));

% plot spikes
a = s(1).times;
for i = 1:numel(a)
    plot([a a]*Fs, [-1 1], 'LineWidth', 2);
end

%% plot spike-field coherence, mean across trials
[C,~,~,~,~,fAxis] = coherencycpt(y, s, params);

figure;
plot(fAxis, C);

%% %%%%%%%%%%%%%%%%%%%%%
% if there is no underlying oscillation, the coherence is more or less
% random, but still can have high values, even if the spikes are perfectly
% oscillatory

%% set up 
% chronux parameters
params.tapers = [4 7];
params.fpass = [5 100];
params.pad = 1;
params.Fs = 1000;
params.trialave = 1;

% data parameters
Fs = params.Fs;
f = 50;
nTime = 401;
nTrial = 10;
sigma = 0.1;
spikeSigma = 0; %%%

%% create lfp data, random noise
t = 1:nTime;
y = rand(nTime, nTrial)*2-1; %%%
r = rand(1, nTrial)*(Fs/f); % phase shift 

%% create spikes, one per cycle with jitter
s = struct();
for i = 1:nTrial
    s(i).times = (r(i) + 0:Fs/f:(nTime-r(i)))'/Fs; % based on phase shift
    s(i).times = s(i).times + rand(size(s(i).times))*spikeSigma;
end

%% plot one trial
figure;
hold on;

% plot LFP
plot(y(:,1));

% plot spikes
a = s(1).times;
for i = 1:numel(a)
    plot([a a]*Fs, [-1 1], 'LineWidth', 2);
end

%% plot spike-field coherence, mean across trials
[C,~,~,~,~,fAxis] = coherencycpt(y, s, params);

figure;
plot(fAxis, C);
