%% setup
N = 3;
T = 10;
maxLag = T-1;

% this was for testing when my shuffle predictor was computed using
%     for k = 1:N and j ~= k
% now, SPNA is symmetric about 0 

%% no time-locked behavior
spikesInTrials = zeros(N, T);

% trial 1: spikes at every time point
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,:) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> all 0

% trials 1-2: spikes at every time point
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1:2,:) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> all 0.65

% trial 1: spikes at every time point + one spike in trial 2
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,:) = 1;
spikesInTrialsCopy(2,4) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> large peak at 1, above 0 at 1-6

% pairs of spikes, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,1:2) = 1;
spikesInTrialsCopy(2,3:4) = 1;
spikesInTrialsCopy(3,5:6) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1

% pairs of spikes, different times, reverse
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,5:6) = 1;
spikesInTrialsCopy(2,3:4) = 1;
spikesInTrialsCopy(3,1:2) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, same as above

% pairs of spikes, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,4:5) = 1;
spikesInTrialsCopy(2,6:7) = 1;
spikesInTrialsCopy(3,9:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, higher than above

% triplets of spikes, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,2:4) = 1;
spikesInTrialsCopy(2,5:7) = 1;
spikesInTrialsCopy(3,8:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, above 0 at 2

% triplets of spikes, mostly different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,3:5) = 1;
spikesInTrialsCopy(2,5:7) = 1;
spikesInTrialsCopy(3,7:9) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, less than above

% pairs of spikes, separated by a time point, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,[1 3]) = 1;
spikesInTrialsCopy(2,[4 6]) = 1;
spikesInTrialsCopy(3,[8 10]) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 2

% pairs of spikes, separated by a time point, overlapping windows
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,[1 3]) = 1;
spikesInTrialsCopy(2,[2 4]) = 1;
spikesInTrialsCopy(3,[5 7]) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 2

% pairs of spikes, separated by a time point, overlapping times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,[1 3]) = 1;
spikesInTrialsCopy(2,[3 5]) = 1;
spikesInTrialsCopy(3,[5 7]) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 2




%% 1st two time points have spikes on all trials -- time-locked behavior
spikesInTrials = zeros(N, T);
spikesInTrials(1:3,1:2) = 1;

spikesInTrialsCopy = spikesInTrials;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> all 0

% trial 1: spikes at every time point
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,:) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> large peak at 1, above 0 at 1-7
% compare to constant spiking in trial1 and no spiking in trial2, trial3 -> 0
% compare to spiking in first two time points only on all trials -> 0

% trial 1: spikes immediately after time-locked behavior
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,1:5) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, above 0 at 2

% trial 1: spikes sometime after time-locked behavior
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,8:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> all 0 
% compare with above

% trial 1: lots of spikes sometime after time-locked behavior
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,4:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 2, above 0 at 2-6

% pairs of spikes after time-locked behavior, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,4:5) = 1;
spikesInTrialsCopy(2,6:7) = 1;
spikesInTrialsCopy(3,9:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1

% triplets of spikes after time-locked behavior, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,3:5) = 1;
spikesInTrialsCopy(2,6:8) = 1;
spikesInTrialsCopy(3,8:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, above 0 at 2

% triplets of spikes after time-locked behavior, different times
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,4:6) = 1;
spikesInTrialsCopy(2,6:8) = 1;
spikesInTrialsCopy(3,8:10) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> peak at 1, less than above

%% last two time points have spikes on all trials -- time-locked behavior
spikesInTrials = zeros(N, T);
spikesInTrials(1:3,9:10) = 1;

spikesInTrialsCopy = spikesInTrials;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% -> all 0

% trial 1: spikes at every time point
spikesInTrialsCopy = spikesInTrials;
spikesInTrialsCopy(1,:) = 1;
runBurstTest(spikesInTrialsCopy, N, maxLag);
% peak at 1, above 0 1-7

%% random
spikesInTrials = zeros(N, T);

spikesInTrialsCopy = round(rand(N, T));
runBurstTest(spikesInTrialsCopy, N, maxLag);
