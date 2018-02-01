%% setup
N = 3;
T = 10;
maxLag = T-1;

% this was for testing when my shuffle predictor was computed using
%     for k = j+1:N
% but cross-correlation is NOT symmetric so i need to run the other side
% too for a proper mean and SD.

%% no time-locked behavior
spikesInTrials = zeros(N, T);

% trial 1: spikes at every time point
% spikesInTrials(1,:) = 1;
% -> all 0

% trials 1-2: spikes at every time point
% spikesInTrials(1:2,:) = 1;
% -> all 0.58

% pairs of spikes, different times
% spikesInTrials(1,1:2) = 1;
% spikesInTrials(2,3:4) = 1;
% spikesInTrials(3,5:6) = 1;
% -> peak at lag = -1

% pairs of spikes, different times, reverse
spikesInTrials(1,5:6) = 1;
spikesInTrials(2,3:4) = 1;
spikesInTrials(3,1:2) = 1;
% -> peak at lag = 1

% pairs of spikes, different times
% spikesInTrials(1,4:5) = 1;
% spikesInTrials(2,6:7) = 1;
% spikesInTrials(3,9:10) = 1;
% -> peak at lag = -1

% triplets of spikes, different times
% spikesInTrials(1,2:4) = 1;
% spikesInTrials(2,5:7) = 1;
% spikesInTrials(3,8:10) = 1;
% -> peak at lag = -1

% triplets of spikes, mostly different times
% spikesInTrials(1,3:5) = 1;
% spikesInTrials(2,5:7) = 1;
% spikesInTrials(3,7:9) = 1;
% -> peak at lag = -1

% pairs of spikes, separated by a time point, different times
% spikesInTrials(1,[1 3]) = 1;
% spikesInTrials(2,[4 6]) = 1;
% spikesInTrials(3,[8 10]) = 1;
% -> peak at lag = -2

% pairs of spikes, separated by a time point, overlapping windows
% spikesInTrials(1,[1 3]) = 1;
% spikesInTrials(2,[2 4]) = 1;
% spikesInTrials(3,[5 7]) = 1;
% -> peak at lag = -2

% pairs of spikes, separated by a time point, overlapping times
% spikesInTrials(1,[1 3]) = 1;
% spikesInTrials(2,[3 5]) = 1;
% spikesInTrials(3,[5 7]) = 1;
% -> no peak (other than lag=0)




%% 1st two time points have spikes on all trials -- time-locked behavior
% spikesInTrials = zeros(N, T);
% spikesInTrials(1:3,1:2) = 1;
% -> all 0

% trial 1: spikes at every time point
% spikesInTrials(1,:) = 1;
% -> peak at lag = 1

% trial 1: spikes immediately after time-locked behavior
% spikesInTrials(1,1:5) = 1;
% -> peak at lag = 1

% trial 1: spikes sometime after time-locked behavior
% spikesInTrials(1,8:10) = 1;
% -> all 0 or negative

% trial 1: lots of spikes sometime after time-locked behavior
% spikesInTrials(1,4:10) = 1;
% -> peak at lag = 2, small peak at lag=3,4, but 0 or negative
% after

% pairs of spikes after time-locked behavior, different times
% spikesInTrials(1,4:5) = 1;
% spikesInTrials(2,6:7) = 1;
% spikesInTrials(3,9:10) = 1;
% -> peak at lag = -1

% triplets of spikes after time-locked behavior, different times
% spikesInTrials(1,3:5) = 1;
% spikesInTrials(2,6:8) = 1;
% spikesInTrials(3,8:10) = 1;
% -> above 0 at lag = -6, -5, -1, 0, 1, 2, 3, 4, 7. peak at -1

%% last two time points have spikes on all trials -- time-locked behavior
% spikesInTrials = zeros(N, T);
% spikesInTrials(1:3,9:10) = 1;
% -> all 0

% trial 1: spikes at every time point
% spikesInTrials(1,:) = 1;
% -> peak at lag = -1

%% compute shuffle-predictor normalized autocorrelation
% computeSPNA(spikesInTrials, N, maxLag);

computeSPNA;
