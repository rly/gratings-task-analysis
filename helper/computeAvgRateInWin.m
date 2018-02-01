function [avgRateInWin,sdRateInWin,spikeRateTrial] = computeAvgRateInWin(spikeTimes, win)
% compute the average firing rate in the window

% compute the average firing rate in the window specified by computing the
% total number of spikes in all trials and dividing by window length and 
% number of trials
%   spikeTimes - struct array of spike times. if empty, returns NaN
%   win - two-element vector of minimum and maximum time for computing the 
%   average rate, inclusive

if numel(win) ~= 2
    error('flanker_task_analysis:incorrectWindowFormat',...
            'Window must have two elements');
end

if win(2) <= win(1) 
    error('flanker_task_analysis:nonIncreasingWindow',...
            'Window [%0.3f,%0.3f] must be increasing', win(1), win(2));
end

nTrials = numel(spikeTimes);
nSpikesInWinTrial = nan(nTrials,1);
for j = 1:nTrials
    %if ~isempty(spikeTimes(j).times)
        nSpikesInWin = sum(win(1) <= spikeTimes(j).times & ...
                spikeTimes(j).times <= win(2));
        nSpikesInWinTrial(j) = nSpikesInWin;
    %end
end

spikeRateTrial = nSpikesInWinTrial / (win(2)-win(1));
avgRateInWin = mean(spikeRateTrial);
sdRateInWin = std(spikeRateTrial);
