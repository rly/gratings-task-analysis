function [SPNA, meanAutoCorr, meanCrossCorr, sdCrossCorr] = computeSPNA(spikesInTrials, maxLag)
% see Anderson, Mitchell, Reynolds. (2011). Attentional Modulation of
% Firing Rate Varies with Burstiness across Putative Pyramidal Neurons in
% Macaque Visual Area V4.

% by normalizing by the cross-correlation across trials, this measure 
% accounts for trial-locked fluctuations in spiking.

maxNumTrialsBeforeRandomizing = 100;

N = size(spikesInTrials, 1);
assert(N > 1);

autoCorrPerTrial = nan(N, 2 * maxLag + 1);
for j = 1:N
    autoCorrPerTrial(j,:) = xcorr(spikesInTrials(j,:), maxLag);
end
meanAutoCorr = mean(autoCorrPerTrial);
clear autoCorrPerTrial

if N <= maxNumTrialsBeforeRandomizing
    % compute shuffle predictor - every pair of trials O(N^2)
    crossCorrPerTrial = nan(N^2 - N, 2 * maxLag + 1);
    pairCounter = 1;
    for j = 1:N
        for k = 1:N
            if j == k
                continue;
            end
            if mod(pairCounter, 1000) == 0
                fprintf('\tshuffle %d/%d (%d%%)...\n', j, N^2-N, round(j/(N^2-N)*100));
            end
            crossCorrPerTrial(pairCounter, :) = xcorr(spikesInTrials(j,:), ...
                    spikesInTrials(k,:), maxLag);
            pairCounter = pairCounter + 1;
        end
    end
else
    % compute shuffle predictor - random pairs of trials for computational
    % efficiency
    crossCorrPerTrial = nan(maxNumTrialsBeforeRandomizing^2 - maxNumTrialsBeforeRandomizing, 2 * maxLag + 1);
    pairCounter = 1;
    for j = 1:size(crossCorrPerTrial, 1)
        if mod(j, 1000) == 0
            fprintf('\tshuffle %d/%d (%d%%)...\n', j, size(crossCorrPerTrial, 1), ...
                    round(j/size(crossCorrPerTrial, 1)*100));
        end
        randj = randi(N);
        randk = randj;
        while randk == randj % don't pick randk == randj
            randk = randi(N);
        end
        crossCorrPerTrial(pairCounter, :) = xcorr(spikesInTrials(randj,:), ...
                spikesInTrials(randk,:), maxLag);
        pairCounter = pairCounter + 1;
    end
end
meanCrossCorr = mean(crossCorrPerTrial);
sdCrossCorr = std(crossCorrPerTrial);
clear crossCorrPerTrial;

% adjust std so that values don't go to Inf TODO verify correctness
sdCrossCorrToNull = 1e-10;
sdCrossCorr(sdCrossCorr < sdCrossCorrToNull) = Inf;

% SPNA: shuffle-predictor-normalized autocorrelation
SPNA = (meanAutoCorr - meanCrossCorr) ./ sdCrossCorr;

% zero out the values that are clearly zero, for easier display
SPNAToNull = 1e-10;
SPNA(abs(SPNA) < SPNAToNull) = 0;
