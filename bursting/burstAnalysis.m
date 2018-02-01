%% Load workspace

disp('-----------------------------------------------------')
clear


%% burst analysis
% a group of action potentials: 2 to 7 APs with interspike intervals <= 4
% ms and preceded by a silent period >= 100 ms. (Lu et al., 1992; Guido et
% al., 1995; Reinagel et al., 1999)

for i = 1:numel(S.spikeFiles5D)
    spikeFile = S.spikeFiles5D{i};
    sessionName = S.spikeFiles5D{i}(1:7);
    cellName = ['sig00' S.spikeFiles5D{i}(12:13)];

    allSpikesFileName = [dataDir sessionName '\' sessionName '_allSpikes.mat'];
    stimOnsetFileName = [dataDir sessionName '\StimuliOnset_5D.mat'];

    load(allSpikesFileName, cellName);
    load(stimOnsetFileName);
    
    burstSessionFigTitle = [makeSpikeFileTitle(spikeFile) ' -- Whole Session'];
    burstInTrialFigTitle = [makeSpikeFileTitle(spikeFile) ' -- Only Around Trials'];
    
    burstSessionSaveFile = [spikeFile(1:end-4) '_burstSession.png'];
    burstInTrialSaveFile = [spikeFile(1:end-4) '_burstInTrial.png'];
    
    spikeTimes = eval(cellName);
    Fs = 1000;
    minSilentPeriod = 0.02; % 0.1;
    minISIWithinBurst = 0.004;

    %% look at entire recording
    [diffSpikeTimesAll, silentPeriodDurationsPrecedingSpikeAll, numSpikesTotalAll, spikesInBurstAll, averageISIWithinBurstAll] = ...
            burstAnalysis2(spikeTimes, min(spikeTimes), max(spikeTimes), minSilentPeriod, minISIWithinBurst);

    diffSpikeTimes = diffSpikeTimesAll{1};
    fprintf('Note that there are %d spikes with ISI < 1ms.\n', sum(diffSpikeTimes < 0.001));

    plotBurstAnalysis(silentPeriodDurationsPrecedingSpikeAll, numel(spikeTimes), ...
            spikesInBurstAll, diffSpikeTimes, Fs, minSilentPeriod, ...
            minISIWithinBurst, burstSessionFigTitle, burstSessionSaveFile)
    close

    %% look in task periods

    % define some task-relevant periods before cue and after array onsets
    preCuePeriod = 0.2;
    postArrPeriod = 0.4; 

    allCueTimes = CueP(:);
    allCueTimes(isnan(allCueTimes)) = [];
    allCueTimes = sort(allCueTimes);
    allArrTimes = ArrP(:);
    allArrTimes(isnan(allArrTimes)) = [];
    allArrTimes = sort(allArrTimes);

    maxCueToArrInterval = 1; % really should be 0.7 but there is some leeway
    assert(numel(allCueTimes) == numel(allArrTimes));
    assert(max(allArrTimes - allCueTimes) < maxCueToArrInterval);

    % iterate through trials 
    startWindowTimes = allCueTimes - preCuePeriod;
    endWindowTimes = allArrTimes + postArrPeriod;

    [diffSpikeTimesInTrial, silentPeriodDurationsPrecedingSpikeInTrial, numSpikesTotalInTrial, spikesInBurstInTrial, averageISIWithinBurstInTrial] = ...
            burstAnalysis2(spikeTimes, startWindowTimes, endWindowTimes, minSilentPeriod, minISIWithinBurst);

    diffSpikeTimesInTrialAll = vertcat(diffSpikeTimesInTrial{:});


    plotBurstAnalysis(silentPeriodDurationsPrecedingSpikeInTrial, numSpikesTotalInTrial, ...
            spikesInBurstInTrial, diffSpikeTimesInTrialAll, Fs, minSilentPeriod, ...
            minISIWithinBurst, burstInTrialFigTitle, burstInTrialSaveFile)
    close
    
    %% compute burstiness using anderson et al 2011
    
    % need equal sized time periods
    preCuePeriodAnderson = 0.1;
    postCuePeriodAnderson = 0.6;
    % 700ms total
    spikeFs = 1000;
    T = (postCuePeriodAnderson - preCuePeriodAnderson) * spikeFs;
    N = numel(allCueTimes); % number of trials
    
    % generate binary trial x time matrix (1 = spike, 0 = not)
    spikesInTrials = zeros(N, T);
    spikeTimesInTrialAdj = cell(N, 1);
    
    for j = 1:N
        spikeTimesInTrial = spikeTimes(spikeTimes >= allCueTimes(j) - preCuePeriodAnderson + 1/spikeFs & ...
                spikeTimes <= allCueTimes(j) + postCuePeriodAnderson);
        spikeTimesInTrialAdj{j} = round((spikeTimesInTrial - ...
                (allCueTimes(j) - preCuePeriodAnderson)) * spikeFs);
        spikesInTrials(j, spikeTimesInTrialAdj{j}) = 1;
    end
    
    % compute auto-correlation with up to 10ms lag in either direction
    % TODO: use coeff scaling so that autocorr at zero lag = 1? this will
    % typically reduce all the values
    maxLag = spikeFs / 100; % 10ms
    lagT = (1:maxLag) * Fs / spikeFs;
    
    [SPNA, meanAutoCorr, meanCrossCorr, stdCrossCorr] = computeSPNA(spikesInTrials, maxLag);
    
    figure_tr_inch(18, 9);
    set(gcf,'Color','white');
    set(gcf,'renderer','painters');
    
    % raster plot by time
    subplot(2,3,[1 4]);
    hold on;
    rasterY = 0;
    plotParams = {'Color', [0 0 0], 'MarkerSize', 1.5, 'MarkerFaceColor', [0 0 0]};

    % plot diamond at each spike, one row (y-coord) per trial
    for j = 1:N
        rasterY = rasterY + 1;
        if ~isempty(spikeTimesInTrialAdj{j})
            plot(spikeTimesInTrialAdj{j} - preCuePeriodAnderson * spikeFs, ...
                    rasterY*ones(1, length(spikeTimesInTrialAdj{j})),...
                    'd', plotParams{:});
        end
    end
    title('Raster plot by time');
    xlabel('Time from cue onset (ms)');
    ylabel('Trial number');
    ylim([0 rasterY + 1]);

    subplot(2,3,2);
    plot(lagT, meanAutoCorr(maxLag+2:end));
    title('Auto-correlation');
    xlabel('Lag (ms)');
    
    subplot(2,3,3);
    plot(lagT, meanCrossCorr(maxLag+2:end));
    title('Mean cross-correlation');
    xlabel('Lag (ms)');
    
    subplot(2,3,5);
    plot(lagT, stdCrossCorr(maxLag+2:end));
    title('Stdev cross-correlation');
    xlabel('Lag (ms)');
    
    subplot(2,3,6);
    hold on;
    plot(lagT, SPNA(maxLag+2:end));
    xlabel('Lag (ms)');
    ylabel('Shuffle Predictor SD Units');
    title('Shuffle Predictor Normalized Autocorrelation');
    
    % BRI: burstiness-refractoriness index:
    % average of SPNA over lag of 1-4 ms
    BRILagStartInd = maxLag + 1 + spikeFs * 0.001;
    BRILagEndInd = maxLag + 1 + spikeFs * 0.004;
    BRI = mean(SPNA(BRILagStartInd:BRILagEndInd));
    
    currYLim = ylim;
    fillT = 1:Fs/spikeFs:4;
    jbfill(fillT, ones(size(fillT)) * currYLim(1), ...
        ones(size(fillT)) * currYLim(2), ...
        [0 1 0], [0 1 0], 1, 0.2);
    ylim(currYLim);
    
    hold on;
    plot(fillT, ones(size(fillT)) * BRI, 'g--');
    
    briFigTitle = sprintf('%s -- %d Trials -- BRI: %0.2f', ...
            makeSpikeFileTitle(spikeFile), numel(allCueTimes), BRI);
%     title(briFigTitle);
    
    briSaveFile = [spikeFile(1:end-4) '_BRI.png'];
    suptitle(briFigTitle);

    export_fig(briSaveFile, '-nocrop'); %, '-r300');
end
