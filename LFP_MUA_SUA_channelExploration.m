clear all; close all
processedDataRootDir = '/Volumes/scratch/rly/gratings-task-analysis/processed_data/';
dataDirRoot = '/Volumes/kastner/ryanly/Ferdy/merged';
muaDataDirRoot = '/Volumes/scratch/rly/simple-mua-detection/processed_data/';
recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
recordingInfo = readRecordingInfo(recordingInfoFileName);

for sessioni = 38:numel(recordingInfo)
    sessionInd = sessioni; %12 18 15 22
    if ~isnan(recordingInfo(sessioni).aepmIndices)        
        channelsToLoad = recordingInfo(sessioni).spikeChannelsToLoad;
        [Raud, Daud, processedDataDirAud, blockNameAud] = loadRecordingData(...
                processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
                sessionInd, channelsToLoad, 'AEPM', 'LFP_AEPM', 1, 1, 1, 0);
        origAudioEvents = Daud.events{3};
        preAudioEvents = Daud.events{2};
        
        sessionName = Raud.sessionName;
        areaName = Raud.areaName;
        
        periAudioWindowOffset = [-0.25 0.3]; % seconds around sound
        hiCutoffFreq = 100;

        Fs = Daud.lfpFs;
        Daud.adjLfpsClean = Daud.adjLfps; 
        [channelDataCARNormAud,channelDataNormAud,commonAverageNormAud,isNoisyChannelAud] = preprocessLfps(...
                Daud.adjLfpsClean, Fs, Daud.lfpNames, processedDataDirAud, [], hiCutoffFreq, 1, []);
        Daud.adjLfpsClean = [];
        audioEventsClean = detectOutlierLfpEvents(channelDataCARNormAud, Fs, Daud.lfpNames, origAudioEvents, ...
                periAudioWindowOffset, processedDataDirAud, [], 1, []);
        
        nAudio = numel(audioEventsClean);
        
        ref = 'CAR';
        fprintf('\nExtracting sound-aligned responses, ref: %s...\n', ref);

        startIndices = round((audioEventsClean + periAudioWindowOffset(1)) * Fs); % time to index conversion
        endIndices = startIndices + round(diff(periAudioWindowOffset) * Fs) - 1;
        t = periAudioWindowOffset(1):1/Fs:periAudioWindowOffset(2)-1/Fs;
        nTime = numel(t);
        assert(all(nTime == (endIndices - startIndices + 1)));
        
        nChannels = 32;
        responses = nan(nChannels, nTime, nAudio);
        for j = 1:nChannels
            % can vectorize???
            for i = 1:nAudio
                if strcmp(ref, 'RAW')
                    responses(j,:,i) = channelDataNorm(j,startIndices(i):endIndices(i));
                elseif strcmp(ref, 'CAR')
                    responses(j,:,i) = channelDataCARNormAud(j,startIndices(i):endIndices(i));
                elseif strcmp(ref, 'BIP')
                    responses(j,:,i) = channelDataBIPNorm(j,startIndices(i):endIndices(i));
                end            
            end
        end
        %% subtract out baseline
        % units are standard deviations from baselined mean
        averageResponse = mean(responses, 3); % average across flashes

        % subtract mean baseline activity (-0.1, 0] seconds before audio
        audioBaselineWindowOffset = [-0.1 0];

        % TODO make so that baseline can have arbitrary end time
        indexAudioTime = -round(periAudioWindowOffset(1) * Fs);
        indexStartBaselineAudioTime = -round((periAudioWindowOffset(1) - audioBaselineWindowOffset(1)) * Fs) + 1;
        for j = 1:nChannels
            averageResponse(j,:) = averageResponse(j,:) - mean(averageResponse(j,indexStartBaselineAudioTime:indexAudioTime));
        end

        %% staggered channel line plot of average visually evoked LFP
        xBounds = [-0.1 0.3];

        responsePlotYOffset = 1;
        responsePlotYScale = 5/max(max(abs(averageResponse)));

        figure_tr_inch(8, 10);
        hold on;
        maxYLim = 3;
        minYLim = -(nChannels + 2);
        for j = 1:nChannels
            plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
            if j == 1
                maxYLim = max([maxYLim max(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) + 1]);
            elseif j == nChannels
                minYLim = min([minYLim min(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) - 1]);
            end
        end
        xlim(periAudioWindowOffset);
        set(gca, 'YTickMode', 'manual');
        set(gca, 'YTick', -1*nChannels:-1);
        set(gca, 'YTickLabelMode', 'manual');
        set(gca, 'YTickLabel', nChannels:-1:1);
        set(gca, 'FontSize', 16);
        ax = ancestor(gca, 'axes');
        ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
        ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
        xlabel('Time from Sound Onset (s)');
        title(sprintf('%s %s - Response to Audio Mapping (N=%d) (%s)', sessionName, areaName, nAudio, ref));

        origYLim = [minYLim maxYLim]; % ylim();
        plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
        xlim(xBounds);
        ylim(origYLim);

        % plot early latency line at 35 ms, 45 ms
        plot([0.035 0.035], [-1000 1000], 'm-');
        plot([0.045 0.045], [-1000 1000], 'm-');
        text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');

         plotFileName = sprintf('%s-%s-%s-AEPMlfpLines.png', recordingInfo(sessioni).sessionName, sessioni, ref);
        
        cd('/Users/labmanager/Documents/MATLAB/aepm_lfp')
        saveas(gcf,plotFileName,'jpg')
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis')
        
        %% staggered channel color plot of average visually evoked LFP
        figure_tr_inch(8, 10);
        subaxis(1, 1, 1, 'ML', 0.1);
        hold on;
        imagesc(t, 1:nChannels, averageResponse);
        set(gca, 'YDir', 'reverse');
        plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
        xlim(xBounds);
        ylim([0.5 nChannels+0.5]);
        xlabel('Time from Sound Onset (s)');
        ylabel('Channel Number (1 = topmost)');
        title(sprintf('%s %s - Response to Audio Mapping (N=%d) (%s)', sessionName, areaName, nAudio, ref));
        maxCAxis = max(abs(caxis));
        caxis([-maxCAxis maxCAxis]);
        colormap(getCoolWarmMap());
        colorbar;
        set(gca, 'FontSize', 16);

        if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
            plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
        end

        % plot early latency line at 35 ms, 45 ms
        plot([0.035 0.035], [-1000 1000], 'm-');
        plot([0.045 0.045], [-1000 1000], 'm-');
        text(0.050, nChannels, '35-45 ms', 'Color', 'm');

        %% plot boxes around areas abs() > thresh over last plot and re-save
        boxAbsThresh = 0.25;
        strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
        strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
        for k = 1:length(strongResponseGroupBoundaries)
           boundary = strongResponseGroupBoundaries{k};
           plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)

           minBY = min(boundary(:,1));
           maxBY = max(boundary(:,1));
           text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
           text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
        end

        plotFileName = sprintf('%s-%s-%s-AEPMlfpColorBounds.png', recordingInfo(sessioni).sessionName, sessioni, ref);
        
        cd('/Users/labmanager/Documents/MATLAB/aepm_lfp')
        saveas(gcf,plotFileName,'jpg')
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis')

    end
end
%[fixationAndLeverTimes,isMissingData] = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc)

%%
% flash related activity 200ms windows
postFlashWindowOffset = [0.025 0.225]; % seconds after flash
baselineWindowOffset = [-.175 .025]; % seconds around preflashevent
saccadeWindowOffset = [-.425 -.225]; % fixation onset at 325ms before preflashevent
postAudioWindowOffset = postFlashWindowOffset;

% preprocess LFPs
Fs = D.lfpFs;
D.adjLfpsClean = D.adjLfps; %interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);
hiCutoffFreq = 100;
[channelDataCARNorm,channelDataNorm,commonAverageNorm,isNoisyChannel] = preprocessLfps(...
        D.adjLfpsClean, Fs, D.lfpNames, processedDataDir, [], hiCutoffFreq, 1, []);
D.adjLfpsClean = [];
flashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, origFlashEvents, ...
        postFlashWindowOffset, processedDataDir, [], 1, []);
preFlashEventsCleanWindowOffset = [min([baselineWindowOffset saccadeWindowOffset]) max([baselineWindowOffset saccadeWindowOffset])];
preFlashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, preFlashesEvents, ...
        preFlashEventsCleanWindowOffset, processedDataDir, [], 1, []);

Daud.adjLfpsClean = Daud.adjLfps; 
[channelDataCARNormAud,channelDataNormAud,commonAverageNormAud,isNoisyChannelAud] = preprocessLfps(...
        Daud.adjLfpsClean, Fs, Daud.lfpNames, processedDataDir, [], hiCutoffFreq, 1, []);
Daud.adjLfpsClean = [];
audioEventsClean = detectOutlierLfpEvents(channelDataCARNormAud, Fs, Daud.lfpNames, origAudioEvents, ...
        postAudioWindowOffset, processedDataDir, [], 1, []);
preAudioEventsCleanWindowOffset = baselineWindowOffset;
preAudioEventsClean = detectOutlierLfpEvents(channelDataCARNormAud, Fs, Daud.lfpNames, preAudioEvents, ...
        preAudioEventsCleanWindowOffset, processedDataDir, [], 1, []);
      
nFlashes = numel(flashEventsClean);
nTrials = numel(preFlashEventsClean);
startIndicesStim = round((flashEventsClean + postFlashWindowOffset(1)) * Fs); % time to index conversion
endIndicesStim = startIndicesStim + round(diff(postFlashWindowOffset) * Fs) - 1;
startIndicesBl = round((preFlashEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
endIndicesBl = startIndicesBl + round(diff(baselineWindowOffset) * Fs) - 1;
startIndicesSac = round((preFlashEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
endIndicesSac = startIndicesSac + round(diff(saccadeWindowOffset) * Fs) - 1;
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
% assert(all(nTime == (endIndicesStim - startIndicesStim + 1)));

% to look at Raw data
% adjLfps = D.adjLfps;
% isNanOrig = isnan(adjLfps);
% bFirLowPass = fir1(3*fix(Fs/hiCutoffFreq), hiCutoffFreq/(Fs/2), 'low');
% adjLfps(isnan(adjLfps)) = 0; % zero out nans
% adjLfpsLP = filtfilt(bFirLowPass, 1, adjLfps')';
% adjLfpsLP(isNanOrig) = NaN; % set missing vals back to nan

nChannels = length(channelsToLoad);
responses = nan(nChannels, nTime, nFlashes);
baseline = nan(nChannels, nTime, nTrials);
saccade = nan(nChannels, nTime, nTrials);
% responsesRaw = nan(nChannels, nTime, nFlashes);
% baselineRaw = nan(nChannels, nTime, nTrials);
% saccadeRaw = nan(nChannels, nTime, nTrials);
for j = 1:nChannels
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataCARNorm(j,startIndicesStim(i):endIndicesStim(i));
%         responsesRaw(j,:,i) = adjLfpsLP(j,startIndicesStim(i):endIndicesStim(i));
    end
    for k = 1:nTrials
        baseline(j,:,k) = channelDataCARNorm(j,startIndicesBl(k):endIndicesBl(k));
        saccade(j,:,k) = channelDataCARNorm(j,startIndicesSac(k):endIndicesSac(k));
%         baselineRaw(j,:,k) = adjLfpsLP(j,startIndicesBl(k):endIndicesBl(k));
%         saccadeRaw(j,:,k) = adjLfpsLP(j,startIndicesSac(k):endIndicesSac(k));
    end
end
allResponses = cat(3,responses,saccade);

% repeat for auditory task
nAudio = numel(audioEventsClean);
nTrialsAud = numel(preAudioEventsClean);
startIndicesStimAud = round((audioEventsClean + postAudioWindowOffset(1)) * Fs); % time to index conversion
endIndicesStimAud = startIndicesStimAud + round(diff(postAudioWindowOffset) * Fs) - 1;
startIndicesBlAud = round((preAudioEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
endIndicesBlAud = startIndicesBlAud + round(diff(baselineWindowOffset) * Fs) - 1;
startIndicesSacAud = round((preAudioEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
endIndicesSacAud = startIndicesSacAud + round(diff(saccadeWindowOffset) * Fs) - 1;
tAud = postAudioWindowOffset(1):1/Fs:postAudioWindowOffset(2)-1/Fs;
nTimeAud = numel(tAud);

nChannels = length(channelsToLoad);
responsesAud = nan(nChannels, nTimeAud, nAudio);
baselineAud = nan(nChannels, nTimeAud, nTrialsAud);
saccadeAud = nan(nChannels, nTimeAud, nTrialsAud);
for j = 1:nChannels
    % can vectorize???
    for i = 1:nAudio
        responsesAud(j,:,i) = channelDataCARNormAud(j,startIndicesStimAud(i):endIndicesStimAud(i));
    end
    for k = 1:nTrialsAud
        baselineAud(j,:,k) = channelDataCARNormAud(j,startIndicesBlAud(k):endIndicesBlAud(k));
        saccadeAud(j,:,k) = channelDataCARNormAud(j,startIndicesSacAud(k):endIndicesSacAud(k));
    end
end

% plot LFP for flashes, baseline and saccade of flash mapping task
ref = 'CAR';
epochType = {'Full Flashes'; 'Baseline FM'; 'Saccade FM'};
for condi = 1:3
    if condi == 1
        averageResponse = mean(responses, 3); % average across flashes
        t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
    elseif condi == 2
        averageResponse = mean(baseline, 3);
        t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
    else
        averageResponse = mean(saccade, 3);
        t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
    end

    xBounds = [-0.5 0.3];

    responsePlotYOffset = 1;
    responsePlotYScale = 5/max(max(abs(averageResponse)));

    % staggered channel color plot of average visually evoked LFP
    figure_tr_inch(8, 10);
    subaxis(1, 1, 1, 'ML', 0.1);
    hold on;
    imagesc(t, 1:nChannels, averageResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
    xlim(xBounds);
    ylim([0.5 nChannels+0.5]);
    xlabel('Time from Flash Onset (s)');
    ylabel('Channel Number (1 = topmost)');
    title(sprintf(['%s %s - Response to ' epochType{condi} ' (N=%d) (%s)'], sessionName, areaName, nFlashes, ref));
    maxCAxis = max(abs(caxis));
    caxis([-maxCAxis maxCAxis]);
    colormap(getCoolWarmMap());
    colorbar;
    set(gca, 'FontSize', 16);
    if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
        plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
    end

    % plot early latency line at 35 ms, 45 ms
    plot([0.035 0.035], [-1000 1000], 'm-');
    plot([0.045 0.045], [-1000 1000], 'm-');
    text(0.050, nChannels, '35-45 ms', 'Color', 'm');

    %% plot boxes around areas abs() > thresh over last plot and re-save
    boxAbsThresh = 0.25;
    strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
    strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
    for k = 1:length(strongResponseGroupBoundaries)
       boundary = strongResponseGroupBoundaries{k};
       plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)

       minBY = min(boundary(:,1));
       maxBY = max(boundary(:,1));
       text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
       text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
    end
    cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
    saveas(gcf,[sessionName '_LFP_' epochType{condi}],'jpg')
    cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')
end

% plot LFP for tones, baseline and saccade of auditory mapping task
ref = 'CAR';
epochType = {'Tones'; 'Baseline AM'; 'Saccade AM'};
for condi = 1:3
    if condi == 1
        averageResponse = mean(responsesAud, 3); % average across flashes
        t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
    elseif condi == 2
        averageResponse = mean(baselineAud, 3);
        t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
    else
        averageResponse = mean(saccadeAud, 3);
        t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
    end

    xBounds = [-0.5 0.3];

    responsePlotYOffset = 1;
    responsePlotYScale = 5/max(max(abs(averageResponse)));

    % staggered channel color plot of average visually evoked LFP
    figure_tr_inch(8, 10);
    subaxis(1, 1, 1, 'ML', 0.1);
    hold on;
    imagesc(t, 1:nChannels, averageResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
    xlim(xBounds);
    ylim([0.5 nChannels+0.5]);
    xlabel('Time from Flash Onset (s)');
    ylabel('Channel Number (1 = topmost)');
    title(sprintf(['%s %s - Response to ' epochType{condi} ' (N=%d) (%s)'], sessionName, areaName, nFlashes, ref));
    maxCAxis = max(abs(caxis));
    caxis([-maxCAxis maxCAxis]);
    colormap(getCoolWarmMap());
    colorbar;
    set(gca, 'FontSize', 16);
    if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
        plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
    end

    % plot early latency line at 35 ms, 45 ms
    plot([0.035 0.035], [-1000 1000], 'm-');
    plot([0.045 0.045], [-1000 1000], 'm-');
    text(0.050, nChannels, '35-45 ms', 'Color', 'm');

    %% plot boxes around areas abs() > thresh over last plot and re-save
    boxAbsThresh = 0.25;
    strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
    strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
    for k = 1:length(strongResponseGroupBoundaries)
       boundary = strongResponseGroupBoundaries{k};
       plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)

       minBY = min(boundary(:,1));
       maxBY = max(boundary(:,1));
       text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
       text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
    end
    cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
    saveas(gcf,[sessionName '_LFP_' epochType{condi}],'jpg')
    cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')
end

% extract power for each trial, than average
responsesPow = nan(size(responses,3),1025,nChannels);
for triali = 1:size(responses, 3)
    [responsesFreq,responsesPow(triali,:,:)] = PSD(responses(:,:,triali)');
end
responsesIndPow = nan(size(responses,3),1025,nChannels);
for triali = 1:size(responses, 3)
    [responsesFreq,responsesIndPow(triali,:,:)] = PSD((responses(:,:,triali)-mean(responses,3))');
end
baselinePow = nan(size(baseline,3),1025,nChannels);
for triali = 1:size(baseline, 3)
    [baselineFreq,baselinePow(triali,:,:)] = PSD(baseline(:,:,triali)');
end
saccadePow = nan(size(saccade,3),1025,nChannels);
for triali = 1:size(saccade, 3)
    [saccadeFreq,saccadePow(triali,:,:)] = PSD(saccade(:,:,triali)');
end

responsesAudPow = nan(size(responsesAud,3),1025,nChannels);
for triali = 1:size(responsesAud, 3)
    [responsesFreq,responsesAudPow(triali,:,:)] = PSD(responsesAud(:,:,triali)');
end
responsesAudIndPow = nan(size(responsesAud,3),1025,nChannels);
for triali = 1:size(responsesAud, 3)
    [responsesFreq,responsesAudIndPow(triali,:,:)] = PSD((responsesAud(:,:,triali)-mean(responsesAud,3))');
end
baselineAudPow = nan(size(baselineAud,3),1025,nChannels);
for triali = 1:size(baselineAud, 3)
    [baselineFreq,baselineAudPow(triali,:,:)] = PSD(baselineAud(:,:,triali)');
end
saccadeAudPow = nan(size(saccadeAud,3),1025,nChannels);
for triali = 1:size(saccadeAud, 3)
    [saccadeFreq,saccadeAudPow(triali,:,:)] = PSD(saccadeAud(:,:,triali)');
end

figure; 
subplot(141)
imagesc(squeeze(mean(responses,3)))
caxis([-2.5 2.5])
title('Full-flashes')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
subplot(143)
imagesc(squeeze(mean(baseline,3)))
caxis([-2.5 2.5])
title('Baseline')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
subplot(144)
imagesc(squeeze(mean(saccade,3)))
caxis([-2.5 2.5])
title('Saccade')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

figure; 
subplot(141)
imagesc(squeeze(mean(responsesAud,3)))
caxis([-2.5 2.5])
title('Audio')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
subplot(143)
imagesc(squeeze(mean(baselineAud,3)))
caxis([-2.5 2.5])
title('Baseline AM')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
subplot(144)
imagesc(squeeze(mean(saccadeAud,3)))
caxis([-2.5 2.5])
title('Saccade AM')
ylabel('Electrodes (1 = topmost)'); xlabel('Time')
colormap(getCoolWarmMap());
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% calculate CSD Ryan's way
yVals = 2:nChannels-1;
xBounds = [0.025 0.225];
yBounds = yVals([1 end]) + [-0.5 0.5];
cBounds = [-0.5 0.5];
meanFlashOnsetLfp = squeeze(mean(responses, 3));
flashOnsetCSD = nan(numel(yVals), size(meanFlashOnsetLfp, 2));
meanBaselineLfp = squeeze(mean(baseline, 3));
baselineCSD = nan(numel(yVals), size(meanBaselineLfp, 2));
meanSaccadeLfp = squeeze(mean(saccade, 3));
saccadeCSD = nan(numel(yVals), size(meanSaccadeLfp, 2));
meanAudioOnsetLfp = squeeze(mean(responsesAud, 3));
audioOnsetCSD = nan(numel(yVals), size(meanAudioOnsetLfp, 2));
meanBaselineAudLfp = squeeze(mean(baselineAud, 3));
baselineAudCSD = nan(numel(yVals), size(meanBaselineAudLfp, 2));
meanSaccadeAudLfp = squeeze(mean(saccadeAud, 3));
saccadeAudCSD = nan(numel(yVals), size(meanSaccadeAudLfp, 2));
for j = 1:numel(yVals)
    ji = yVals(j);
    flashOnsetCSD(j,:) = meanFlashOnsetLfp(ji+1,:) - 2 * meanFlashOnsetLfp(ji,:) + meanFlashOnsetLfp(ji-1,:);
    baselineCSD(j,:) = meanBaselineLfp(ji+1,:) - 2 * meanBaselineLfp(ji,:) + meanBaselineLfp(ji-1,:);
    saccadeCSD(j,:) = meanSaccadeLfp(ji+1,:) - 2 * meanSaccadeLfp(ji,:) + meanSaccadeLfp(ji-1,:);
    audioOnsetCSD(j,:) = meanAudioOnsetLfp(ji+1,:) - 2 * meanAudioOnsetLfp(ji,:) + meanAudioOnsetLfp(ji-1,:);
    baselineAudCSD(j,:) = meanBaselineAudLfp(ji+1,:) - 2 * meanBaselineAudLfp(ji,:) + meanBaselineAudLfp(ji-1,:);
    saccadeAudCSD(j,:) = meanSaccadeAudLfp(ji+1,:) - 2 * meanSaccadeAudLfp(ji,:) + meanSaccadeAudLfp(ji-1,:);
end

t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
figure
subplot(141)
imagesc(t, yVals, flashOnsetCSD);
ylim(yBounds);
caxis(cBounds);
colormap(getCoolWarmMap());
subplot(143)
t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
imagesc(t, yVals, baselineCSD);
ylim(yBounds);
caxis(cBounds);
subplot(144)
t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
imagesc(t, yVals, saccadeCSD);
ylim(yBounds);
caxis(cBounds);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_CSD_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
figure
subplot(141)
imagesc(t, yVals, audioOnsetCSD);
ylim(yBounds);
caxis(cBounds);
colormap(getCoolWarmMap());
subplot(143)
t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
imagesc(t, yVals, baselineAudCSD);
ylim(yBounds);
caxis(cBounds);
subplot(144)
t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
imagesc(t, yVals, saccadeAudCSD);
ylim(yBounds);
caxis(cBounds);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_CSD_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% use mod_iCSD from Kacie to get at CSD
flashOnsetCSDKacie = mod_iCSD(squeeze(mean(responses,3)));
baselineCSDKacie = mod_iCSD(squeeze(mean(baseline,3)));
saccadeCSDKacie = mod_iCSD(squeeze(mean(saccade,3)));
audioOnsetCSDKacie = mod_iCSD(squeeze(mean(responsesAud,3)));
baselineAudCSDKacie = mod_iCSD(squeeze(mean(baselineAud,3)));
saccadeAudCSDKacie = mod_iCSD(squeeze(mean(saccadeAud,3)));
cBounds = [-5 5];
figure
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
subplot(141)
imagesc(t, yVals, flashOnsetCSDKacie);
ylim(yBounds);
caxis(cBounds);
colormap(getCoolWarmMap());
subplot(143)
t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
imagesc(t, yVals, baselineCSDKacie);
ylim(yBounds);
caxis(cBounds);
subplot(144)
t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
imagesc(t, yVals, saccadeCSDKacie);
ylim(yBounds);
caxis(cBounds);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_CSDKacie_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

cBounds = [-12 12];
figure
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
subplot(141)
imagesc(t, yVals, audioOnsetCSDKacie);
ylim(yBounds);
caxis(cBounds);
colormap(getCoolWarmMap());
subplot(143)
t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2)-1/Fs;
imagesc(t, yVals, baselineAudCSDKacie);
ylim(yBounds);
caxis(cBounds);
subplot(144)
t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;
imagesc(t, yVals, saccadeAudCSDKacie);
ylim(yBounds);
caxis(cBounds);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_AllEpochs_CSDKacie_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% plot power of LFP, CSD from Ryan and CSD from Kacie
figure; 
subplot(141)
imagesc(squeeze(mean(responsesPow,1))')
caxis([0 2.5])
xlim([0 150])
title('Full-flashes')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(142)
imagesc(squeeze(mean(responsesIndPow,1))')
caxis([0 2.5])
xlim([0 150])
title('Full-flashes Induced')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(143)
imagesc(squeeze(mean(baselinePow,1))')
caxis([0 2.5])
xlim([0 150])
title('baseline FM')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(144)
imagesc(squeeze(mean(saccadePow,1))')
caxis([0 2.5])
xlim([0 150])
title('saccade FM')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')
    
figure; 
subplot(141)
imagesc(squeeze(mean(responsesAudPow,1))')
caxis([-100 100])
xlim([0 150])
title('Audio')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(142)
imagesc(squeeze(mean(responsesAudIndPow,1))')
caxis([-100 100])
xlim([0 150])
title('Audio Induced')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(143)
imagesc(squeeze(mean(baselineAudPow,1))')
caxis([-100 100])
xlim([0 150])
title('baseline AM')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
subplot(144)
imagesc(squeeze(mean(saccadeAudPow,1))')
caxis([-100 100])
xlim([0 150])
title('saccade AM')
ylabel('Electrodes (1 = topmost)'); xlabel('Frequency (Hz)')
colormap(getCoolWarmMap());
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')


% extract power from CSD data
[responsesFreq,responsesPowAvg] = PSD(mean(responses,3)');
[responsesFreq,baselinePowAvg] = PSD(mean(baseline,3)');
[responsesFreq,saccadePowAvg] = PSD(mean(saccade,3)');
[responsesFreq,responsesCSDPowAvg] = PSD(mean(flashOnsetCSD,3)');
[responsesFreq,baselineCSDPowAvg] = PSD(mean(baselineCSD,3)');
[responsesFreq,saccadeCSDPowAvg] = PSD(mean(saccadeCSD,3)');
figure
subplot(141)
imagesc(responsesPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(143)
imagesc(baselinePowAvg')
xlim([0 150])
caxis([-500 500])
subplot(144)
imagesc(saccadePowAvg')
xlim([0 150])
caxis([-500 500])
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_FM_evoked'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

[responsesFreq,responsesAudPowAvg] = PSD(mean(responsesAud,3)');
[responsesFreq,baselineAudPowAvg] = PSD(mean(baselineAud,3)');
[responsesFreq,saccadeAudPowAvg] = PSD(mean(saccadeAud,3)');
[responsesFreq,responsesAudCSDPowAvg] = PSD(mean(audioOnsetCSD,3)');
[responsesFreq,baselineAudCSDPowAvg] = PSD(mean(baselineAudCSD,3)');
[responsesFreq,saccadeAudCSDPowAvg] = PSD(mean(saccadeAudCSD,3)');
figure
subplot(141)
imagesc(responsesAudPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(143)
imagesc(baselineAudPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(144)
imagesc(saccadeAudPowAvg')
xlim([0 150])
caxis([-500 500])
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_AM_evoked'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

figure
subplot(141)
imagesc(responsesCSDPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(143)
imagesc(baselineCSDPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(144)
imagesc(saccadeCSDPowAvg')
xlim([0 150])
caxis([-500 500])
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_FM_CSD'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

figure
subplot(141)
imagesc(responsesAudCSDPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(143)
imagesc(baselineAudCSDPowAvg')
xlim([0 150])
caxis([-500 500])
subplot(144)
imagesc(saccadeAudCSDPowAvg')
xlim([0 150])
caxis([-500 500])
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_NormPower_AM_CSD'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')


% note that PSD for the different CSD data are almost identical! So only
% plotted for Ryan's CSD.
% [responsesFreq,responsesPowCSDAvg] = PSD(flashOnsetCSD');
% [responsesFreq,responsesPowCSDKacieAvg] = PSD(flashOnsetCSDKacie');
% imagesc(responsesPowCSDAvg')
% imagesc(responsesPowCSDKacieAvg')




% repeat for MUA activity
nUnits = numel(D.allMUAStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d MUAs...\n', nUnits);

psthWindow = [-0.025 0.225]; % secs before, secs after
kernelSigma = 0.01;
% nTime = fix(5*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
nTime = fix(10*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime); % starts at 0
t = psthT - psthWindow(1);
isFiringNonsparseByBlock = findNonsparseBlocks(D, D.allMUAStructs, R.vepmIndices);

psthResponse = nan(nChannels,size(t,2));
psthBootstrapErr = nan(nChannels, size(t,2));
psthBaseline = nan(nChannels,size(t,2));
psthBaselineBootstrapErr = nan(nChannels, size(t,2));
psthSaccade = nan(nChannels,size(t,2));
psthSaccadeBootstrapErr = nan(nChannels, size(t,2));
for j = 1:nChannels
    muaStruct = D.allMUAStructs{j};
    unitName = muaStruct.name;
    spikeTimes = muaStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));

    %% compute psth for each flash condition
    % this does not subtract out per-condition baseline activity

    % remove spike and event times during sparse blocks
    spikeTimesClean = spikeTimes;
    origFlashEventsClean = origFlashEvents;
    sparseBlocks = find(~isFiringNonsparseByBlock(:,j));
    rangeSparseRemovedFlashEvents = nan(numel(sparseBlocks), 2); % start, end
    for k = 1:numel(sparseBlocks)
        blockStartTime = D.blockStartTimes(R.vepmIndices(sparseBlocks(k)));
        blockStopTime = D.blockStopTimes(R.vepmIndices(sparseBlocks(k)));
        spikeTimesClean(spikeTimesClean >= blockStartTime & spikeTimesClean <= blockStopTime) = [];
        origFlashEventsClean(origFlashEventsClean >= blockStartTime & origFlashEventsClean <= blockStopTime) = [];
        rangeSparseRemovedFlashEvents(k,:) = [find(origFlashEvents >= blockStartTime, 1, 'first') find(origFlashEvents <= blockStopTime, 1, 'last')];
    end

%     nFlashes = numel(flashEventsClean);
%     nTrials = numel(preFlashEventsClean);
%     startIndicesStim = round((flashEventsClean + postFlashWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesStim = startIndicesStim + round(diff(postFlashWindowOffset) * Fs) - 1;
%     startIndicesBl = round((preFlashEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesBl = startIndicesBl + round(diff(baselineWindowOffset) * Fs) - 1;
%     startIndicesSac = round((preFlashEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesSac = startIndicesSac + round(diff(saccadeWindowOffset) * Fs) - 1;
%     t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
%     nTime = numel(t);
    allAlignedSpikeTimes = createdatamatpt(spikeTimes, origFlashEvents, psthWindow);
    allAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, origFlashEventsClean, psthWindow);
    baselineAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, preFlashEventsClean, abs(baselineWindowOffset));
    saccadeAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, preFlashEventsClean, [abs(saccadeWindowOffset(1)) saccadeWindowOffset(2)]);

    % compute psth response
    [psthResponse(j,:),~,psthBootstrapErr(j,:)] = fixedPsth(allAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    [psthBaseline(j,:),~,psthBaselineBootstrapErr(j,:)] = fixedPsth(baselineAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    [psthSaccade(j,:),~,psthSaccadeBootstrapErr(j,:)] = fixedPsth(saccadeAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    
end

% staggered channel line plot of average visually evoked MUA
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthResponse)));

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - MUA of Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Flashes_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for baseline period
xBounds = [-0.2 0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthBaseline)));

t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2);
% t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Flash event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_baseline_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for saccade period
xBounds = [-0.5 -0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthSaccade)));

t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Flash event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_saccade_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/boundedline/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/catuneven/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/Inpaint_nans/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/singlepatch/')
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthResponse)));
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2);
figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset,psthBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - MUA of Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Flashes_wBootstrapErr_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for baseline period
xBounds = [-0.2 0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthBaseline)));

t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset,psthBaselineBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Flash event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Baseline_wBootstrapErr_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for saccade period
xBounds = [-0.5 -0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthSaccade)));

t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset,psthSaccadeBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Flash event (s)');
title(sprintf('%s %s - MUA of saccade (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Saccade_wBootstrapErr_FM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% repeat MUA for auditory mapping

nUnits = numel(Daud.allMUAStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d MUAs...\n', nUnits);

psthWindow = [-0.025 0.225]; % secs before, secs after
kernelSigma = 0.01;
% nTime = fix(5*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
nTime = fix(10*sum(psthWindow)/kernelSigma) + 1; % number of time steps. 
psthT = linspace(0, sum(psthWindow), nTime); % starts at 0
t = psthT - psthWindow(1);
isFiringNonsparseByBlock = findNonsparseBlocks(Daud, Daud.allMUAStructs, Raud.vepmIndices);

psthResponse = nan(nChannels,size(t,2));
psthBootstrapErr = nan(nChannels, size(t,2));
psthBaseline = nan(nChannels,size(t,2));
psthBaselineBootstrapErr = nan(nChannels, size(t,2));
psthSaccade = nan(nChannels,size(t,2));
psthSaccadeBootstrapErr = nan(nChannels, size(t,2));
for j = 1:nChannels
    muaStruct = Daud.allMUAStructs{j};
    unitName = muaStruct.name;
    spikeTimes = muaStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
            nUnits, round(j/nUnits*100));

    %% compute psth for each flash condition
    % this does not subtract out per-condition baseline activity

    % remove spike and event times during sparse blocks
    spikeTimesClean = spikeTimes;
    origAudioEventsClean = origAudioEvents;
    sparseBlocks = find(~isFiringNonsparseByBlock(:,j));
    rangeSparseRemovedAudioEvents = nan(numel(sparseBlocks), 2); % start, end
    for k = 1:numel(sparseBlocks)
        blockStartTime = D.blockStartTimes(R.vepmIndices(sparseBlocks(k)));
        blockStopTime = D.blockStopTimes(R.vepmIndices(sparseBlocks(k)));
        spikeTimesClean(spikeTimesClean >= blockStartTime & spikeTimesClean <= blockStopTime) = [];
        origAudioEventsClean(origAudioEventsClean >= blockStartTime & origAudioEventsClean <= blockStopTime) = [];
        rangeSparseRemovedAudioEvents(k,:) = [find(origAudioEvents >= blockStartTime, 1, 'first') find(origAudioEvents <= blockStopTime, 1, 'last')];
    end

%     nFlashes = numel(flashEventsClean);
%     nTrials = numel(preFlashEventsClean);
%     startIndicesStim = round((flashEventsClean + postFlashWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesStim = startIndicesStim + round(diff(postFlashWindowOffset) * Fs) - 1;
%     startIndicesBl = round((preFlashEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesBl = startIndicesBl + round(diff(baselineWindowOffset) * Fs) - 1;
%     startIndicesSac = round((preFlashEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
%     endIndicesSac = startIndicesSac + round(diff(saccadeWindowOffset) * Fs) - 1;
%     t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
%     nTime = numel(t);
    allAlignedSpikeTimes = createdatamatpt(spikeTimes, origAudioEvents, psthWindow);
    allAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, origAudioEventsClean, psthWindow);
    baselineAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, preAudioEventsClean, abs(baselineWindowOffset));
    saccadeAlignedSpikeTimesClean = createdatamatpt(spikeTimesClean, preAudioEventsClean, [abs(saccadeWindowOffset(1)) saccadeWindowOffset(2)]);

    % compute psth response
    [psthResponse(j,:),~,psthBootstrapErr(j,:)] = fixedPsth(allAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    [psthBaseline(j,:),~,psthBaselineBootstrapErr(j,:)] = fixedPsth(baselineAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    [psthSaccade(j,:),~,psthSaccadeBootstrapErr(j,:)] = fixedPsth(saccadeAlignedSpikeTimesClean, kernelSigma, 2, psthT);
    
end

% staggered channel line plot of average visually evoked MUA
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthResponse)));

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Audio Onset (s)');
title(sprintf('%s %s - MUA of Audio events (N=%d) (%s)', sessionName, areaName, nFlashes));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Audio_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for baseline period
xBounds = [-0.2 0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthBaseline)));

t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2);
% t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2)-1/Fs;

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Audio event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_baseline_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for saccade period
xBounds = [-0.5 -0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthSaccade)));

t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Audio event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_saccade_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/boundedline/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/catuneven/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/Inpaint_nans/')
addpath('/Users/labmanager/Documents/MATLAB/kakearney-boundedline-pkg-50f7e4b/singlepatch/')
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthResponse)));
t = postAudioWindowOffset(1):1/Fs:postAudioWindowOffset(2);
figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset,psthBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - MUA of Audio events (N=%d) (%s)', sessionName, areaName, nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Audio_wBootstrapErr_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for baseline period
xBounds = [-0.2 0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthBaseline)));

t = baselineWindowOffset(1):1/Fs:baselineWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset,psthBaselineBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthBaseline(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Audio event (s)');
title(sprintf('%s %s - MUA of baseline (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Baseline_wBootstrapErr_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')

% repeat for saccade period
xBounds = [-0.5 -0.1];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(psthSaccade)));

t = saccadeWindowOffset(1):1/Fs:saccadeWindowOffset(2);

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    boundedline(t, responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset,psthSaccadeBootstrapErr(j,:),'alpha');
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*psthSaccade(j,:) - j*responsePlotYOffset) - 1]);
    end
end
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from pre-Audio event (s)');
title(sprintf('%s %s - MUA of saccade (N=%d) (%s)', sessionName, areaName, nTrials));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);
cd('/Users/labmanager/Documents/MATLAB/channelExploration/')
saveas(gcf,[sessionName '_MUA_Saccade_wBootstrapErr_AM'],'jpg')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/')









% eof