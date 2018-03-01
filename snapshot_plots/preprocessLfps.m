function [channelDataNorm,eventOnsetClean,isEventOutlier,isNoisyChannel] = preprocessLfps(adjLfps, Fs, channelNames, eventOnset, ...
        processedDataDir, blockName, rfMappingNewMode, v)
% remove events with outlier data and apply CAR and low-pass filtering to adjLfps

doOutlierCheckPlot = 1;

maxAbsNonOutlier = 0.5;
warningSDByChannelLarge = 2.5;
outlierMaxSDStep1 = 6; % first-pass threshold before LPF
outlierMaxSDStep2 = 4.5; % second-pass threshold after LPF
outlierCheckWindowOffset = [-0.25 0.3]; % seconds around event

isNanOrig = isnan(adjLfps); % track nans in original data
nEvents = numel(eventOnset);
nChannels = size(adjLfps, 1);
%assert(nChannels == 32);
isNoisyChannel = false(nChannels, 1);


%% do common average reference first
adjLfpsCAR = adjLfps - nanmean(adjLfps, 1);
clear adjLfps;

%% plot histogram of SD of each channel's LFP
origSDByChannel = nanstd(adjLfpsCAR, 0, 2);
upperThreshold = mean(origSDByChannel) + warningSDByChannelLarge * std(origSDByChannel);

figure_tr_inch(6, 6);
bar(1:nChannels, origSDByChannel);
hold on;
plot([1 nChannels], upperThreshold * [1 1], 'm--');
text(1, upperThreshold + 0.001, ...
        sprintf('Warning threshold (%0.1f SDs = %0.3f)', warningSDByChannelLarge, ...
        upperThreshold), 'Color', 'm', ...
        'VerticalAlignment', 'bottom');
title('LFPs SD by Channel');
xlim([1 nChannels] + [-0.5 0.5]);
xlabel('Channel Number');
ylabel('Standard Deviation of LFP (mV)');
plotFileName = sprintf('%s/allFP-%s-%s-%s-rfm_mode%d-CARdata-SDByChannel_v%d.png', ...
        processedDataDir, channelNames{[1 end]}, blockName, rfMappingNewMode, v);
fprintf('Saving SD by channel plot to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;

%% print and sanity check
for j = 1:nChannels
    if origSDByChannel(j) > upperThreshold
        warning('SD of Channel %d is abnormally high (%0.2f > %0.2f)', j, ...
                origSDByChannel(j), upperThreshold);
        isNoisyChannel(j) = 1;
        if doOutlierCheckPlot
            figure_tr_inch(6, 6);
            plot(adjLfpsCAR(j,:));
            title(channelNames{j}, 'Interpreter', 'none');
            ylim([-0.5 0.5]);
            xlim([1 size(adjLfpsCAR, 2)]);
            
            plotFileName = sprintf('%s/%s-%s-rfm_mode%d-CARdata_v%d.png', ...
                    processedDataDir, channelNames{j}, blockName, rfMappingNewMode, v);
            fprintf('\tSaving CAR data plot to %s...\n', plotFileName);
            export_fig(plotFileName, '-nocrop');
            savefig([plotFileName(1:end-3) 'fig']);
            close;
        end
    end
end

%% remove obvious outliers and spike/lever artifacts before low-pass 
% filtering further
fprintf('Removing major outliers and spike/lever artifacts...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
    
    extremeLogicalMaxAbs = abs(adjLfpsCAR(j,:)) > maxAbsNonOutlier;
    if any(extremeLogicalMaxAbs)
        fprintf('\tDetected %d outlier time points because voltage in trial is too extreme (> %0.1f). Zeroing them out.\n', ...
        sum(extremeLogicalMaxAbs), maxAbsNonOutlier);
    end
    % replace with mean
    adjLfpsCleanMaxAbs = adjLfpsCAR(j,:);
    adjLfpsCleanMaxAbs(extremeLogicalMaxAbs) = nanmean(adjLfpsCleanMaxAbs(~extremeLogicalMaxAbs));
    
    meanAdjLfp = nanmean(adjLfpsCleanMaxAbs);
    sdAdjLfp = nanstd(adjLfpsCleanMaxAbs);
    extremeLogicalSD = adjLfpsCleanMaxAbs > meanAdjLfp + outlierMaxSDStep1 * sdAdjLfp | ...
            adjLfpsCleanMaxAbs < meanAdjLfp - outlierMaxSDStep1 * sdAdjLfp;
    if any(extremeLogicalSD)
        fprintf('\tDetected %d outlier time points because voltage in trial is too extreme (> %0.1f SDs). Zeroing them out.\n', ...
                sum(extremeLogicalSD), outlierMaxSDStep1);
    end
    % TODO probably bad just to zero it out instead of smoothly dealing
    % with it or NaNing it. fix this with filtering or case by case
    % replace outlier time points with the average for this channel
    adjLfpsCAR(j,extremeLogicalMaxAbs | extremeLogicalSD) = meanAdjLfp; 
    clear adjLfpsCleanMaxAbs;
end

%% low-pass filter data even further
% ideally don't low-pass at 300 Hz and then again at another freq
% use FIR1 filter
fprintf('Low-pass filtering...\n');
hiCutoffFreqCustom = 100; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
adjLfpsCAR(isnan(adjLfpsCAR)) = 0; % zero out nans
adjLfpsLP = filtfilt(bFirLowPassCustom, 1, adjLfpsCAR')';
clear adjLfpsCAR;

% set missing vals back to nan
adjLfpsLP(isNanOrig) = NaN;

%% redo common average reference
adjLfpsLpCAR = adjLfpsLP - nanmean(adjLfpsLP, 1);
clear adjLfpsLP;

%% normalize
% normalize each channel by its mean and standard deviation
% units are now standard deviations away from the mean
for j = 1:nChannels
    fprintf('Channel %d: Mean: %0.3f, SD: %0.2f\n', j, nanmean(adjLfpsLpCAR(j,:)), nanstd(adjLfpsLpCAR(j,:)));
end
channelDataNorm = (adjLfpsLpCAR - nanmean(adjLfpsLpCAR, 2)) ./ nanstd(adjLfpsLpCAR, 0, 2);
assert(all(size(adjLfpsLpCAR) == size(channelDataNorm)));
clear adjLfpsLpCAR;

%% detect events with outlier activity (e.g. amp saturated)
outlierCheckStartIndices = round((eventOnset + outlierCheckWindowOffset(1)) * Fs); 
outlierCheckEndIndices = outlierCheckStartIndices + round(diff(outlierCheckWindowOffset) * Fs) - 1;
outlierCheckT = outlierCheckWindowOffset(1):1/Fs:outlierCheckWindowOffset(2)-1/Fs;
nOutlierCheckTime = numel(outlierCheckT);
assert(all(nOutlierCheckTime == (outlierCheckEndIndices - outlierCheckStartIndices + 1)));

isEventOutlier = false(nEvents, 1);
fprintf('Checking for events with outlier activity (e.g. amp saturated)...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
    
    isEventOutlierThisChannel = false(nEvents, 1);
    for n = 1:nEvents
        dataInOutlierCheckPeriod = channelDataNorm(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
        % mark outlier flashes
        if any(isnan(dataInOutlierCheckPeriod))
            isEventOutlierThisChannel(n) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because of NaN in trial\n', ...
                    n, sum(isnan(dataInOutlierCheckPeriod)));
        elseif any(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2)
            isEventOutlierThisChannel(n) = 1;
            fprintf('\t\tFlash %d: Detected %d outlier time points because voltage in trial is too extreme (> %0.1f SDs)\n', ...
                    n, sum(abs(dataInOutlierCheckPeriod) > outlierMaxSDStep2), outlierMaxSDStep2);
        end
    end
    
    if doOutlierCheckPlot && any(isEventOutlierThisChannel)
        figure_tr_inch(6, 6);
        suptitle(sprintf('Outlier Check on Channel %s', channelNames{j}));
        hold on;
        for n = 1:nEvents
            dataInOutlierCheckPeriod = channelDataNorm(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
            plot(outlierCheckT, dataInOutlierCheckPeriod);
        end
        for n = 1:nEvents
            if isEventOutlierThisChannel(n)
                dataInOutlierCheckPeriod = channelDataNorm(j,outlierCheckStartIndices(n):outlierCheckEndIndices(n));
                % plot over other lines, with thicker line
                plot(outlierCheckT, dataInOutlierCheckPeriod, 'LineWidth', 3);
            end
        end
        
        title('Original Responses');
        xlabel('Time from Event Onset');
        ylabel('SDs from Overall Mean');
        plot([0 0], [-1000 1000], '-', 'Color', [0.5 0.5 0.5]);
        xlim(outlierCheckWindowOffset);
        ylim([-8 8]);
        
        plotFileName = sprintf('%s/%s-%s-rfm_mode%d-outlierCheck_v%d.png', ...
                processedDataDir, channelNames{j}, blockName, rfMappingNewMode, v);
        fprintf('\t\tSaving outlier check plot to %s...\n', plotFileName);
        export_fig(plotFileName, '-nocrop');
        close;
    end
    
    isEventOutlier = isEventOutlier | isEventOutlierThisChannel;
end

% remove outlier flashes
fprintf('Detected %d outlier events out of %d (%d%%).\n', sum(isEventOutlier), ...
        nEvents, round(sum(isEventOutlier)/nEvents * 100));
eventOnsetClean = eventOnset(~isEventOutlier);
if sum(isEventOutlier) > 0.1 * nEvents
    error('More than 10%% of flashes have outlier responses.');
end
