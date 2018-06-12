function [channelDataCARNorm,channelDataNorm,commonAverageNorm,isNoisyChannel] = preprocessLfps(...
        adjLfps, Fs, channelNames, processedDataDir, plotFileNamePrefix, hiCutoffFreq, rfMappingMode, v)
% remove events with outlier data and apply CAR and low-pass filtering to adjLfps

doOutlierCheckPlot = 1;

maxAbsNonOutlier = 0.5;
warningSDByChannelLarge = 2.5;
outlierMaxSDStep1 = 8; % first-pass threshold before LPF

isNanOrig = isnan(adjLfps); % track nans in original data
nChannels = size(adjLfps, 1);

assert(nChannels == 32);
isNoisyChannel = false(nChannels, 1);

%% do common average reference first
adjLfpsCAR = adjLfps - nanmean(adjLfps, 1);

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
plotFileName = sprintf('%s/%s-allFP-rfmMode%d-CARdata-SDByChannel-v%d.png', ...
        processedDataDir, plotFileNamePrefix, rfMappingMode, v);
fprintf('Saving SD by channel plot to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
close;

%% plot and mark channel data if SD abnormally high
for j = 1:nChannels
    if origSDByChannel(j) > upperThreshold
        warning('SD of Channel %d is abnormally high (%0.3f > %0.3f = %0.1f SDs)', j, ...
                origSDByChannel(j), upperThreshold, warningSDByChannelLarge);
        isNoisyChannel(j) = 1; % mark as noisy
        
        if doOutlierCheckPlot
            figure_tr_inch(6, 6);
            plot(adjLfpsCAR(j,:));
            title(sprintf('High Channel SD: %s', channelNames{j}), 'Interpreter', 'none');
            ylim([-0.5 0.5]);
            xlim([1 size(adjLfpsCAR, 2)]);
            
            plotFileName = sprintf('%s/%s-rfmMode%d-%s-CARdata-v%d.png', ...
                    processedDataDir, plotFileNamePrefix, rfMappingMode, channelNames{j}, v);
            fprintf('\tSaving channel %s data (CAR) plot to %s...\n', channelNames{j}, plotFileName);
            export_fig(plotFileName, '-nocrop');
            savefig([plotFileName(1:end-3) 'fig']);
            close;
        end
    end
end

%% remove obvious outliers and spike/lever artifacts before low-pass 
% filtering further. use CAR data for detection but clean up both CAR and
% non-CAR data
fprintf('Removing major outliers and spike/lever artifacts...\n');
for j = 1:nChannels
    fprintf('\tProcessing %s...\n', channelNames{j});
    
    extremeLogicalMaxAbs = abs(adjLfpsCAR(j,:)) > maxAbsNonOutlier;
    if any(extremeLogicalMaxAbs)
        fprintf('\t\tDetected %d outlier time points because voltage is too extreme (> %0.1f). Replacing them with mean.\n', ...
        sum(extremeLogicalMaxAbs), maxAbsNonOutlier);
    end
    % replace with mean for both CAR and non CAR data
    adjLfpCARCleanMaxAbs = adjLfpsCAR(j,:);
    adjLfpCARCleanMaxAbs(extremeLogicalMaxAbs) = nanmean(adjLfpCARCleanMaxAbs(~extremeLogicalMaxAbs));
    adjLfpCleanMaxAbs = adjLfps(j,:);
    adjLfpCleanMaxAbs(extremeLogicalMaxAbs) = nanmean(adjLfpCleanMaxAbs(~extremeLogicalMaxAbs));
    
    meanAdjLfpCAR = nanmean(adjLfpCARCleanMaxAbs);
    sdAdjLfpCAR = nanstd(adjLfpCARCleanMaxAbs);
    extremeLogicalSD = adjLfpCARCleanMaxAbs > meanAdjLfpCAR + outlierMaxSDStep1 * sdAdjLfpCAR | ...
            adjLfpCARCleanMaxAbs < meanAdjLfpCAR - outlierMaxSDStep1 * sdAdjLfpCAR;
    if any(extremeLogicalSD)
        fprintf('\t\tDetected %d outlier time points because voltage is too extreme (> %0.1f SDs). Replacing them with mean.\n', ...
                sum(extremeLogicalSD), outlierMaxSDStep1);
    end
    
    % TODO probably bad just to zero it out instead of smoothly dealing
    % with it or NaNing it. fix this with filtering or case by case
    % replace outlier time points with the average for this channel
    adjLfps(j,extremeLogicalMaxAbs | extremeLogicalSD) = nanmean(adjLfpCleanMaxAbs); 
    clear adjLfpsCleanMaxAbs;
end

clear adjLfpsCAR; % discard CAR data -- recompute after filtering

%% low-pass filter data even further
% SLOW - may be better to concatenate chunks of relevant data.....
fprintf('Low-pass filtering all channels, all data at %d Hz... this will take some time... \n', hiCutoffFreq); 
bFirLowPass = fir1(3*fix(Fs/hiCutoffFreq), hiCutoffFreq/(Fs/2), 'low');

adjLfps(isnan(adjLfps)) = 0; % zero out nans
adjLfpsLP = filtfilt(bFirLowPass, 1, adjLfps')';
adjLfpsLP(isNanOrig) = NaN; % set missing vals back to nan
clear adjLfps;

% TODO better to CAR first, then filter? or filter then CAR?

%% redo common average reference
commonAverageLP = nanmean(adjLfpsLP, 1);
adjLfpsLP_CAR = adjLfpsLP - commonAverageLP;

%% z-score normalize each channel
% units are now standard deviations away from the mean
channelDataCARNorm = (adjLfpsLP_CAR - nanmean(adjLfpsLP_CAR, 2)) ./ nanstd(adjLfpsLP_CAR, 0, 2);
assert(all(size(adjLfpsLP_CAR) == size(channelDataCARNorm)));

channelDataNorm = (adjLfpsLP - nanmean(adjLfpsLP, 2)) ./ nanstd(adjLfpsLP, 0, 2);
assert(all(size(adjLfpsLP) == size(channelDataNorm)));

% do the same to refMeanLP
commonAverageNorm = (commonAverageLP - nanmean(commonAverageLP, 2)) ./ nanstd(commonAverageLP, 0, 2);
assert(all(size(commonAverageLP) == size(commonAverageNorm)));

fprintf('Done preprocessing LFPs.\n');