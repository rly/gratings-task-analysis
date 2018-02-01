

%%
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
sessionName = 'M20170127';
besselPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-bessel.pl2';
butterPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-butterworth.pl2';
bessel1HzPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-bessel1Hz.pl2';

besselFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, besselPl2FileName);
butterFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, butterPl2FileName);
bessel1HzFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, bessel1HzPl2FileName);

% assume besselFile and butterFile have same data

%% load basic info
besselDataInfo = PL2GetFileIndex(besselFileName);
butterDataInfo = PL2GetFileIndex(butterFileName);
bessel1HzDataInfo = PL2GetFileIndex(bessel1HzFileName);

dataInfo = besselDataInfo;
Fs = dataInfo.TimestampFrequency;

startTime = 1; % from 1 second in
endTime = 11; % to 11 seconds in
maxTimeToSave = (endTime - startTime) * Fs;

doesWBChannelHaveData = false(numel(dataInfo.AnalogChannels), 1);
fprintf('WB channels with data: ');
for i = 1:numel(dataInfo.AnalogChannels)
    if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'WB') && ...
            dataInfo.AnalogChannels{i}.NumValues > 0 && ...
            dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
        adInfo = PL2AdTimeSpan(besselFileName, dataInfo.AnalogChannels{i}.Name, 0, 0.01);
        if any(adInfo.Values)
            fprintf('%d ', i);
            doesWBChannelHaveData(i) = 1;
        end
    end
end
fprintf('\n');


doesSPKCChannelHaveData = false(numel(dataInfo.AnalogChannels), 1);
fprintf('SPKC channels with data: ');
for i = 1:numel(dataInfo.AnalogChannels)
    if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'SPKC') && ...
            dataInfo.AnalogChannels{i}.NumValues > 0 && ...
            dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
        adInfo = PL2AdTimeSpan(besselFileName, dataInfo.AnalogChannels{i}.Name, 0, 0.01);
        if any(adInfo.Values)
            fprintf('%d ', i);
            doesSPKCChannelHaveData(i) = 1;
        end
    end
end
fprintf('\n');

doesFiltWBChannelHaveData = false(numel(dataInfo.AnalogChannels), 1);
fprintf('SPKC channels with data: ');
for i = 1:numel(dataInfo.AnalogChannels)
    if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'FILT_WB') && ...
            dataInfo.AnalogChannels{i}.NumValues > 0 && ...
            dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
        adInfo = PL2AdTimeSpan(besselFileName, dataInfo.AnalogChannels{i}.Name, 0, 0.01);
        if any(adInfo.Values)
            fprintf('%d ', i);
            doesFiltWBChannelHaveData(i) = 1;
        end
    end
end
fprintf('\n');

wbChannelsToRun = find(doesWBChannelHaveData);
spkcChannelsToRun = find(doesSPKCChannelHaveData);
filtWBChannelsToRun = find(doesFiltWBChannelHaveData);
assert(all([numel(wbChannelsToRun) numel(filtWBChannelsToRun)] == numel(spkcChannelsToRun)));

%%
% 4-pole butterworth bandpass filter
[butterBandPassB, butterBandPassA] = butter(4, [300 8000] / (Fs/2));

% 4-pole butterworth high-pass filter
[butterHighPassB, butterHighPassA] = butter(4, 300/(Fs/2), 'high');

%%
for i = 3%1:numel(wbChannelsToRun)
    channelID = wbChannelsToRun(i);
    fprintf('\nLoading data from channel %d: t = %0.1f s to t = %0.1f s...\n', channelID, startTime, endTime);
    wbInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{channelID}.Name, startTime, endTime);
    spkcInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{spkcChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
    plexonBesselHighPassFiltWBInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
    plexonButterHighPassFiltWBInfo = PL2AdTimeSpan(butterFileName, butterDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
    plexonBessel1HzHighPassFiltWBInfo = PL2AdTimeSpan(bessel1HzFileName, bessel1HzDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering

    rawData = wbInfo.Values;
    plexonBesselHighPassData = plexonBesselHighPassFiltWBInfo.Values;
    plexonButterHighPassData = plexonButterHighPassFiltWBInfo.Values;
    plexonBessel1HzHighPassData = plexonBessel1HzHighPassFiltWBInfo.Values;
    matlabButterHighPassNLPData = filter(butterHighPassB, butterHighPassA, rawData); % non-linear phase
    matlabButterBandPassNLPData = filter(butterBandPassB, butterBandPassA, rawData);
    matlabButterHighPassZPData = filtfilt(butterHighPassB, butterHighPassA, rawData); % zero phase
    matlabButterBandPassZPData = filtfilt(butterBandPassB, butterBandPassA, rawData);
    matlabWaveletHighPassData =  wavefilter(rawData', 6);
    
    plexonBesselHighPassData1 = spkcInfo.Values;
    assert(max(abs(plexonBesselHighPassData1 - plexonBesselHighPassData)) < 1e-5);
    clear plexonBesselData1;
        
%     figure_tr_inch(18, 10);
%     subaxis(2, 1, 1);
%     hold on;
%     plot(plexonBesselHighPassData);
%     title(sprintf('Channel %d, Bessel', channelID));
%     subaxis(2, 1, 2);
%     plot(rawData);
%     title(sprintf('Channel %d, Raw', channelID));
    
    
    %% use filtered raw data
    
    % extract waveforms -- lower threshold only for now
    [extractedWaveforms,startWaveformInds] = findSpikeWaveformsUsingThreshold(plexonBesselHighPassData, 4, 16, 40, 40, 0, channelID);
    [extractedWaveforms,startWaveformInds] = findSpikeWaveformsUsingThreshold(plexonBesselHighPassData, 4, 80, 160, 40, 0, channelID);
    nWaveforms = size(extractedWaveforms, 1);
    nWaveformTime = size(extractedWaveforms, 2);
    
    % do pca
    [pcaCoeff,pcaScore,~,~,pcaPctExplained] = pca(extractedWaveforms);
    
    fprintf('PCA: %d variables, %d observations\n', nWaveforms, nWaveformTime);
    fprintf('\tPC1 explains %0.1f%% of the variance.\n', pcaPctExplained(1));
    fprintf('\tPC1 + PC2 explain %0.1f%% of the variance.\n', pcaPctExplained(1) + pcaPctExplained(2));
    
%     figure;
%     hold on;
%     plot(pcaScore(:,1), pcaScore(:,2), '.');
%     title(sprintf('Channel %d, PCA on Lower-Threshold Waveforms', wbChannelsToRun(i)));

    nWaveformSamples = 241;
    reExtractedWaveformInds = cell2mat(arrayfun(@(x) x:x+nWaveformSamples-1, startWaveformInds, 'UniformOutput', false));
    
    extractedWaveformsRaw = rawData(reExtractedWaveformInds);
    extractedWaveformsPlexonBesselHighPass = plexonBesselHighPassData(reExtractedWaveformInds);
    extractedWaveformsPlexonButterHighPass = plexonButterHighPassData(reExtractedWaveformInds);
    extractedWaveformsPlexonBessel1HzHighPass = plexonBessel1HzHighPassData(reExtractedWaveformInds);
    extractedWaveformsMatlabButterHighPassNLPData = matlabButterHighPassNLPData(reExtractedWaveformInds);
    extractedWaveformsMatlabButterBandPassNLPData = matlabButterBandPassNLPData(reExtractedWaveformInds);
    extractedWaveformsMatlabButterHighPassZPData = matlabButterHighPassZPData(reExtractedWaveformInds);
    extractedWaveformsMatlabButterBandPassZPData = matlabButterBandPassZPData(reExtractedWaveformInds);
    extractedWaveformsMatlabWaveletHighPassData =  matlabWaveletHighPassData(reExtractedWaveformInds);
%     load('sampleWaveformsAll_ch5.mat', 'note');
%     save('sampleWaveformsAll_ch5.mat', 'extractedWaveforms*', 'note');
    
%     figure_tr_inch(12, 5);
%     subaxis(1, 2, 1);
%     hold on;
%     plot(filtData(reExtractedWaveformInds)');
%     title('Filtered Waveforms');
%     
%     subaxis(1, 2, 2);
%     hold on;
%     plot(wbInfo.Values(reExtractedWaveformInds)');
%     title('Raw Waveforms');
    
    % separate the two groups (ch5) by the post-trough peak at t=35 (bessel
    % only)
    isGroup1 = extractedWaveforms(:,35) > 0.03;
    isGroup2 = ~isGroup1;
    reExtractedWaveformIndsGroup1 = cell2mat(arrayfun(@(x) x:x+nWaveformSamples, startWaveformInds(isGroup1), 'UniformOutput', false));
    reExtractedWaveformIndsGroup2 = cell2mat(arrayfun(@(x) x:x+nWaveformSamples, startWaveformInds(isGroup2), 'UniformOutput', false));
    
    figure_tr_inch(12, 5);
    subaxis(1, 2, 1);
    hold on;
    plot(plexonBesselHighPassData(reExtractedWaveformIndsGroup1)', 'b');
    plot(plexonBesselHighPassData(reExtractedWaveformIndsGroup2)', 'r');
    plot(mean(plexonBesselHighPassData(reExtractedWaveformIndsGroup1),1), 'k', 'LineWidth', 3);
    plot(mean(plexonBesselHighPassData(reExtractedWaveformIndsGroup2),1), 'k', 'LineWidth', 3);
    title('Filtered Waveforms');
    
    subaxis(1, 2, 2);
    hold on;
    plot(wbInfo.Values(reExtractedWaveformIndsGroup1)', 'b');
    plot(wbInfo.Values(reExtractedWaveformIndsGroup2)', 'r');
    plot(mean(wbInfo.Values(reExtractedWaveformIndsGroup1),1), 'k', 'LineWidth', 3);
    plot(mean(wbInfo.Values(reExtractedWaveformIndsGroup2),1), 'k', 'LineWidth', 3);
    title('Raw Waveforms');
    
    %% try t-sne
%     if nWaveforms >= 50
%         % do t-sne, reduce to 2 dimensions
%         fprintf('t-SNE: %d variables, %d observations\n', nWaveforms, nWaveformTime);
%         tsneVals = tsne(extractedWaveforms);
% 
%         figure;
%         hold on;
%         plot(tsneVals(:,1), tsneVals(:,2), '.');
%         title(sprintf('t-SNE on Waveforms: Channel %d', channels(i)));
%     end
    
    %% try motif discovery
%     data = spkcInfo.Values;
%     initNumSDsThresh = 4; % both + and -
%     initPreExtremeSamples = 10;
%     initPostExtremeSamples = 30;
%     isUseMAD = 0;
%     isExtreme = findDataCrossingThreshold(...
%             data, initNumSDsThresh, initPreExtremeSamples, initPostExtremeSamples, isUseMAD);
%     fprintf('Considering %d/%d (%d%%) extreme data points.\n', sum(isExtreme), numel(data), round(sum(isExtreme)/numel(data)*100));
% 
%     % concatenate extreme segments, pulling from unfiltered data
%     extremeData = spkcInfo.Values(isExtreme);
%     
%     tic; 
%     [matrixProfile, profileIndex, motifIdxs] = interactiveMatrixProfileVer2(extremeData, 60); 
%     toc

    
end

stop

%% simulate a waveform, filter it
channelID = 9;
load(sprintf('sampleWaveformsAll_ch%d.mat', channelID));
Fs = 40000;

for waveformInd = 1:25%1:size(extractedWaveforms, 1)
    rawWaveform = extractedWaveformsRaw(waveformInd,:);
    rawWaveform = rawWaveform - mean(rawWaveform);
    
    plexonBesselHighPassWaveform = extractedWaveformsPlexonBesselHighPass(waveformInd,:);
    plexonBesselHighPassWaveform = plexonBesselHighPassWaveform - mean(plexonBesselHighPassWaveform);
    
    plexonButterHighPassWaveform = extractedWaveformsPlexonButterHighPass(waveformInd,:);
    plexonButterHighPassWaveform = plexonButterHighPassWaveform - mean(plexonButterHighPassWaveform);

    plexonBessel1HzHighPassWaveform = extractedWaveformsPlexonBessel1HzHighPass(waveformInd,:);
    plexonBessel1HzHighPassWaveform = plexonBessel1HzHighPassWaveform - mean(plexonBessel1HzHighPassWaveform);
    
    matlabButterHighPassNLPWaveform = extractedWaveformsMatlabButterHighPassNLPData(waveformInd,:);
    matlabButterHighPassNLPWaveform = matlabButterHighPassNLPWaveform - mean(matlabButterHighPassNLPWaveform);
    
    matlabButterBandPassNLPWaveform = extractedWaveformsMatlabButterBandPassNLPData(waveformInd,:);
    matlabButterBandPassNLPWaveform = matlabButterBandPassNLPWaveform - mean(matlabButterBandPassNLPWaveform);
    
    matlabButterHighPassZPWaveform = extractedWaveformsMatlabButterHighPassZPData(waveformInd,:);
    matlabButterHighPassZPWaveform = matlabButterHighPassZPWaveform - mean(matlabButterHighPassZPWaveform);
    
    matlabButterBandPassZPWaveform = extractedWaveformsMatlabButterBandPassZPData(waveformInd,:);
    matlabButterBandPassZPWaveform = matlabButterBandPassZPWaveform - mean(matlabButterBandPassZPWaveform);
    
    matlabWaveletHighPassWaveform = extractedWaveformsMatlabWaveletHighPassData(waveformInd,:);
    matlabWaveletHighPassWaveform = matlabWaveletHighPassWaveform - mean(matlabWaveletHighPassWaveform);
    
    t = (1:numel(rawWaveform)) / Fs * 1000;
    cols = linspecer(7);
    
    figure_tr_inch(13, 8);
    subaxis(1, 1, 1, 'ML', 0.07, 'MB', 0.09);
    hold on;
    plot(t, plexonBesselHighPassWaveform, 'LineWidth', 2, 'Color', cols(1,:));
    plot(t, plexonButterHighPassWaveform, 'LineWidth', 2, 'Color', cols(2,:));
%     plot(t, plexonBessel1HzHighPassWaveform, 'LineWidth', 2, 'Color', cols(3,:));
    plot(t, matlabButterHighPassNLPWaveform, ':', 'LineWidth', 2, 'Color', cols(3,:));
    plot(t, matlabButterBandPassNLPWaveform, 'LineWidth', 2, 'Color', cols(4,:));
    plot(t, matlabButterHighPassZPWaveform, 'LineWidth', 2, 'Color', cols(5,:));
    plot(t, matlabButterBandPassZPWaveform, 'LineWidth', 2, 'Color', cols(6,:));
    plot(t, matlabWaveletHighPassWaveform, 'LineWidth', 2, 'Color', cols(7,:));
    plot(t, rawWaveform, 'k', 'LineWidth', 2);
    
    box off;
    xlim([0 t(end)]);
    ylim([-0.1 0.08]);
    xlabel('Time (ms)');
    ylabel('Voltage (V)');
    set(gca, 'FontSize', 14);
    title(sprintf('Waveform #%d', waveformInd), 'FontSize', 16);

    legend({'Plexon Bessel High-Pass 300 Hz', ...
            'Plexon Butterworth High-Pass 300 Hz', ...
            ... %'Plexon Bessel High-Pass 1 Hz', ...
            'Matlab Butterworth High-Pass 300 Hz NLP', ...
            'Matlab Butterworth Band-Pass 300-8000 Hz NLP', ...
            'Matlab Butterworth High-Pass 300 Hz ZP', ...
            'Matlab Butterworth Band-Pass 300-8000 Hz ZP', ...
            'Matlab Wavelet High-Pass 312.5 Hz', ...
            'Raw', ...
            }, 'Location', 'SouthEast');

    plotFileName = sprintf('filterComparison_ch%d_waveform%d.png', channelID, waveformInd);
    export_fig(plotFileName, '-nocrop');
    close;
    
    % TODO next: run bessel and butterworth filters using OFS on the WB
    % data and compare...
    
    % almost no difference (< 1e-5) between OmniPlex SPKC and OFS WB -> Bessel High-Pass 
    
    % is it because plexon does only forward filtering and I am doing both
    % forward and backward filtering in matlab?
    
    % cannot implement zero-phase filtering online
end