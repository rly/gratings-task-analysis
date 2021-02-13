

%%
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
sessionName = 'M20170127';
besselPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-bessel.pl2';
% butterPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-butterworth.pl2';
% bessel1HzPl2FileName = '20170127-g1_wb1_wb5_wb9_wb25_spkc1_spkc5_spkc9_spkc25-bessel1Hz.pl2';

besselFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, besselPl2FileName);
% butterFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, butterPl2FileName);
% bessel1HzFileName = sprintf('%s/%s/%s', dataDirRoot, sessionName, bessel1HzPl2FileName);

% assume besselFile and butterFile have same data

%% load basic info
besselDataInfo = PL2GetFileIndex(besselFileName);
% butterDataInfo = PL2GetFileIndex(butterFileName);
% bessel1HzDataInfo = PL2GetFileIndex(bessel1HzFileName);

dataInfo = besselDataInfo;
Fs = dataInfo.TimestampFrequency;

startTime = 1; % from 1 second in
endTime = 6; % to 6 seconds in
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
% [butterBandPassB, butterBandPassA] = butter(4, [300 8000] / (Fs/2));

% 4-pole butterworth high-pass filter
% [butterHighPassB, butterHighPassA] = butter(4, 300/(Fs/2), 'high');

hiCutoffFreqCustom = 300; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');

%%
for i = 2%1:numel(wbChannelsToRun)
    channelID = wbChannelsToRun(i);
    fprintf('\nLoading data from channel %d: t = %0.1f s to t = %0.1f s...\n', channelID, startTime, endTime);
    wbInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{channelID}.Name, startTime, endTime);
    spkcInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{spkcChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
    plexonBesselHighPassFiltWBInfo = PL2AdTimeSpan(besselFileName, besselDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
%     plexonButterHighPassFiltWBInfo = PL2AdTimeSpan(butterFileName, butterDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering
%     plexonBessel1HzHighPassFiltWBInfo = PL2AdTimeSpan(bessel1HzFileName, bessel1HzDataInfo.AnalogChannels{filtWBChannelsToRun(i)}.Name, startTime, endTime); % assume same ordering

    rawData = wbInfo.Values;
    plexonBesselHighPassData = plexonBesselHighPassFiltWBInfo.Values;
%     plexonButterHighPassData = plexonButterHighPassFiltWBInfo.Values;
%     plexonBessel1HzHighPassData = plexonBessel1HzHighPassFiltWBInfo.Values;
%     matlabButterHighPassNLPData = filter(butterHighPassB, butterHighPassA, rawData); % non-linear phase
%     matlabButterBandPassNLPData = filter(butterBandPassB, butterBandPassA, rawData);
%     matlabButterHighPassZPData = filtfilt(butterHighPassB, butterHighPassA, rawData); % zero phase
%     matlabButterBandPassZPData = filtfilt(butterBandPassB, butterBandPassA, rawData);
%     matlabWaveletHighPassData =  wavefilter(rawData', 6);
    
    plexonBesselHighPassData1 = spkcInfo.Values;
    assert(max(abs(plexonBesselHighPassData1 - plexonBesselHighPassData)) < 1e-5);
    clear plexonBesselData1;
    
    matlabLowPassData = filtfilt(bFirLowPassCustom, 1, rawData);
    
    t = 0:1/Fs:(endTime-startTime);
    
    figure_tr_inch(4, 2.5);
    subaxis(1, 1, 1, 'ML', 0.05, 'MB', 0.05, 'MT', 0.05, 'MR', 0.05);
    plot(t, rawData);
    ylim([min(rawData) max(rawData)]);
    box off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    
    plotFileName = sprintf('ch%d_sample_raw.png', channelID);
    export_fig(plotFileName, '-nocrop');
    
    figure_tr_inch(4, 2.5);
    subaxis(1, 1, 1, 'ML', 0.05, 'MB', 0.05, 'MT', 0.05, 'MR', 0.05);
    plot(t, plexonBesselHighPassData);
    ylim([min(plexonBesselHighPassData) max(plexonBesselHighPassData)]);
    box off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    
    plotFileName = sprintf('ch%d_sample_highpass.png', channelID);
    export_fig(plotFileName, '-nocrop');
    
    figure_tr_inch(4, 2.5);
    subaxis(1, 1, 1, 'ML', 0.05, 'MB', 0.05, 'MT', 0.05, 'MR', 0.05);
    plot(t, matlabLowPassData);
    ylim([min(matlabLowPassData) max(matlabLowPassData)]);
    box off;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    
    plotFileName = sprintf('ch%d_sample_lowpass.png', channelID);
    export_fig(plotFileName, '-nocrop');
end