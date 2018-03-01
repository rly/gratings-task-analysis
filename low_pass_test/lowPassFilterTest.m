% quick script to compare Plexon low-pass-filtered WB data (FP) with my
% low-pass-filtered WB data at 300 Hz and 100 Hz

% the goal was to try to figure out why channels 7, 10, and 13 of M20170608
% have much higher variance in the LFP. perhaps there is extra noise or
% spikes are contaminating the Plexon filter

%%
clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
sessionName = 'M20170608';
wbPL2FileName = 'd1_34.50mm_d2_37.00mm_rfm1.pl2';
wbPL2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, wbPL2FileName);

channelInd = 30; % indices

params.tapers = [3 5];
params.pad = 0;
params.Fs = 1000;
params.fpass = [2 100];
params.err = [2 0.0500];
params.trialave = 1;

startTime = 10; % seconds
endTime = 210; % seconds 

%% load basic info
dataInfo = PL2GetFileIndex(wbPL2FilePath);
Fs = dataInfo.TimestampFrequency;

doesWBChannelHaveData = false(numel(dataInfo.AnalogChannels), 1);
fprintf('WB channels with data: ');
for i = 1:numel(dataInfo.AnalogChannels)
    if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'WB') && ...
            dataInfo.AnalogChannels{i}.NumValues > 0 && ...
            dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
        adInfo = PL2AdTimeSpan(wbPL2FilePath, dataInfo.AnalogChannels{i}.Name, 0, 0.01);
        if any(adInfo.Values)
            fprintf('%d ', i);
            doesWBChannelHaveData(i) = 1;
        end
    end
end
fprintf('\n');

doesFPChannelHaveData = false(numel(dataInfo.AnalogChannels), 1);
fprintf('FP channels with data: ');
for i = 1:numel(dataInfo.AnalogChannels)
    if dataInfo.AnalogChannels{i}.Enabled && strcmp(dataInfo.AnalogChannels{i}.SourceName, 'FP') && ...
            dataInfo.AnalogChannels{i}.NumValues > 0 && ...
            dataInfo.AnalogChannels{i}.Channel <= 125 % SKIP extra analog signal data
        adInfo = PL2AdTimeSpan(wbPL2FilePath, dataInfo.AnalogChannels{i}.Name, 0, 0.01);
        if any(adInfo.Values)
            fprintf('%d ', i);
            doesFPChannelHaveData(i) = 1;
        end
    end
end
fprintf('\n');

wbChannelsToRun = find(doesWBChannelHaveData);
fpChannelsToRun = find(doesFPChannelHaveData);

%% load data
maxTimeToSave = (endTime - startTime) * Fs;
fprintf('Reading channel %d from %0.1f to %0.1f seconds.\n', channelInd, startTime, endTime);

wbChannelName = sprintf('WB%03d', channelInd);
adInfo = PL2AdTimeSpan(wbPL2FilePath, wbChannelName, startTime, endTime);
wbData = adInfo.Values;
wbFs = adInfo.ADFreq;
tWB = startTime * wbFs : endTime * wbFs;
clear adInfo;

fpChannelName = sprintf('FP%03d', channelInd);
adInfo = PL2AdTimeSpan(wbPL2FilePath, fpChannelName, startTime, endTime);
lfpData = adInfo.Values;
lfpFs = adInfo.ADFreq;
tLFP = startTime * lfpFs : endTime * lfpFs;
if numel(lfpData) == numel(tLFP) - 1
    lfpData(end + 1) = 0; % pad with zero
end
clear adInfo;

%% low pass filter at 300 Hz using FIR1 filter
fprintf('Low-pass filtering...\n');
hiCutoffFreqCustom = 300; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
wbDataLPF300 = filtfilt(bFirLowPassCustom, 1, wbData')';
wbDataLPF300DS = downsample(wbDataLPF300, round(wbFs / lfpFs));

%% low pass filter at 100 Hz using FIR1 filter
fprintf('Low-pass filtering...\n');
hiCutoffFreqCustom = 100; % low-pass filter at 100 Hz
bFirLowPassCustom = fir1(3*fix(Fs/hiCutoffFreqCustom), hiCutoffFreqCustom/(Fs/2), 'low');
wbDataLPF100 = filtfilt(bFirLowPassCustom, 1, wbData')';
wbDataLPF100DS = downsample(wbDataLPF100, round(wbFs / lfpFs));

%% plot comparisons
yBounds = [-0.2 0.2];
yBoundsSpectrum = [-90 -45];

figure_tr_inch(20, 10);

subaxis(2, 3, 1);
plot(tWB, wbData);
xlim(tWB([1 end]));
ylim(yBounds);

subaxis(2, 3, 4);
plot(tLFP, lfpData);
xlim(tLFP([1 end]));
ylim(yBounds);

subaxis(2, 3, 2);
hold on;
[S,f] = mtspectrumc(wbDataLPF300DS, params);
plot(f, 10*log10(S));
plot(movmean(f, round(numel(f)/50)), movmean(10*log10(S), round(numel(f)/50)), 'LineWidth', 3);
xlim(params.fpass);
ylim(yBoundsSpectrum);

subaxis(2, 3, 5);
plot(tLFP, wbDataLPF300DS);
xlim(tLFP([1 end]));
ylim(yBounds);

subaxis(2, 3, 3);
hold on;
[S,f] = mtspectrumc(wbDataLPF100DS, params);
plot(f, 10*log10(S));
plot(movmean(f, round(numel(f)/50)), movmean(10*log10(S), round(numel(f)/50)), 'LineWidth', 3);
xlim(params.fpass);
ylim(yBoundsSpectrum);

subaxis(2, 3, 6);
plot(tLFP, wbDataLPF100DS);
xlim(tLFP([1 end]));
ylim(yBounds);

suptitle(sprintf('Channel %d', channelInd));
