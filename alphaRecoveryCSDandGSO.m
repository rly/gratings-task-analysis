close all

clear
outputDir = '/Users/labmanager/Documents/MATLAB/BurstSep4/';
v = 14;

fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigMcCartney.csv.bak');
% fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigFerdy.csv.bak');
sessionInfo = textscan(fid, '%d8%s%s', 'Delimiter', ',', 'HeaderLines' ,1-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fid);
counterFR = zeros(1,numel(sessionInfo{1}));
nUnitsPerSession = zeros(1,numel(sessionInfo{1}));

burstsession = [3 15 18 19 21 22 23 25 28 30 31 33];
burstunit = [0 0 7;
    %    0 0 2;
    0 0 1;
    0 0 4;
    1 2 7;
    %    0 1 3;
    1 2 5;
    0 0 33;
    0 2 9;
    0 0 9;
    0 0 4;
    0 0 6;
    0 0 2;
    0 0 3];
for sessioni = 31%1:numel(sessionInfo{1})
    if ismember(sessioni,burstsession)
        clearvars -except sessionInfo sessioni outputDir counterFR nUnitsPerSession v burstsession burstunit
        sessionInd = sessionInfo{1}(sessioni);
        processedDataRootDir = '/Volumes/scratch/rly/gratings-task-analysis/processed_data/';
        dataDirRoot = '/Volumes/kastner/ryanly/McCartney/merged';
        %     dataDirRoot = '/Volumes/kastner/ryanly/Ferdy/merged';
        muaDataDirRoot = '/Volumes/scratch/rly/simple-mua-detection/processed_data/';
        suaMuaDataDirRoot = muaDataDirRoot;
        recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
        channels = strsplit(sessionInfo{3}{sessioni},'-');
        channelsToLoad = str2double(channels{1}):str2double(channels{2});
        muaChannelsToLoad = channelsToLoad;
        lfpChannelsToLoad = str2double(channels{1}):str2double(channels{2});
        lfpChannels = lfpChannelsToLoad;
        isZeroDistractors = 0;
        numRandomizations = 2;
        isLoadSortedSua = 1;
        isLoadMua = 0;
        isLoadAllSpikes = 1;
        %%
        if isZeroDistractors
            scriptName = 'SUA_MUA_GRATINGS_0D';
            taskName = 'GRATINGS_0D';
        else
            scriptName = 'SUA_MUA_GRATINGS';
            taskName = 'GRATINGS';
        end
        %% Plot raster, firing rate, pre-post ISI, waveform/hist & BRI
        [R, D, processedDataDir, blockName] = loadRecordingData(processedDataRootDir, ...
            dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ...
            taskName, scriptName, isLoadSortedSua, isLoadMua, 1, 0, ...
            [],[],isLoadAllSpikes);
        
        sessionName = R.sessionName;
        fprintf('Processing %s...\n', sessionName);
        dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
        gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));
        UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D, blockName);
        dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
        gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));
        % totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
        minFiringRateOverall = 5; % Hz
        nUnits = numel(D.allUnitStructs);
        % kCluster = 4;
        % centroid = zeros(nUnits,kCluster,2);
        
        % apply CSD over 3 channels at the same time
        fileNamePrefix = sprintf('%s-ch%d-ch%d-%s', sessionName, lfpChannels([1 end]), blockName);
        v = 13;
        processedDataDirLFP = '/Volumes/scratch/rly/gratings-task-analysis/processed_data//M20170618/LFP_GRATINGS';
        saveFileName = sprintf('%s/%s-evokedLfps-v%d.mat', ...
                processedDataDirLFP, fileNamePrefix, v);
        fprintf('Loading file %s...\n', saveFileName);
        EL = load(saveFileName);

        nLoc = 4;
        nChannels = D.nLfpCh;
        yVals = 2:nChannels-1;
        xBounds = [-0.4 0.4];
        yBounds = yVals([1 end]) + [-0.5 0.5];
        cBounds = [-0.5 0.5];


        meanCueOnsetLfp = squeeze(mean(EL.cueOnsetLfp.lfp, 2));
        avgCueOnsetCSD = nan(numel(yVals), size(meanCueOnsetLfp, 2));
        for j = 1:numel(yVals)
            ji = yVals(j);
            avgCueOnsetCSD(j,:) = meanCueOnsetLfp(ji+1,:) - 2 * meanCueOnsetLfp(ji,:) + meanCueOnsetLfp(ji-1,:);
        %     cueOnsetCSD(j,:) = cueOnsetCSD(j,:) - mean(cueOnsetCSD(j,preCueBaselineIndices));
        end

        % average 3 neighboring channels
        subCueOnsetLfp = EL.cueOnsetLfp.lfp(yVals,:,:);
        endChannel = 3:3:30; startChannel = 1:3:28;
        chanAvgCueOnset = nan(length(startChannel),size(subCueOnsetLfp,2),size(subCueOnsetLfp,3));
        for ci = 1:length(startChannel)
            chanAvgCueOnset(ci,:,:) = mean(subCueOnsetLfp(startChannel(ci):endChannel(ci),:,:),1);
        end
        yValsChanAvg = 1:10;
        chanAvgMeanCueOnsetLfp = squeeze(mean(chanAvgCueOnset, 2));
        chanAvgCueOnsetCSD = nan(numel(yValsChanAvg), size(chanAvgMeanCueOnsetLfp, 2));
        for j = 2:numel(yValsChanAvg)-1
            chanAvgCueOnsetCSD(j,:) = chanAvgMeanCueOnsetLfp(j+1,:) - 2 * chanAvgMeanCueOnsetLfp(j,:) + chanAvgMeanCueOnsetLfp(j-1,:);
        end

        % try gram-schmidt orthonormalization
        avgCueOnsetGSO = nan(numel(yVals), size(meanCueOnsetLfp, 2));
        for j = 1:numel(yVals)
            ji = yVals(j);
            avgCueOnsetGSO(j,:) = orthogonalizeTimeDomain(meanCueOnsetLfp(ji,:),meanCueOnsetLfp(ji-1,:) + meanCueOnsetLfp(ji+1,:));
        end

        chanAvgCueOnsetGSO = nan(numel(yValsChanAvg), size(meanCueOnsetLfp, 2));
        for j = 1:numel(yValsChanAvg)-1
            ji = yVals(endChannel(j));
            chanAvgCueOnsetGSO(j,:) = orthogonalizeTimeDomain(meanCueOnsetLfp(ji,:),mean(meanCueOnsetLfp(ji-3:ji,:),1) + mean(meanCueOnsetLfp(ji+3,:),1) );
        end


        figure
        subplot(311)
        imagesc(EL.cueOnsetLfp.t, yVals,squeeze(mean(EL.cueOnsetLfp.lfp(yVals,:,:),2)))
        title('raw LFP')
        subplot(323)
        imagesc(EL.cueOnsetLfp.t, yVals,avgCueOnsetCSD)
        title('CSD')
        subplot(324)
        imagesc(EL.cueOnsetLfp.t, yValsChanAvg(2:end-1),chanAvgCueOnsetCSD(2:end-1,:))
        title('CSD over 3Ch')
        subplot(325)
        imagesc(EL.cueOnsetLfp.t, yVals,avgCueOnsetGSO)
        title('Gram-Schmidt normalized')
        subplot(326)
        imagesc(EL.cueOnsetLfp.t, yValsChanAvg(2:end-1),chanAvgCueOnsetGSO(2:end-1,:))
        title('Gram-Schmidt normalized over 3Ch')


% meanCueOnsetLfpByLoc = cell(nLoc, 1);
% preCueBaselineIndices = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, EL.preCueBaselineWindowOffset);
% 
% figure_tr_inch(20, 10); 
% clf;
% set(gcf, 'renderer', 'painters');
% 
% for m = 1:nLoc
%     meanCueOnsetLfpByLoc{m} = squeeze(mean(EL.cueOnsetLfp.lfp(:,UE.cueLoc == m,:), 2));
%     cueOnsetCSD = nan(numel(yVals), size(meanCueOnsetLfpByLoc{m}, 2));
%     for j = 1:numel(yVals)
%         ji = yVals(j);
%         cueOnsetCSD(j,:) = meanCueOnsetLfpByLoc{m}(ji+1,:) - 2 * meanCueOnsetLfpByLoc{m}(ji,:) + meanCueOnsetLfpByLoc{m}(ji-1,:);
%         cueOnsetCSD(j,:) = cueOnsetCSD(j,:) - mean(cueOnsetCSD(j,preCueBaselineIndices));
%     end
%     
%     subaxis(nLoc, 4, (m-1)*4 + 1, 'SH', 0.02, 'MT', 0.06);
%     hold on;
%     imagesc(EL.cueOnsetLfp.t, yVals, cueOnsetCSD);
%     plot([0 0], yBounds, 'Color', 0.3*ones(3, 1));
%     set(gca, 'YDir', 'reverse');
%     xlim(xBounds);
%     ylim(yBounds);
%     caxis(cBounds);
%     if m == 1
%         title('Cue Onset');
%     end
%     if m == nLoc
%         xlabel('Time from Cue Onset (s)');
%     end
% end

        fig_h1 = figure;
        fig_h2 = figure;
        for chani = 1:32
            cueOnset.windowOffset = [-0.7 0.7];
            cueOnset = alignLfpToEvents(D.adjLfps(chani,:), UE.cueOnset, D.lfpFs, cueOnset);
            
            for triali = 1:size(cueOnset.lfp(1,:,:),2)
                lfpFFT(triali,:) = fft(cueOnset.lfp(1,triali,:));
            end
            
            meanpow(chani,:) = squeeze(mean(imag(lfpFFT).^2 + real(lfpFFT).^2,1));
            figure(fig_h1); hold on; plot(meanpow(chani,:))
            xlim([0 50])
            [Frac{chani}] = getFractalOscillatoryComponents(squeeze(cueOnset.lfp));
            figure(fig_h2); hold on; plot(mean(Frac{chani}.osci,2))
        end      
    end
end






%eof