clear
close all
outputDir = '/Users/labmanager/Documents/MATLAB/BurstSep5/';
v = 14;

fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigMcCartney.csv.bak');
% fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigFerdy.csv.bak');
sessionInfo = textscan(fid, '%d8%s%s', 'Delimiter', ',', 'HeaderLines' ,2-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fid);
for sessioni = 1:3%numel(sessionInfo{1})
    clear -sessionInfo -sessioni
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
        taskName, scriptName, isLoadSortedSua, isLoadMua, 0, 0, ...
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
    
    for uniti = 1:nUnits
        unitStruct = D.allUnitStructs{uniti};
        unitName = unitStruct.name;
        fprintf('-----------------------------\n');
        fprintf('\n');
        if ~isempty(unitStruct.ts)
            
            spikeTimes = D.allUnitStructs{uniti}.ts;
            totalTimeOverall = spikeTimes(end) - spikeTimes(1);
            
            %                     spikeTimes(50000:end) = spikeTimes(50000:end) + 10000;
            
            firingRateOverall = numel(spikeTimes) / totalTimeOverall;
            fprintf('SPNA %s (%d/%d = %d%%)... \n', unitName, uniti, ...
                nUnits, round(uniti/nUnits*100));
            
            %     saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            %             processedDataDir, unitName, blockName, v);
            if firingRateOverall >= minFiringRateOverall
                fprintf('\tOverall firing rate = %0.2f Hz > minimum firing rate = %0.2f Hz in these blocks.\n', ...
                    firingRateOverall, minFiringRateOverall);
                %         fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);
                
                allDiffSpikeTimes = diff(spikeTimes); % get the ISIs
                allDiffSpikeTimesPre = [NaN; allDiffSpikeTimes]; % miss align ISIs by adding NaN to create pre and postISI
                allDiffSpikeTimesPost = [allDiffSpikeTimes; NaN];
                allDiffSpikeTimesPre = log10(allDiffSpikeTimesPre*1000); % logISI in ms
                allDiffSpikeTimesPost = log10(allDiffSpikeTimesPost*1000);
                
                FigH = figure;
                % return plot
                subplot(221)
                plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','MarkerSize', 0.05)
                preTimesOverThreshLong = find(allDiffSpikeTimesPre > log10(100));
                postTimesOverThreshLong = find(allDiffSpikeTimesPost > log10(100));
                %             loglog(allDiffSpikeTimes(:,1),allDiffSpikeTimes(:,2),'.')
                hold on
                x1 = allDiffSpikeTimesPre(preTimesOverThreshLong);
                y1 = allDiffSpikeTimesPost(preTimesOverThreshLong);
                %             plot(x1,y1,'r.')
                x2 = allDiffSpikeTimesPre(postTimesOverThreshLong);
                y2 = allDiffSpikeTimesPost(postTimesOverThreshLong);
                %             plot(x2,y2,'k.')
                bothTimesOverThresh = find(allDiffSpikeTimesPre(postTimesOverThreshLong) > log10(100));
                x3 = x2(bothTimesOverThresh);
                y3 = y2(bothTimesOverThresh);
                plot(x3,y3,'k.','MarkerSize', 0.05)
                preTimesOverThreshShort = find(allDiffSpikeTimesPre < log10(4));
                postTimesOverThreshShort = find(allDiffSpikeTimesPost < log10(4));
                x4 = allDiffSpikeTimesPre(preTimesOverThreshShort);
                y4 = allDiffSpikeTimesPost(preTimesOverThreshShort);
                %             plot(x4,y4,'b.')
                x5 = allDiffSpikeTimesPre(postTimesOverThreshShort);
                y5 = allDiffSpikeTimesPost(postTimesOverThreshShort);
                %             plot(x5,y5,'m.')
                bothTimesOverThreshShort = find(allDiffSpikeTimesPre(postTimesOverThreshShort) < log10(4));
                x6 = x5(bothTimesOverThreshShort);
                y6 = y5(bothTimesOverThreshShort);
                plot(x6,y6,'m.','MarkerSize', 0.05)
                xlabel('Pre ISI')
                ylabel('Post ISI')
                title('Return plot: All data')
                hold on
                plot([2 2],[log10(1) log10(4)],'--k')
                plot([2 4],[log10(4) log10(4)],'--k')
                %             plot([3 3],[log10(1) log10(4)],'--k')
                plot([log10(4) log10(4)],[0 10],'--k')
                xlim([0 4])
                ylim([0 4])
                xticks(0:4)
                xticklabels({'1', '10', '100', '1000'})
                xlabel('PreISI in ms on log scale')
                yticks(0:4)
                yticklabels({'1', '10', '100', '1000'})
                ylabel('PostISI in ms on log scale')
                
                % calculate % of bursting activity
                %             figure
                %             plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','MarkerSize', 1)
                %             hold on
                %             plot(allDiffSpikeTimesPre(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4)),allDiffSpikeTimesPost(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4)),'m*','MarkerSize', 1)
                firstSpikeBurst = sum(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4),1);
                %             hold on
                %             plot(allDiffSpikeTimesPre(allDiffSpikeTimesPre < log10(4)),allDiffSpikeTimesPost(allDiffSpikeTimesPre < log10(4)),'k*','MarkerSize', 1)
                overFirstSpikeBurst = sum(allDiffSpikeTimesPre < log10(4),1);
                
                allISIs = size(allDiffSpikeTimesPre,1);
                burstISIs = overFirstSpikeBurst + firstSpikeBurst;
                percBurst(uniti) = (burstISIs*100)/allISIs;
                fprintf('\tBurst activity = %0.2f %% of %d spikes.\n', ...
                    percBurst(uniti), allISIs);
                
                set(FigH, 'NumberTitle', 'off', ...
                    'Name', sprintf('Burst activity = %0.2f %% of %d spikes.', percBurst(uniti), allISIs));
                
                timeXaxis = linspace(min(allDiffSpikeTimesPre), max(allDiffSpikeTimesPre),100);
                % plot logISI histogram
                subplot(222)
                %histogram(log10(allDiffSpikeTimesPre(2:end)),timeXaxis)
                histogram(allDiffSpikeTimesPre,timeXaxis)
                xticks(-1:3)
                xticklabels({'0.1','1', '10', '100', '1000'})
                xlabel('ISI in ms on log scale')
                title('Histogram logISI')
                
                % smoothing procedure to detect peaks
                %outputHistogram = histcounts(log10(allDiffSpikeTimesPre(2:end)),timeXaxis);
                outputHistogram = histcounts(allDiffSpikeTimesPre,timeXaxis);
                smoothISIhist = smooth(outputHistogram,'lowess');
                [pks,locs,widths] = findpeaks(smoothISIhist,'MinPeakHeight',2);
                hold on; plot(timeXaxis(1:end-1),smoothISIhist','m'); plot(timeXaxis(locs), pks', 'm*')
                
                subplot(234)
                plot(D.allUnitStructs{uniti}.meanWf,'k','LineWidth',3)
                hold on
                plot(D.allUnitStructs{uniti}.meanWf - (std(D.allUnitStructs{uniti}.wf,1)*3),'--k');
                plot(D.allUnitStructs{uniti}.meanWf + (std(D.allUnitStructs{uniti}.wf,1)*3),'--k');
                title('Spike waveform')
                
                binnedSpikes = histcounts(spikeTimes,'BinWidth',0.001);
                subplot(438)
                [autocorrSpikeTimes,lags] = xcorr(binnedSpikes,1000,'coeff');
                plot(lags,autocorrSpikeTimes)
                title('Normalized Autocorrelogram')
                subplot(4,6,21)
                [autocorrSpikeTimes,lags] = xcorr(binnedSpikes,1000,'coeff');
                plot(lags,autocorrSpikeTimes)
                xlim([-75 75])
                title('Normalized Autocorrelogram Zoomed In')
                
                tcorr=0:0.001:1; % time-vector for the calculation of correlations
                nbins = length(binnedSpikes);
                ncorr = length(tcorr); %number of points of time at which correlation is calculated
                corr=zeros(size(tcorr)); % set up the vector that will contain the correlations
                
                for j=1:ncorr % this will generate ncorr values
                    corr(j) = binnedSpikes(j:nbins)*binnedSpikes(1:nbins+1-j).';
                end
                subplot(4,6,22)
                plot(tcorr,corr)
                title('Autocorrelogram')
                
                Fs = 1000;
                
                spikeTimesAdj = spikeTimes - spikeTimes(1);
                binarySignalInd = spikeTimesAdj(1):(1/Fs):spikeTimesAdj(end);
                binarySignal = zeros(1,size(binarySignalInd,2));
                binarySignal(round(spikeTimesAdj'*Fs)+1) = 1;
                subplot(236)
                welchwin = Fs*8.192;
                [PSDSpikeTimes, freq] = pwelch(binarySignal,welchwin,[],[],Fs);
                plot(freq,PSDSpikeTimes)
                xlim([1 60])
                title('PSD')
                
                clear spikeTimes allDiffSpikeTimes allDiffSpikeTimesPost allDiffSpikeTimesPre ...
                    preTimesOverThreshLong preTimesOverThreshShort bothTimesOverThresh ...
                    postTimesOverThreshLong postTimesOverThreshShort bothTimesOverThreshShort
                
                plotFileName = sprintf('%s/%s-sessionInd%d-Uniti%d-BurstingTvV-v%d.png', outputDir, sessionName, sessionInd, uniti, v);
                fprintf('\tSaving figure to file %s...\n', plotFileName);
                saveas(FigH, plotFileName);
                close all
            end
        end
    end
end



