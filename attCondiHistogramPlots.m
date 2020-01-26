clear
close all
outputDir = '/Users/labmanager/Documents/MATLAB/BurstSep4/';
v = 14;

fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigMcCartney.csv.bak');
% fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigFerdy.csv.bak');
sessionInfo = textscan(fid, '%d8%s%s', 'Delimiter', ',', 'HeaderLines' ,1-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fid);
counterFR = zeros(1,numel(sessionInfo{1}));
nUnitsPerSession = zeros(1,numel(sessionInfo{1}));

binWidth4hist = 10; % bin width of histograms in ms
idx = 1;

for sessioni = 1:numel(sessionInfo{1})
    clearvars -except sessionInfo sessioni outputDir counterFR ...
        nUnitsPerSession v burstsession burstunit cueArrayAttDiff...
        binWidth4hist cueTargetAttDiff idx
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
    
    %         fig_h1 = figure;
    %         fig_h2 = figure;
    %         for chani = 1:32
    %             cueOnset.windowOffset = [-0.7 0.7];
    %             cueOnset = alignLfpToEvents(D.adjLfps(chani,:), UE.cueOnset, D.lfpFs, cueOnset);
    %
    %             for triali = 1:size(cueOnset.lfp(1,:,:),2)
    %                 lfpFFT(triali,:) = fft(cueOnset.lfp(1,triali,:));
    %             end
    %
    %             meanpow(chani,:) = squeeze(mean(imag(lfpFFT).^2 + real(lfpFFT).^2,1));
    %             figure(fig_h1); hold on; plot(meanpow(chani,:))
    %             xlim([0 50])
    %             [Frac{chani}] = getFractalOscillatoryComponents(squeeze(cueOnset.lfp));
    %             figure(fig_h2); hold on; plot(mean(Frac{chani}.osci,2))
    %         end
    
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
        clear spikeTimes binarySpikeTrain targetDimTrialAtt targetDimTrialAttOut ...
            arrayOnsetTrialAttIn arrayOnsetTrialAttOut cueOnsetTrialAttIn cueOnsetTrialAttOut
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
            
            nUnitsPerSession(sessioni) = nUnits;
            
            %     saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            %             processedDataDir, unitName, blockName, v);
            if firingRateOverall >= minFiringRateOverall
                fprintf('\tOverall firing rate = %0.2f Hz > minimum firing rate = %0.2f Hz in these blocks.\n', ...
                    firingRateOverall, minFiringRateOverall);
                %         fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);
                cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
                
                % Match spikeTimes with first fixation event
                %                         if isnan(unitStruct.unitStartTime)
                startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>spikeTimes(1),1,'first'));
                %                         else
                %                             startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>unitStruct.unitStartTime,1,'first'));
                %                         end
                
                % Match spikeTimes with last lever release
                %                         if isnan(unitStruct.unitEndTime)
                endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<spikeTimes(end),1,'last'));
                %                         else
                %                             endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<unitStruct.unitEndTime,1,'last'));
                %                         end
                
                spikeTimes = spikeTimes(spikeTimes >= startTime & spikeTimes <= endTime);
                spikeTimes2use = (spikeTimes - spikeTimes(1))';
                
                % calculate percentage of spikes >= 200 Hz firing rate
                percentageBurst =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);
                
                
                if percentageBurst > 1
                    targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
                    targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
                    targetDimTrialAttOut = targetDimTrialAttOut - spikeTimes(1);
                    targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut>0);

                    targetDimTrialAttIn = UE.targetDimByLoc{3};
                    targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
                    targetDimTrialAttIn = targetDimTrialAttIn - spikeTimes(1);
                    targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn>0);

                    arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
                    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
                    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - spikeTimes(1);
                    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut>0);

                    arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
                    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
                    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - spikeTimes(1);
                    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn>0);

                    cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
                    cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
                    cueOnsetTrialAttOut = cueOnsetTrialAttOut - spikeTimes(1);
                    cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut>0);

                    cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
                    cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
                    cueOnsetTrialAttIn = cueOnsetTrialAttIn - spikeTimes(1);
                    cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn>0);
                    
                    data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
                    binarySpikeTrain = zeros(1,length(data_pts));
                    binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;
                    
                    % 200ms post-cue until TARGET DIM
                    data = binarySpikeTrain;
                    Fs = 1000;
                    NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                    nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                    datatmpAttIn=[];datatmpAttOut=[];
                    for n=1:NEAttIn;
                    %     nwinl=round(0.200*Fs);
                        nwinr=round(targetDimTrialAttIn(n)*Fs);
                        indx=nEAttIn(n):nwinr-1;
                        if length(indx) >1 && indx(end) < size(data,2)
                            datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                            nSpikesAttInPerTrial(n) = sum(data(indx));
                            clear indx
                        end
                    end
                    for n=1:NEAttOut;
                    %     nwinl=round(0.200*Fs);
                        nwinr=round(targetDimTrialAttOut(n)*Fs);
                        indx=nEAttOut(n):nwinr-1;
                        if length(indx) >1 && indx(end) < size(data,2)
                            datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                        end
                    end

%                         sum(datatmpAttIn(find(datatmpAttIn(1:end-1)>=100) + 1) <= 5) * 100 / size(datatmpAttIn,2)
%                         sum(datatmpAttOut(find(datatmpAttOut(1:end-1)>=100) + 1) <= 5) * 100 / size(datatmpAttOut,2)
%                         sum(datatmpAttIn(find(datatmpAttIn(1:end-1)>=50) + 1) <= 5) * 100 / size(datatmpAttIn,2)
%                         sum(datatmpAttOut(find(datatmpAttOut(1:end-1)>=50) + 1) <= 5) * 100 / size(datatmpAttOut,2)
                    if ~isempty(datatmpAttIn) & ~isempty(datatmpAttOut)
                        figure
                        subplot(131)
                        histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
                        ylimsub1 = ylim;
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
                        ylimsub2 = ylim;
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
                        title('overlay both conditions')
                        xlim([-10 500])
                        subplot(131)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(132)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(133)
                        ylim([0 max([ylimsub1 ylimsub2])])

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   Repeat but with trial based SCM
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        nSpikesAttIn = length(datatmpAttIn);
                        nSpikesAttOut = length(datatmpAttOut);
                        counter = 0;
                        trialSel = randperm(NEAttIn);
                        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
                            counter = counter + 1;
                            NEAttInNw = NEAttIn - counter;
                            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); targetDimTrialAttInNw = targetDimTrialAttIn(trialSel(1:NEAttInNw));
                            for n=1:NEAttInNw;
                                %     nwinl=round(0.200*Fs);
                                nwinr=round(targetDimTrialAttInNw(n)*Fs);
                                indx=nEAttInNw(n):nwinr-1;
                                if length(indx) >1 && indx(end) < size(data,2)
                                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                                end
                            end
                            datatmpAttIn2use = datatmpAttIn;
                            nSpikesAttIn = length(datatmpAttIn2use);
                            datatmpAttIn = [];
                        end
                        if exist('datatmpAttIn2use','var') == 0
                            datatmpAttIn2use = datatmpAttIn;
                        end
                        figure
                        subplot(131)
                        histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        ylimsub1 = ylim;
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        ylimsub2 = ylim;
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        attIn4statsCT = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        hold on
                        attOut4statsCT = histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        title('overlay both conditions')
                        xlim([-10 500])
                        subplot(131)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(132)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(133)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        text(20,max([ylimsub1 ylimsub2])*0.9, ...
                        {sprintf('Womelsdorf method (200Hz) \n burst %%: %0.2f', percentageBurst),...
                        sprintf('First bin %d ms \n AttIn %0.4f - AttOut %0.4f \n = %0.4f', ...
                        binWidth4hist, attIn4statsCT.Values(1), attOut4statsCT.Values(1), attIn4statsCT.Values(1) - attOut4statsCT.Values(1))});                    
                        saveas(gcf,['HistogramISIrealCueTargetDimSCM' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '_' num2str(binWidth4hist) 'ms.png'])

                        cueTargetAttDiff(idx) = attIn4statsCT.Values(1) - attOut4statsCT.Values(1);

                        % 200ms post-cue until ARRAY ONSET
                        data = binarySpikeTrain;
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        for n=1:NEAttIn;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttIn(n)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1 && indx(end) < size(data,2)
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                                nSpikesAttInPerTrialao(n) = sum(data(indx));
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1 && indx(end) < size(data,2)
                                datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                            end
                        end
                        figure
                        subplot(131)
                        histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
                        ylimsub1 = ylim;
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
                        ylimsub2 = ylim;
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
                        title('overlay both conditions')
                        xlim([-10 500])
                        subplot(131)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(132)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(133)
                        ylim([0 max([ylimsub1 ylimsub2])])

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   Repeat but with trial based SCM
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        nSpikesAttIn = length(datatmpAttIn);
                        nSpikesAttOut = length(datatmpAttOut);
                        counter = 0;
                        trialSel = randperm(NEAttIn);
                        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
                            counter = counter + 1;
                            NEAttInNw = NEAttIn - counter;
                            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); arrayOnsetTrialAttInNw = arrayOnsetTrialAttIn(trialSel(1:NEAttInNw));
                            for n=1:NEAttInNw;
                                %     nwinl=round(0.200*Fs);
                                nwinr=round(arrayOnsetTrialAttInNw(n)*Fs);
                                indx=nEAttInNw(n):nwinr-1;
                                if length(indx) >1 && indx(end) < size(data,2)
                                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                                end
                            end
                            datatmpAttIn2use = datatmpAttIn;
                            nSpikesAttIn = length(datatmpAttIn2use);
                            datatmpAttIn = [];
                        end
                        figure
                        subplot(131)
                        histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        ylimsub1 = ylim;
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        ylimsub2 = ylim;
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        attIn4stats = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        hold on
                        attOut4stats = histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
                        title('overlay both conditions')
                        xlim([-10 500])
                        subplot(131)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(132)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        subplot(133)
                        ylim([0 max([ylimsub1 ylimsub2])])
                        text(20,max([ylimsub1 ylimsub2])*0.9, ...
                        {sprintf('Womelsdorf method (200Hz) \n burst %%: %0.2f', percentageBurst),...
                        sprintf('First bin %d ms \n AttIn %0.4f - AttOut %0.4f \n = %0.4f', ...
                        binWidth4hist, attIn4stats.Values(1), attOut4stats.Values(1), attIn4stats.Values(1) - attOut4stats.Values(1))});                    
                        saveas(gcf,['HistogramISIrealCueArrayOnsetSCM' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '_' num2str(binWidth4hist) 'ms.png'])
                    
                        cueArrayAttDiff(idx) = attIn4stats.Values(1) - attOut4stats.Values(1);

                        attIn4statsAll(idx).Values = attIn4stats.Values;
                        attOut4statsAll(idx).Values = attOut4stats.Values;
                        attIn4statsCTAll(idx).Values = attIn4statsCT.Values;
                        attOut4statsCTAll(idx).Values = attOut4statsCT.Values;
                        wdPercentageAll(idx) = percentageBurst;

                        idx = idx + 1;
                        close all
                    end
                    clear binarySpikeTrain
                end
            end
        end
    end
end

save(['AttInMinusAttOut_' num2str(binWidth4hist) 'ms'],'cueArrayAttDiff','cueTargetAttDiff')
[H,P,CI,STATS] = ttest(cueArrayAttDiff)
[H,P,CI,STATS] = ttest(cueTargetAttDiff)
