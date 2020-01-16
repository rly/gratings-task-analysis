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
for sessioni = 1:numel(sessionInfo{1})
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
            [~,iidx] = ismember(sessioni,burstsession);
            if ismember(uniti,burstunit(iidx,:))
                
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
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        
                        % Match spikeTimes with first fixation event
                        if isnan(unitStruct.unitStartTime)
                            startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>spikeTimes(1),1,'first'));
                        else
                            startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>unitStruct.unitStartTime,1,'first'));
                        end
                        
                        % Match spikeTimes with last lever release
                        if isnan(unitStruct.unitEndTime)
                            endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<spikeTimes(end),1,'last'));
                        else
                            endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<unitStruct.unitEndTime,1,'last'));
                        end
                            
                        spikeTimes = spikeTimes(spikeTimes >= startTime & spikeTimes <= endTime);
                        spikeTimes2use = (spikeTimes - spikeTimes(1))';
                        data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
                        binarySpikeTrain = zeros(1,length(data_pts)); 
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

                        % 200ms post-cue until TARGET DIM
                        targetDimTrialAttOut = UE.targetDimByLoc{1};% - spikeTimes(1); 
                        targetDimTrialAttIn = UE.targetDimByLoc{3};% - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;% - spikeTimes(1); 
                        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;% - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; 
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        datatmpAttInGAKS=[];datatmpAttOutGAKS=[];
                        datatmpAttInGAKSsep=[];datatmpAttOutGAKSsep=[];
                        datatmpShuffledAttInGAKS=[];datatmpShuffledAttOutGAKS=[];
                        datatmpShuffledAttInGAKSsep=[];datatmpShuffledAttOutGAKSsep=[];
                        datatmpAttInShuffled = []; datatmpAttOutShuffled = [];
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
                        
                        figure
                        subplot(131)
                        histogram(datatmpAttIn,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn,'BinWidth',10,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        title('overlay both conditions')
                        ylim([0 .22])
                        xlim([-10 500])
                        saveas(gcf,['HistogramISIrealCueTargetDim' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        
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
                        histogram(datatmpAttIn2use,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn2use,'BinWidth',10,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        title('overlay both conditions')
                        ylim([0 .22])
                        xlim([-10 500])
                        saveas(gcf,['HistogramISIrealCueTargetDimSCM' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])

                        
                        % 200ms post-cue until ARRAY ONSET
                        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);% - spikeTimes(1); 
                        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);% - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;% - spikeTimes(1); 
                        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;% - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; 
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        datatmpAttInGAKS=[];datatmpAttOutGAKS=[];
                        datatmpAttInGAKSsep=[];datatmpAttOutGAKSsep=[];
                        datatmpShuffledAttInGAKS=[];datatmpShuffledAttOutGAKS=[];
                        datatmpShuffledAttInGAKSsep=[];datatmpShuffledAttOutGAKSsep=[];
                        datatmpAttInShuffled = []; datatmpAttOutShuffled = [];
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
                        histogram(datatmpAttIn,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn,'BinWidth',10,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        title('overlay both conditions')
                        ylim([0 .22])
                        xlim([-10 500])
                        saveas(gcf,['HistogramISIrealCueArrayOnset' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                       
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
                        histogram(datatmpAttIn2use,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att In')
                        subplot(132)
                        plot(0,0)
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        ylim([0 .22])
                        xlim([-10 500])
                        title('real data Att Out')
                        subplot(133)
                        histogram(datatmpAttIn2use,'BinWidth',10,'Normalization','probability')
                        hold on
                        histogram(datatmpAttOut,'BinWidth',10,'Normalization','probability')
                        title('overlay both conditions')
                        ylim([0 .22])
                        xlim([-10 500])
                        saveas(gcf,['HistogramISIrealCueArrayOnsetSCM' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])

                        clear binarySpikeTrain

                    else
                        counterFR(sessioni) = 1;
                    end
                end
            end
        end
    end
end

