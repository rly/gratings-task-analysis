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
for sessioni = 15%31%1:numel(sessionInfo{1})
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
                        targetDimTrialAttOut = UE.targetDimByLoc{1};
                        targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
                        targetDimTrialAttOut = targetDimTrialAttOut - spikeTimes(1); 

                        targetDimTrialAttIn = UE.targetDimByLoc{3};
                        targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
                        targetDimTrialAttIn = targetDimTrialAttIn - spikeTimes(1); 
    
                        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
                        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
                        arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - spikeTimes(1); 

                        arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
                        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
                        arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - spikeTimes(1); 
                        
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
                        cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
                        cueOnsetTrialAttOut = cueOnsetTrialAttOut - spikeTimes(1); 

                        cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
                        cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
                        cueOnsetTrialAttIn = cueOnsetTrialAttIn - spikeTimes(1); 
                        
                        spikeTimes2use = (spikeTimes - spikeTimes(1))';
                        data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
                        binarySpikeTrain = zeros(1,length(data_pts)); 
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

                        % 200ms post-cue until TARGET DIM                      
                        data = binarySpikeTrain; 
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        for n=1:NEAttIn;
                            nwinl=round(arrayOnsetTrialAttIn(n)*Fs) - 1200;
                            nwinr=round(arrayOnsetTrialAttIn(n)*Fs) + 1200;
                            indx=nwinl:nwinr-1;
                            if length(indx) >1 && indx(end) < size(data,2)
                                binarySTarrayLocked(n,:)=data(indx);
%                                 nSpikesAttInPerTrial(n) = sum(data(indx));
                                binarySTarrayLocked2plot(n,:) = zeros(1,2400);
                                allSpikes = find(binarySTarrayLocked(n,:));
                                tenmsISI = find(diff(find(binarySTarrayLocked(n,:)))<=10);
                                hundredmsISI = find(diff(find(binarySTarrayLocked(n,:)))>10 & diff(find(binarySTarrayLocked(n,:))) <= 50);
                                for isi = 1:length(tenmsISI)
                                    binarySTarrayLocked2plot(n,allSpikes(tenmsISI(isi)):allSpikes(tenmsISI(isi)+1)-1) = 1;
                                end
                                for isii = 1:length(hundredmsISI)
                                    binarySTarrayLocked2plot(n,allSpikes(hundredmsISI(isii)):allSpikes(hundredmsISI(isii)+1)-1) = 0.5;
                                end
                                clear indx
                            end
                            cueOn(n) = 1200 - (round(arrayOnsetTrialAttIn(n)*Fs) - round(cueOnsetTrialAttIn(n)*Fs));
                            tdOn(n) = 1200 + (round(targetDimTrialAttIn(n)*Fs) - round(arrayOnsetTrialAttIn(n)*Fs));
                            binarySTarrayLocked2plot(n,1:cueOn(n)) = -0.5; binarySTarrayLocked2plot(n,tdOn(n):end) = -0.5;
                        end
                        
                        for n=1:NEAttOut;
                            nwinl=round(arrayOnsetTrialAttOut(n)*Fs) - 1200;
                            nwinr=round(arrayOnsetTrialAttOut(n)*Fs) + 1200;
                            indx=nwinl:nwinr-1;
                            if length(indx) >1 && indx(end) < size(data,2)
                                binarySTarrayLockedAttOut(n,:)=data(indx);
%                                 nSpikesAttInPerTrial(n) = sum(data(indx));
                                binarySTarrayLocked2plotAttOut(n,:) = zeros(1,2400);
                                allSpikes = find(binarySTarrayLockedAttOut(n,:));
                                tenmsISI = find(diff(find(binarySTarrayLockedAttOut(n,:)))<=10);
                                hundredmsISI = find(diff(find(binarySTarrayLockedAttOut(n,:)))>10 & diff(find(binarySTarrayLockedAttOut(n,:))) <= 50);
                                for isi = 1:length(tenmsISI)
                                    binarySTarrayLocked2plotAttOut(n,allSpikes(tenmsISI(isi)):allSpikes(tenmsISI(isi)+1)-1) = 1;
                                end
                                for isii = 1:length(hundredmsISI)
                                    binarySTarrayLocked2plotAttOut(n,allSpikes(hundredmsISI(isii)):allSpikes(hundredmsISI(isii)+1)-1) = 0.5;
                                end
                                clear indx
                            end
                            cueOnAttOut(n) = 1200 - (round(arrayOnsetTrialAttOut(n)*Fs) - round(cueOnsetTrialAttOut(n)*Fs));
                            tdOnAttOut(n) = 1200 + (round(targetDimTrialAttOut(n)*Fs) - round(arrayOnsetTrialAttOut(n)*Fs));
                            binarySTarrayLocked2plotAttOut(n,1:cueOnAttOut(n)) = -0.5; binarySTarrayLocked2plotAttOut(n,tdOnAttOut(n):end) = -0.5;
                        end
                        figure
                        timeaxis = -1200:1199;
                        subplot(121)
                        imagesc(timeaxis,1:size(binarySTarrayLocked2plot,1),binarySTarrayLocked2plot)
                        title('Attend In')
                        subplot(122)
                        imagesc(timeaxis,1:size(binarySTarrayLocked2plotAttOut,1),binarySTarrayLocked2plotAttOut)
                        title('Attend Out')
% cmap = [0 0 0; 0.7 0.7 0.7; 0.3 0.3 0.3; 1 0 0];
% colormap(cmap)
                        % white indicates short ISIs of <10ms
                        % light gray are ISIs of tonic spiking >10ms and
                        % <100ms
                        % dark gray are ISIs of quiet periods >100ms
                        % black indicates time before 200ms post-cue and
                        % after target dim

                        saveas(gcf,['trialHeatmapISI' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])

                        clear burstVect burstEventOnset spikeTimes AllBurstEvents
                    else
                        counterFR(sessioni) = 1;
                    end
                end
            end
        end
    end
end

