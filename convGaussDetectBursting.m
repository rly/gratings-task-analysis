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
                        
                        spikeTimes2use = (spikeTimes - spikeTimes(1))';
                        data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
                        binarySpikeTrain = zeros(1,length(data_pts)); 
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

                        % convolve with gaussian
                        kernelSigma = 0.1;
                        x = -5 * kernelSigma * 1000 : 5 * kernelSigma * 1000;
                        gaussian = normpdf(x, 0, kernelSigma * 1000);  % does the same thing as computing a gaussian with an equation as you had done
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;
                        allSpikesGauss = conv(binarySpikeTrain, gaussian,'same');
                        
                        for shuffi = 1:10;
                            shuffledBinarySpikeTrain = binarySpikeTrain(randperm(length(binarySpikeTrain)));
                            shuffledSpikesGauss = conv(shuffledBinarySpikeTrain, gaussian,'same');
                            highestShuf(shuffi) = max(shuffledSpikesGauss);
                        end
                        sepLowHigh = mean(highestShuf);

                        % threshold
                        allSpikesGaussSepTmp = ((allSpikesGauss<sepLowHigh) * 0.5);
                        allSpikesGaussSep = allSpikesGaussSepTmp + (allSpikesGauss==0);
                        
                        shuffledSpikesGaussSepTmp = ((shuffledSpikesGauss<sepLowHigh) * 0.5);
                        shuffledSpikesGaussSep = shuffledSpikesGaussSepTmp + (shuffledSpikesGauss==0);
                        
                        % cut spike times that fall between enter and end
                        % fixation
                        endTimeTrial = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice - spikeTimes(1);
                        startTimeTrial = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue - spikeTimes(1);

                        % heatmap return plot
                        data = binarySpikeTrain; Fs = 1000;
                        NE=length(startTimeTrial);
                        nE=floor(startTimeTrial*Fs)+1;
                        datatmp=[];
                        for n=1:NE;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(endTimeTrial(n)*Fs);
                            indx=nE(n):nwinr-1;
                            if length(indx) >1
                                datatmp=[datatmp diff(find(data(indx)))];
                            end
                        end

                        dataSil = binarySpikeTrain; Fs = 1000;
                        NE=length(endTimeTrial)-1;
                        nE=floor(endTimeTrial*Fs)+1;
                        datatmpSil=[];
                        for n=1:NE;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(startTimeTrial(n+1)*Fs);
                            indx=nE(n):nwinr-1;
                            if length(indx) >1
                                datatmpSil=[datatmpSil diff(find(dataSil(indx)))];
                            end
                        end

                        figure
                        subplot(131)
                        loglog([NaN datatmp],[datatmp NaN],'.','MarkerSize',0.1)
                        title('return plot with spikes between enter and exit fix')
                        subplot(132)
                        loglog([NaN datatmpSil],[datatmpSil NaN],'.','MarkerSize',0.1)
                        title('return plot with spikes outside trial')
                        subplot(133)
                        loglog([NaN diff(find(binarySpikeTrain))],[diff(find(binarySpikeTrain)) NaN],'.','MarkerSize',0.1)                        
                        title('return plot all spikes')
                        saveas(gcf,'returnPlot.png')
                        
                        figure
                        subplot(131)
                        edges = linspace(0,4,100);
                        h = histogram2([NaN log10(datatmp)],[log10(datatmp) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.002 .002])
                        title('return plot with spikes between enter and exit fix')
                        subplot(132)
                        s = histogram2([NaN log10(datatmpSil)],[log10(datatmpSil) NaN], edges,edges,'Normalization','probability')
                        s.FaceColor = 'flat';
                        s.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.002 .002])
                        title('return plot with spikes outside trial')
                        subplot(133)
                        s = histogram2([NaN log10(diff(find(binarySpikeTrain)))],[log10(diff(find(binarySpikeTrain))) NaN], edges,edges,'Normalization','probability')
                        s.FaceColor = 'flat';
                        s.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.002 .002])
                        title('return plot all spikes')
                        saveas(gcf,'returnPlotHeatlikemapNormalized.png')


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %       Repeat heatmap for attention sep Locs
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % use end and start times to select att in and out
                        % trials
                        endTimeTrialAttOut = endTimeTrial(UE.cueLoc == 1); endTimeTrialAttIn = endTimeTrial(UE.cueLoc == 3);
                        startTimeTrialAttOut = startTimeTrial(UE.cueLoc == 1); startTimeTrialAttIn = startTimeTrial(UE.cueLoc == 3);
                        
                        % heatmap return plot
                        data = binarySpikeTrain; Fs = 1000;
                        NEAttIn=length(startTimeTrialAttIn);NEAttOut=length(startTimeTrialAttOut);
                        nEAttIn=floor(startTimeTrialAttIn*Fs)+1;nEAttOut=floor(startTimeTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        for n=1:NEAttIn;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(endTimeTrialAttIn(n)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(endTimeTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                            end
                        end
                        dataSil = binarySpikeTrain; Fs = 1000;
                        NEAttIn=length(endTimeTrialAttIn)-1;NEAttOut=length(endTimeTrialAttOut)-1;
                        nEAttIn=floor(endTimeTrialAttIn*Fs)+1;nEAttOut=floor(endTimeTrialAttOut*Fs)+1;
                        datatmpSilAttIn=[];datatmpSilAttOut=[];
                        for n=1:NEAttIn;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(startTimeTrialAttIn(n+1)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpSilAttIn=[datatmpSilAttIn diff(find(dataSil(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(startTimeTrialAttOut(n+1)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpSilAttOut=[datatmpSilAttOut diff(find(dataSil(indx)))];
                            end
                        end

                        figure
                        subplot(211)
                        loglog([NaN datatmpAttIn],[datatmpAttIn NaN],'.','MarkerSize',0.1)
                        title('return plot with spikes between enter and exit fix Att In')
                        subplot(212)
                        loglog([NaN datatmpAttOut],[datatmpAttOut NaN],'.','MarkerSize',0.1)
                        title('return plot with spikes between enter and exit fix Att Out')
%                         subplot(223)
%                         loglog([NaN datatmpSilAttIn],[datatmpSilAttIn NaN],'.','MarkerSize',0.1)
%                         title('return plot with spikes outside trial Att In')
%                         subplot(224)
%                         loglog([NaN datatmpSilAttOut],[datatmpSilAttOut NaN],'.','MarkerSize',0.1)
%                         title('return plot with spikes outside trial Att Out')
                        saveas(gcf,'returnPlotAttLocs.png')
                        
                        figure
                        subplot(211)
                        edges = linspace(0,4,100);
                        h = histogram2([NaN log10(datatmpAttIn)],[log10(datatmpAttIn) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.002 .002])
                        title('return plot with spikes between enter and exit fix Att In')
                        subplot(212)
                        edges = linspace(0,4,100);
                        h = histogram2([NaN log10(datatmpAttOut)],[log10(datatmpAttOut) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.002 .002])
                        title('return plot with spikes between enter and exit fix Att Out')
%                         subplot(223)
%                         s = histogram2([NaN log10(datatmpSilAttIn)],[log10(datatmpSilAttIn) NaN], edges,edges,'Normalization','probability')
%                         s.FaceColor = 'flat';
%                         s.DisplayStyle = 'tile';
%                         view(2)
%                         xlim([0 4]); ylim([0 4])
%                         caxis([-.002 .002])
%                         title('return plot with spikes outside trial Att In')
%                         subplot(224)
%                         s = histogram2([NaN log10(datatmpSilAttOut)],[log10(datatmpSilAttOut) NaN], edges,edges,'Normalization','probability')
%                         s.FaceColor = 'flat';
%                         s.DisplayStyle = 'tile';
%                         view(2)
%                         xlim([0 4]); ylim([0 4])
%                         caxis([-.002 .002])
%                         title('return plot with spikes outside trial Att Out')
                        saveas(gcf,'returnPlotHeatlikemapNormalizedAttLocs.png')

                        % 200ms post-cue until array onset
                        arrayOnsetTrialAttOut = UE.arrayOnset(UE.cueLoc == 1) - spikeTimes(1); arrayOnsetTrialAttIn = UE.arrayOnset(UE.cueLoc == 3) - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                        
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        for n=1:NEAttIn;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttIn(n)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                            end
                        end
                        dataSil = binarySpikeTrain; Fs = 1000;
                        NEAttIn=length(arrayOnsetTrialAttIn)-1;NEAttOut=length(arrayOnsetTrialAttOut)-1;
                        nEAttIn=floor(arrayOnsetTrialAttIn*Fs)+1;nEAttOut=floor(arrayOnsetTrialAttOut*Fs)+1;
                        datatmpSilAttIn=[];datatmpSilAttOut=[];
                        for n=1:NEAttIn;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(cueOnsetTrialAttIn(n+1)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpSilAttIn=[datatmpSilAttIn diff(find(dataSil(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(cueOnsetTrialAttOut(n+1)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpSilAttOut=[datatmpSilAttOut diff(find(dataSil(indx)))];
                            end
                        end
                        figure
                        subplot(211)
                        edges = linspace(0,4,25);
                        h = histogram2([NaN log10(datatmpAttIn)],[log10(datatmpAttIn) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.04 .04])
                        title('return plot with spikes 200ms post-cue and array onset Att In')
                        subplot(212)
                        h = histogram2([NaN log10(datatmpAttOut)],[log10(datatmpAttOut) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.04 .04])
                        title('return plot with spikes 200ms post-cue and array onset Att Out')
                        saveas(gcf,'returnPlotPostCueArrayOnAttLocs.png')


                        % 200ms post-cue until TARGET DIM
                        targetDimTrialAttOut = UE.targetDimByLoc{1} - spikeTimes(1); targetDimTrialAttIn = UE.targetDimByLoc{3} - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss; 
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
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
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                                datatmpAttInGAKS=[datatmpAttInGAKS dataGAKS(indx)];
                                datatmpAttInGAKSsep=[datatmpAttInGAKSsep dataGAKSsep(indx)];
                                datatmpShuffledAttInGAKS=[datatmpShuffledAttInGAKS dataGAKS2(indx)];
                                datatmpShuffledAttInGAKSsep=[datatmpShuffledAttInGAKSsep dataGAKSsep2(indx)];
                                indx = indx(randperm(length(indx)));
                                datatmpAttInShuffled=[datatmpAttInShuffled diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(targetDimTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                                datatmpAttOutGAKS=[datatmpAttOutGAKS dataGAKS(indx)];
                                datatmpAttOutGAKSsep=[datatmpAttOutGAKSsep dataGAKSsep(indx)];
                                datatmpShuffledAttOutGAKS=[datatmpShuffledAttOutGAKS dataGAKS2(indx)];
                                datatmpShuffledAttOutGAKSsep=[datatmpShuffledAttOutGAKSsep dataGAKSsep2(indx)];
                                indx = indx(randperm(length(indx)));
                                datatmpAttOutShuffled=[datatmpAttOutShuffled diff(find(data(indx)))];
                            end
                        end
                        figure
                        subplot(211)
                        edges = linspace(0,4,25);
                        h = histogram2([NaN log10(datatmpAttIn)],[log10(datatmpAttIn) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.04 .04])
                        title('return plot with spikes 200ms post-cue to target Dim Att In')
                        subplot(212)
                        h = histogram2([NaN log10(datatmpAttOut)],[log10(datatmpAttOut) NaN], edges,edges,'Normalization','probability')
                        h.FaceColor = 'flat';
                        h.DisplayStyle = 'tile';
                        view(2)
                        xlim([0 4]); ylim([0 4])
                        caxis([-.04 .04])
                        title('return plot with spikes 200ms post-cue to target Dim Att Out')
                        saveas(gcf,'returnPlotPostCueTargetDimAttLocs.png')

                        % return plot of all spikes between cue onset and
                        % target dim next to shuffled data for same time
                        % periods
                        % create real and shuffled data for all trials
                        targetDimTrial = UE.targetDim - spikeTimes(1); 
                        cueOnsetTrial = UE.cueOnset(UE.isHoldTrial) + 0.200 - spikeTimes(1);
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss;
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
                        Fs = 1000;
                        NE=length(cueOnsetTrial);
                        nE=floor(cueOnsetTrial*Fs)+1;
                        datatmp=[];datatmpShuffled=[];
                        datatmpGAKS=[];datatmpShuffledGAKS=[];
                        datatmpGAKSsep=[];datatmpShuffledGAKSsep=[];
                        for n=1:NE;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(targetDimTrial(n)*Fs);
                            indx=nE(n):nwinr-1;
                            if length(indx) >1
                                datatmp=[datatmp diff(find(data(indx)))];
                                datatmpGAKS=[datatmpGAKS dataGAKS(indx)];
                                datatmpGAKSsep=[datatmpGAKSsep dataGAKSsep(indx)];
                                datatmpShuffledGAKS=[datatmpShuffledGAKS dataGAKS2(indx)];
                                datatmpShuffledGAKSsep=[datatmpShuffledGAKSsep dataGAKSsep2(indx)];
                                indx = indx(randperm(length(indx)));
                                datatmpShuffled=[datatmpShuffled diff(find(data(indx)))];
                            end
                        end
                        figure
                        subplot(231)
                        histogram(datatmpAttIn,linspace(0,600,100))
                        ylim([0 90])
                        title('real data Att In')
                        subplot(234)
                        histogram(datatmpAttInShuffled,linspace(0,600,100))
                        ylim([0 90])
                        title('shuffled data Att In')
                        subplot(232)
                        histogram(datatmpAttOut,linspace(0,600,100))
                        ylim([0 130])
                        title('real data Att Out')
                        subplot(235)
                        histogram(datatmpAttOutShuffled,linspace(0,600,100))
                        ylim([0 130])
                        title('shuffled data Att Out')
                        subplot(233)
                        histogram(datatmp,linspace(0,600,100))
                        ylim([0 500])
                        title('real data all')
                        subplot(236)
                        histogram(datatmpShuffled,linspace(0,600,100))
                        ylim([0 500])
                        title('shuffled data all')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramISIrealShuffCueTargetDim.png')

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
                        saveas(gcf,'HistogramISIrealCueTargetDim.png')

                        % 200ms post-cue until ARRAY ONSET
                        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1) - spikeTimes(1); arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3) - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss; 
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
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
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                                datatmpAttInGAKS=[datatmpAttInGAKS dataGAKS(indx)];
                                datatmpAttInGAKSsep=[datatmpAttInGAKSsep dataGAKSsep(indx)];
                                datatmpShuffledAttInGAKS=[datatmpShuffledAttInGAKS dataGAKS2(indx)];
                                datatmpShuffledAttInGAKSsep=[datatmpShuffledAttInGAKSsep dataGAKSsep2(indx)];
                                indx = indx(randperm(length(indx)));
                                datatmpAttInShuffled=[datatmpAttInShuffled diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                                datatmpAttOutGAKS=[datatmpAttOutGAKS dataGAKS(indx)];
                                datatmpAttOutGAKSsep=[datatmpAttOutGAKSsep dataGAKSsep(indx)];
                                datatmpShuffledAttOutGAKS=[datatmpShuffledAttOutGAKS dataGAKS2(indx)];
                                datatmpShuffledAttOutGAKSsep=[datatmpShuffledAttOutGAKSsep dataGAKSsep2(indx)];
                                indx = indx(randperm(length(indx)));
                                datatmpAttOutShuffled=[datatmpAttOutShuffled diff(find(data(indx)))];
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
                        saveas(gcf,'HistogramISIrealCueArrayOnset.png')

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   Repeat but with trial based SCM
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % 200ms post-cue until TARGET DIM
                        targetDimTrialAttOut = UE.targetDimByLoc{1} - spikeTimes(1); targetDimTrialAttIn = UE.targetDimByLoc{3} - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss; 
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        datatmpAttInGAKS=[];datatmpAttOutGAKS=[];
                        datatmpAttInGAKSsep=[];datatmpAttOutGAKSsep=[];
                        datatmpShuffledAttInGAKS=[];datatmpShuffledAttOutGAKS=[];
                        datatmpShuffledAttInGAKSsep=[];datatmpShuffledAttOutGAKSsep=[];
                        datatmpAttInShuffled = []; datatmpAttOutShuffled = [];

                        NEAttInNw = NEAttIn - 20;
                        trialSel = randperm(NEAttIn);
                        nEAttIn = nEAttIn(trialSel(1:NEAttInNw)); targetDimTrialAttIn = targetDimTrialAttIn(trialSel(1:NEAttInNw));
                        for n=1:NEAttInNw;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(targetDimTrialAttIn(n)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(targetDimTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
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
                        saveas(gcf,'HistogramISIrealCueTargetDimSCM.png')

                        % 200ms post-cue until ARRAY ONSET
                        arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1) - spikeTimes(1); arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3) - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                       
                        clear datatmp datatmpAttIn datatmpAttOut datatmpSil datatmpSilAttIn datatmpSilAttOut
                        % heatmap return plot
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss; 
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        datatmpAttIn=[];datatmpAttOut=[];
                        datatmpAttInGAKS=[];datatmpAttOutGAKS=[];
                        datatmpAttInGAKSsep=[];datatmpAttOutGAKSsep=[];
                        datatmpShuffledAttInGAKS=[];datatmpShuffledAttOutGAKS=[];
                        datatmpShuffledAttInGAKSsep=[];datatmpShuffledAttOutGAKSsep=[];
                        datatmpAttInShuffled = []; datatmpAttOutShuffled = [];

                        NEAttInNw = NEAttIn - 20;
                        trialSel = randperm(NEAttIn);
                        nEAttIn = nEAttIn(trialSel(1:NEAttInNw)); arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(trialSel(1:NEAttInNw));
                        for n=1:NEAttInNw;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttIn(n)*Fs);
                            indx=nEAttIn(n):nwinr-1;
                            if length(indx) >1
                                datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                            end
                        end
                        for n=1:NEAttOut;
                        %     nwinl=round(0.200*Fs);
                            nwinr=round(arrayOnsetTrialAttOut(n)*Fs);
                            indx=nEAttOut(n):nwinr-1;
                            if length(indx) >1
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
                        saveas(gcf,'HistogramISIrealCueArrayOnsetSCM.png')




                        figure
                        subplot(231)
                        histogram(datatmpAttInGAKS,linspace(0,.1,10000))
                        ylim([0 1400])
                        title('GAKS real data Att In')
                        subplot(234)
                        histogram(datatmpShuffledAttInGAKS,linspace(0,.1,10000))
                        ylim([0 1400])
                        title('GAKS shuffled data Att In')
                        subplot(232)
                        histogram(datatmpAttOutGAKS,linspace(0,.1,10000))
                        ylim([0 1400])
                        title('GAKS real data Att Out')
                        subplot(235)
                        histogram(datatmpShuffledAttOutGAKS,linspace(0,.1,10000))
                        ylim([0 1400])
                        title('GAKS shuffled data Att Out')
                        subplot(233)
                        histogram(datatmpGAKS,linspace(0,.1,10000))
                        ylim([0 7000])
                        title('GAKS real data all')
                        subplot(236)
                        histogram(datatmpShuffledGAKS,linspace(0,.1,10000))
                        ylim([0 7000])
                        title('GAKS shuffled data all')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSrealShuffCueTargetDim.png')

%                         figure
%                         subplot(231)
%                         histogram(datatmpAttInGAKSsep,linspace(-0.5,1.5,10000))
% %                         ylim([0 1400])
%                         title('GAKS real data Att In')
%                         subplot(234)
%                         histogram(datatmpShuffledAttInGAKSsep,linspace(-0.5,1.5,10000))
% %                         ylim([0 1400])
%                         title('GAKS shuffled data Att In')
%                         subplot(232)
%                         histogram(datatmpAttOutGAKSsep,linspace(0,.1,10000))
% %                         ylim([0 1400])
%                         title('GAKS real data Att Out')
%                         subplot(235)
%                         histogram(datatmpShuffledAttOutGAKSsep,linspace(0,.1,10000))
% %                         ylim([0 1400])
%                         title('GAKS shuffled data Att Out')
%                         subplot(233)
%                         histogram(datatmpGAKSsep,linspace(0,.1,10000))
% %                         ylim([0 7000])
%                         title('GAKS real data all')
%                         subplot(236)
%                         histogram(datatmpShuffledGAKSsep,linspace(0,.1,10000))
% %                         ylim([0 7000])
%                         title('GAKS shuffled data all')
%                         cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
%                         saveas(gcf,'HistogramGAKSThresholdedrealShuffCueTargetDim.png')
                        
                        % equate spike count for att in and att out
                        % heatmap return plot
                        targetDimTrialAttOut = UE.targetDimByLoc{1} - spikeTimes(1); targetDimTrialAttIn = UE.targetDimByLoc{3} - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                        data = binarySpikeTrain; dataGAKS = allSpikesGauss; dataGAKS2 = shuffledSpikesGauss; 
                        dataGAKSsep = allSpikesGaussSep; dataGAKSsep2 = shuffledSpikesGaussSep;
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        datatmpAttInSCM=[];datatmpAttOutSCM=[];
                        datatmpAttInShuffledSCM=[];datatmpAttOutShuffledSCM=[];
                        datatmpAttInGAKSSCM=[];datatmpAttOutGAKSSCM=[];
                        datatmpAttInGAKSsepSCM=[];datatmpAttOutGAKSsepSCM=[];
                        datatmpShuffledAttInGAKSSCM=[];datatmpShuffledAttOutGAKSSCM=[];
                        datatmpShuffledAttInGAKSsepSCM=[];datatmpShuffledAttOutGAKSsepSCM=[];
                        nTrials = max(NEAttIn, NEAttOut);
                        if NEAttIn ~= nTrials 
                            upsampledTrials = randperm(NEAttIn,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttIn = [cueOnsetTrialAttIn; cueOnsetTrialAttIn(upsampledTrials)];
                            targetDimTrialAttIn = [targetDimTrialAttIn; targetDimTrialAttIn(upsampledTrials)];
                            NEAttIn=length(cueOnsetTrialAttIn);
                        elseif NEAttOut ~= nTrials 
                            upsampledTrials = randperm(NEAttOut,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttOut = [cueOnsetTrialAttOut; cueOnsetTrialAttOut(upsampledTrials)];
                            targetDimTrialAttOut = [targetDimTrialAttOut; targetDimTrialAttOut(upsampledTrials)];
                            NEAttOut=length(cueOnsetTrialAttOut);
                        end
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        
                        for n=1:nTrials;
                        %     nwinl=round(0.200*Fs);
                            nwinrAttIn=round(targetDimTrialAttIn(n)*Fs);
                            nwinrAttOut=round(targetDimTrialAttOut(n)*Fs);
                            indxAttIn=nEAttIn(n):nwinrAttIn-1;
                            indxAttOut=nEAttOut(n):nwinrAttOut-1;
                            if length(indxAttIn) >1 && length(indxAttOut) >1
                                if length(indxAttIn) == length(indxAttOut)
                                    datatmpAttInSCM=[datatmpAttInSCM diff(find(data(indxAttIn)))];
                                    datatmpAttInGAKSSCM=[datatmpAttInGAKSSCM dataGAKS(indxAttIn)];
                                    datatmpAttInGAKSsepSCM=[datatmpAttInGAKSsepSCM dataGAKSsep(indxAttIn)];
                                    datatmpShuffledAttInGAKSSCM=[datatmpShuffledAttInGAKSSCM dataGAKS2(indxAttIn)];
                                    datatmpShuffledAttInGAKSsepSCM=[datatmpShuffledAttInGAKSsepSCM dataGAKSsep2(indxAttIn)];
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn)));
                                    datatmpAttInShuffledSCM=[datatmpAttInShuffledSCM diff(find(data(indxAttIn)))];
                                    
                                    datatmpAttOutSCM=[datatmpAttOutSCM diff(find(data(indxAttOut)))];
                                    datatmpAttOutGAKSSCM=[datatmpAttOutGAKSSCM dataGAKS(indxAttOut)];
                                    datatmpAttOutGAKSsepSCM=[datatmpAttOutGAKSsepSCM dataGAKSsep(indxAttOut)];
                                    datatmpShuffledAttOutGAKSSCM=[datatmpShuffledAttOutGAKSSCM dataGAKS2(indxAttOut)];
                                    datatmpShuffledAttOutGAKSsepSCM=[datatmpShuffledAttOutGAKSsepSCM dataGAKSsep2(indxAttOut)];
                                    indxAttOut = indxAttOut(randperm(length(indxAttOut)));
                                    datatmpAttOutShuffledSCM=[datatmpAttOutShuffledSCM diff(find(data(indxAttOut)))];
                                elseif length(indxAttIn) > length(indxAttOut)
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                    datatmpAttInSCM=[datatmpAttInSCM diff(find(data(indxAttIn)))];
                                    datatmpAttInGAKSSCM=[datatmpAttInGAKSSCM dataGAKS(indxAttIn)];
                                    datatmpAttInGAKSsepSCM=[datatmpAttInGAKSsepSCM dataGAKSsep(indxAttIn)];
                                    datatmpShuffledAttInGAKSSCM=[datatmpShuffledAttInGAKSSCM dataGAKS2(indxAttIn)];
                                    datatmpShuffledAttInGAKSsepSCM=[datatmpShuffledAttInGAKSsepSCM dataGAKSsep2(indxAttIn)];
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn)));
                                    datatmpAttInShuffledSCM=[datatmpAttInShuffledSCM diff(find(data(indxAttIn)))];
                                    
                                    datatmpAttOutSCM=[datatmpAttOutSCM diff(find(data(indxAttOut)))];
                                    datatmpAttOutGAKSSCM=[datatmpAttOutGAKSSCM dataGAKS(indxAttOut)];
                                    datatmpAttOutGAKSsepSCM=[datatmpAttOutGAKSsepSCM dataGAKSsep(indxAttOut)];
                                    datatmpShuffledAttOutGAKSSCM=[datatmpShuffledAttOutGAKSSCM dataGAKS2(indxAttOut)];
                                    datatmpShuffledAttOutGAKSsepSCM=[datatmpShuffledAttOutGAKSsepSCM dataGAKSsep2(indxAttOut)];
                                    indxAttOut = indxAttOut(randperm(length(indxAttOut)));
                                    datatmpAttOutShuffledSCM=[datatmpAttOutShuffledSCM diff(find(data(indxAttOut)))];
                                elseif length(indxAttOut) > length(indxAttIn)
                                    indxAttOut = indxAttOut(randperm(length(indxAttOut),length(indxAttIn)));
                                    datatmpAttInSCM=[datatmpAttInSCM diff(find(data(indxAttIn)))];
                                    datatmpAttInGAKSSCM=[datatmpAttInGAKSSCM dataGAKS(indxAttIn)];
                                    datatmpAttInGAKSsepSCM=[datatmpAttInGAKSsepSCM dataGAKSsep(indxAttIn)];
                                    datatmpShuffledAttInGAKSSCM=[datatmpShuffledAttInGAKSSCM dataGAKS2(indxAttIn)];
                                    datatmpShuffledAttInGAKSsepSCM=[datatmpShuffledAttInGAKSsepSCM dataGAKSsep2(indxAttIn)];
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn)));
                                    datatmpAttInShuffledSCM=[datatmpAttInShuffledSCM diff(find(data(indxAttIn)))];
                                    
                                    datatmpAttOutSCM=[datatmpAttOutSCM diff(find(data(indxAttOut)))];
                                    datatmpAttOutGAKSSCM=[datatmpAttOutGAKSSCM dataGAKS(indxAttOut)];
                                    datatmpAttOutGAKSsepSCM=[datatmpAttOutGAKSsepSCM dataGAKSsep(indxAttOut)];
                                    datatmpShuffledAttOutGAKSSCM=[datatmpShuffledAttOutGAKSSCM dataGAKS2(indxAttOut)];
                                    datatmpShuffledAttOutGAKSsepSCM=[datatmpShuffledAttOutGAKSsepSCM dataGAKSsep2(indxAttOut)];
                                    indxAttOut = indxAttOut(randperm(length(indxAttOut)));
                                    datatmpAttOutShuffledSCM=[datatmpAttOutShuffledSCM diff(find(data(indxAttOut)))];
                                end
                            end
                        end
                        figure
                        subplot(221)
                        histogram(datatmpAttInGAKSSCM,linspace(0,.1,10000))
                        ylim([0 2000])
                        title('SCM GAKS real data Att In')
                        subplot(223)
                        histogram(datatmpShuffledAttInGAKSSCM,linspace(0,.1,10000))
                        ylim([0 2000])
                        title('SCM GAKS shuffled data Att In')
                        subplot(222)
                        histogram(datatmpAttOutGAKSSCM,linspace(0,.1,10000))
                        ylim([0 2000])
                        title('SCM GAKS real data Att Out')
                        subplot(224)
                        histogram(datatmpShuffledAttOutGAKSSCM,linspace(0,.1,10000))
                        ylim([0 2000])
                        title('SCM GAKS shuffled data Att Out')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSrealShuffCueTargetDimSCM.png')

                        % Calculate 10ms firing rate for SCM data and slide
                        % 2ms + repeat GAKS with 10ms Gaussian
                        Fs = 1000;
                        targetDimTrialAttOut = UE.targetDimByLoc{1} - spikeTimes(1); targetDimTrialAttIn = UE.targetDimByLoc{3} - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                        data = binarySpikeTrain;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nTrials = max(NEAttIn, NEAttOut);
                        if NEAttIn ~= nTrials 
                            upsampledTrials = randperm(NEAttIn,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttIn = [cueOnsetTrialAttIn; cueOnsetTrialAttIn(upsampledTrials)];
                            targetDimTrialAttIn = [targetDimTrialAttIn; targetDimTrialAttIn(upsampledTrials)];
                            NEAttIn=length(cueOnsetTrialAttIn);
                        elseif NEAttOut ~= nTrials 
                            upsampledTrials = randperm(NEAttOut,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttOut = [cueOnsetTrialAttOut; cueOnsetTrialAttOut(upsampledTrials)];
                            targetDimTrialAttOut = [targetDimTrialAttOut; targetDimTrialAttOut(upsampledTrials)];
                            NEAttOut=length(cueOnsetTrialAttOut);
                        end
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        
                        attInSCM2plot = [];
                        attOutSCM2plot = [];
                        for n = 1:nTrials
                            nwinrAttIn=round(targetDimTrialAttIn(n)*Fs);
                            nwinrAttOut=round(targetDimTrialAttOut(n)*Fs);
                            indxAttIn=nEAttIn(n):nwinrAttIn-1;
                            indxAttOut=nEAttOut(n):nwinrAttOut-1;
                            if length(indxAttIn) == length(indxAttOut)
                                attInSCMall = movmean(data(indxAttIn),10);
                                attInSCM{n} = attInSCMall(1:2:length(indxAttIn));
                                attOutSCMall = movmean(data(indxAttOut),10);
                                attOutSCM{n} = attOutSCMall(1:2:length(indxAttOut));
                                attInSCM2plot = [attInSCM2plot attInSCM{n}];
                                attOutSCM2plot = [attOutSCM2plot attOutSCM{n}];
                            elseif length(indxAttIn) > length(indxAttOut)
                                indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                attInSCMall = movmean(data(indxAttIn),10);
                                attInSCM{n} = attInSCMall(1:2:length(indxAttIn));
                                attOutSCMall = movmean(data(indxAttOut),10);
                                attOutSCM{n} = attOutSCMall(1:2:length(indxAttOut));
                                attInSCM2plot = [attInSCM2plot attInSCM{n}];
                                attOutSCM2plot = [attOutSCM2plot attOutSCM{n}];
                            elseif length(indxAttIn) < length(indxAttOut)
                                indxAttOut = indxAttOut(randperm(length(indxAttOut),length(indxAttIn)));
                                attInSCMall = movmean(data(indxAttIn),10);
                                attInSCM{n} = attInSCMall(1:2:length(indxAttIn));
                                attOutSCMall = movmean(data(indxAttOut),10);
                                attOutSCM{n} = attOutSCMall(1:2:length(indxAttOut));
                                attInSCM2plot = [attInSCM2plot attInSCM{n}];
                                attOutSCM2plot = [attOutSCM2plot attOutSCM{n}];
                            end
                        end

                        figure
                        subplot(141)
                        hAttIn = histogram(attInSCM2plot,linspace(0,0.3,1000))
                        title('Att In')
                        subplot(142)
                        histogram(attInSCM2plot,linspace(0,0.3,1000))
                        ylim([0 50])
                        title('Att In Zoomed in')
                        subplot(143)
                        hAttOut = histogram(attOutSCM2plot,linspace(0,0.3,1000))
                        title('Att Out')
                        subplot(144)
                        histogram(attOutSCM2plot,linspace(0,0.3,1000))
                        ylim([0 50])
                        title('Att Out Zoomed in')
                        saveas(gcf,'movMean10msAttInAttOutSCM.png')
                        hAttIn.Values(find(hAttIn.Values>0))
                        hAttOut.Values(find(hAttOut.Values>0))

                        % create poisson distribution for the comparison
                        targetDimTrialAttOut = UE.targetDimByLoc{1} - spikeTimes(1); targetDimTrialAttIn = UE.targetDimByLoc{3} - spikeTimes(1);
                        cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200 - spikeTimes(1); cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200 - spikeTimes(1);
                        data = binarySpikeTrain; 
                        Fs = 1000;
                        NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
                        nTrials = max(NEAttIn, NEAttOut);
                        if NEAttIn ~= nTrials 
                            upsampledTrials = randperm(NEAttIn,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttIn = [cueOnsetTrialAttIn; cueOnsetTrialAttIn(upsampledTrials)];
                            targetDimTrialAttIn = [targetDimTrialAttIn; targetDimTrialAttIn(upsampledTrials)];
                            NEAttIn=length(cueOnsetTrialAttIn);
                        elseif NEAttOut ~= nTrials 
                            upsampledTrials = randperm(NEAttOut,diff([NEAttIn,NEAttOut]));
                            cueOnsetTrialAttOut = [cueOnsetTrialAttOut; cueOnsetTrialAttOut(upsampledTrials)];
                            targetDimTrialAttOut = [targetDimTrialAttOut; targetDimTrialAttOut(upsampledTrials)];
                            NEAttOut=length(cueOnsetTrialAttOut);
                        end
                        nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                        poissonBinary = []; attInBinary = []; attOutBinary = [];
                        for n=1:nTrials;
                        %     nwinl=round(0.200*Fs);
                            nwinrAttIn=round(targetDimTrialAttIn(n)*Fs);
                            nwinrAttOut=round(targetDimTrialAttOut(n)*Fs);
                            indxAttIn=nEAttIn(n):nwinrAttIn-1;
                            indxAttOut=nEAttOut(n):nwinrAttOut-1;
                            if length(indxAttIn) >1 && length(indxAttOut) >1
                                if length(indxAttIn) == length(indxAttOut)
                                    durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                    spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                    timeStepS = 0.001;                  % 1 msec
                                    times = 0:timeStepS:durationS;	% a vector with each time step		
                                    vt = rand(size(times)); % use rand to set probobility of spike
                                    spikes = (spikesPerS*timeStepS) > vt;
                                    poissonBinary = [poissonBinary spikes];
                                    attInBinary = [attInBinary data(indxAttIn)];
                                    attOutBinary = [attOutBinary data(indxAttOut)];
                                elseif length(indxAttIn) > length(indxAttOut)
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                    durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                    spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                    timeStepS = 0.001;                  % 1 msec
                                    times = 0:timeStepS:durationS;	% a vector with each time step		
                                    vt = rand(size(times)); % use rand to set probobility of spike
                                    spikes = (spikesPerS*timeStepS) > vt;
                                    poissonBinary = [poissonBinary spikes];
                                    attInBinary = [attInBinary data(indxAttIn)];
                                    attOutBinary = [attOutBinary data(indxAttOut)];
                                elseif length(indxAttOut) > length(indxAttIn)
                                    indxAttOut = indxAttOut(randperm(length(indxAttOut),length(indxAttIn)));
                                    indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                    durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                    spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                    timeStepS = 0.001;                  % 1 msec
                                    times = 0:timeStepS:durationS;	% a vector with each time step		
                                    vt = rand(size(times)); % use rand to set probobility of spike
                                    spikes = (spikesPerS*timeStepS) > vt;
                                    poissonBinary = [poissonBinary spikes];
                                    attInBinary = [attInBinary data(indxAttIn)];
                                    attOutBinary = [attOutBinary data(indxAttOut)];
                                end
                            end
                        end
                        
                        % convolve with gaussian
                        attInSpikesGauss = conv(attInBinary, gaussian,'same');
                        attOutSpikesGauss = conv(attOutBinary, gaussian,'same');
                        poissonSpikesGauss = conv(poissonBinary, gaussian,'same');
                        figure
                        subplot(231)
                        loglog([NaN diff(find(attInBinary))],[diff(find(attInBinary)) NaN],'.','MarkerSize',1)
                        title('returnplot Att In')
                        xlim([0 10000]); ylim([0 10000]);
                        subplot(232)
                        loglog([NaN diff(find(attOutBinary))],[diff(find(attOutBinary)) NaN],'.','MarkerSize',1)
                        title('returnplot Att Out')
                        xlim([0 10000]); ylim([0 10000]);
                        subplot(233)
                        loglog([NaN diff(find(poissonBinary))],[diff(find(poissonBinary)) NaN],'.','MarkerSize',1)
                        title('returnplot Poisson')
                        xlim([0 10000]); ylim([0 10000]);
                        subplot(234)
                        histogram(attInSpikesGauss,linspace(0,.1,10000))
                        title('SCM GAKS Att In')
                        ylim([0 1400])
                        subplot(235)
                        histogram(attOutSpikesGauss,linspace(0,.1,10000))
                        title('SCM GAKS Att Out')
                        ylim([0 1400])
                        subplot(236)
                        histogram(poissonSpikesGauss,linspace(0,.1,10000))
                        title('SCM GAKS Poisson')
                        ylim([0 1400])
                        saveas(gcf,'RPHistGAKSrealPoissonCueTargetDimSCM.png')

                        % create several poisson distributions
                        figure
                        for pi = 1:10
                            nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
                            poissonBinary = []; attInBinary = []; attOutBinary = [];
                            for n=1:nTrials;
                            %     nwinl=round(0.200*Fs);
                                nwinrAttIn=round(targetDimTrialAttIn(n)*Fs);
                                nwinrAttOut=round(targetDimTrialAttOut(n)*Fs);
                                indxAttIn=nEAttIn(n):nwinrAttIn-1;
                                indxAttOut=nEAttOut(n):nwinrAttOut-1;
                                if length(indxAttIn) >1 && length(indxAttOut) >1
                                    if length(indxAttIn) == length(indxAttOut)
                                        durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                        spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                        timeStepS = 0.001;                  % 1 msec
                                        times = 0:timeStepS:durationS;	% a vector with each time step		
                                        vt = rand(size(times)); % use rand to set probobility of spike
                                        spikes = (spikesPerS*timeStepS) > vt;
                                        poissonBinary = [poissonBinary spikes];
                                    elseif length(indxAttIn) > length(indxAttOut)
                                        indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                        durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                        spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                        timeStepS = 0.001;                  % 1 msec
                                        times = 0:timeStepS:durationS;	% a vector with each time step		
                                        vt = rand(size(times)); % use rand to set probobility of spike
                                        spikes = (spikesPerS*timeStepS) > vt;
                                        poissonBinary = [poissonBinary spikes];
                                    elseif length(indxAttOut) > length(indxAttIn)
                                        indxAttOut = indxAttOut(randperm(length(indxAttOut),length(indxAttIn)));
                                        indxAttIn = indxAttIn(randperm(length(indxAttIn),length(indxAttOut)));
                                        durationS = targetDimTrialAttIn(n) - cueOnsetTrialAttIn(n); % length of simulation
                                        spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
                                        timeStepS = 0.001;                  % 1 msec
                                        times = 0:timeStepS:durationS;	% a vector with each time step		
                                        vt = rand(size(times)); % use rand to set probobility of spike
                                        spikes = (spikesPerS*timeStepS) > vt;
                                        poissonBinary = [poissonBinary spikes];
                                    end
                                end
                            end

                            % convolve with gaussian
                            poissonSpikesGauss = conv(poissonBinary, gaussian,'same');
                            
                            subplot(2,10,pi)
                            loglog([NaN diff(find(poissonBinary))],[diff(find(poissonBinary)) NaN],'.','MarkerSize',1)
                            xlim([0 10000]); ylim([0 10000]);
                            subplot(2,10,pi+10)
                            histogram(poissonSpikesGauss,linspace(0,.05,1000))
                            title('SCM GAKS Poisson')
                            ylim([0 8000]);
                            clear poissonBinary poissonSpikesGauss
                        end
                        saveas(gcf,'RPHistGAKSSeveralPoissonCueTargetDimSCM.png')






                        allSpikesGaussSepCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepCueLocked = createEventLockedGAKS(allSpikesGaussSep,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussSepCueLocked.window);

                        figure
                        imagesc(allSpikesGaussSepCueLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('cue locked t=1000')
                        saveas(gcf,'cueLockedGAKS.png')

                        allSpikesGaussSepArrayLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepArrayLocked = createEventLockedGAKS(allSpikesGaussSep,UE.arrayOnset- spikeTimes(1),D.directFs,allSpikesGaussSepArrayLocked.window);

                        figure
                        imagesc(allSpikesGaussSepArrayLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('array locked t=1000')
                        saveas(gcf,'arrayLockedGAKS.png')
                        
                        % fixation locked
                        allSpikesGaussSepFixLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepFixLocked = createEventLockedGAKS(allSpikesGaussSep,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue- spikeTimes(1),D.directFs,allSpikesGaussSepFixLocked.window);

                        figure
                        imagesc(allSpikesGaussSepFixLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('fixation locked t=1000')

                        hold on
                        cueTimeReltoFix = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
                        timeOnxAxis = linspace(-1,2,size(allSpikesGaussSepFixLocked.eventLockedGAKS,2));
                        cueTimeReltoFix = (cueTimeReltoFix*1000) + 1000; cueTimeReltoFix = round(cueTimeReltoFix);
                        for ci = 1:length(cueTimeReltoFix)
                            plot([cueTimeReltoFix(ci) cueTimeReltoFix(ci)],[ci ci],'*r');%,'LineWidth',2,'MarkerSize',10)
                        end
                        saveas(gcf,'fixLockedGAKS.png')

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %           Loosen boundaries
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % cut in trials
                        sepLowHighLoose = 0.03;
                        allSpikesGaussSepTmp = ((allSpikesGauss<sepLowHighLoose & allSpikesGauss>.0015) * 0.5);
                        allSpikesGaussSep = allSpikesGaussSepTmp + (allSpikesGauss<=0.0015);

%                         allSpikesGaussSepCueLocked2.window = [1 2]; % seconds before, after
%                         allSpikesGaussSepCueLocked2 = createdatamatc(allSpikesGaussSep',UE.cueOnset,D.directFs,allSpikesGaussSepCueLocked2.window);

                        allSpikesGaussSepCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepCueLocked = createEventLockedGAKS(allSpikesGaussSep,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussSepCueLocked.window);
                        
                        allSpikesGaussCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussCueLocked = createEventLockedGAKS(allSpikesGauss,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussCueLocked.window);

                        figure
%                         subplot(211)
                        imagesc(allSpikesGaussSepCueLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('cue locked t=1000')
%                         ylim([2.5 4.5])
%                         subplot(212)
%                         hold on
%                         for ti = 3:4
%                             plot(allSpikesGaussCueLocked.eventLockedGAKS(ti,:))
%                         end
                        saveas(gcf,'cueLockedGAKSloose2.png')

                        allSpikesGaussSepArrayLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepArrayLocked = createEventLockedGAKS(allSpikesGaussSep,UE.arrayOnset- spikeTimes(1),D.directFs,allSpikesGaussSepArrayLocked.window);

                        figure
                        imagesc(allSpikesGaussSepArrayLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('array locked t=1000')
                        saveas(gcf,'arrayLockedGAKSloose2.png')
                        
                        % fixation locked
                        allSpikesGaussSepFixLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepFixLocked = createEventLockedGAKS(allSpikesGaussSep,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue- spikeTimes(1),D.directFs,allSpikesGaussSepFixLocked.window);

                        figure
                        imagesc(allSpikesGaussSepFixLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('fixation locked t=1000')
                        hold on
                        cueTimeReltoFix = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
                        timeOnxAxis = linspace(-1,2,size(allSpikesGaussSepFixLocked.eventLockedGAKS,2));
                        cueTimeReltoFix = (cueTimeReltoFix*1000) + 1000; cueTimeReltoFix = round(cueTimeReltoFix);
                        for ci = 1:length(cueTimeReltoFix)
                            plot([cueTimeReltoFix(ci) cueTimeReltoFix(ci)],[ci ci],'*r');%,'LineWidth',2,'MarkerSize',10)
                        end
                        saveas(gcf,'fixLockedGAKSloose2.png')

                        figure
                        subplot(121)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepFixLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 1),allSpikesGaussSepFixLocked.eventLockedGAKS(UE.cueLoc == 1,:))
                        colormap(gray)
                        title('Fix locked attend away')
                        subplot(122)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepFixLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 3),allSpikesGaussSepFixLocked.eventLockedGAKS(UE.cueLoc == 3,:))
                        colormap(gray)
                        title('Fix locked attend in')
                        saveas(gcf,'GAKSFixLockedAttLocs.png')

                        figure
                        subplot(121)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepCueLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 1),allSpikesGaussSepCueLocked.eventLockedGAKS(UE.cueLoc == 1,:))
                        colormap(gray)
                        title('Cue locked attend away')
                        subplot(122)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepCueLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 3),allSpikesGaussSepCueLocked.eventLockedGAKS(UE.cueLoc == 3,:))
                        colormap(gray)
                        title('Cue locked attend in')
                        saveas(gcf,'GAKSCueLockedAttLocs.png')

                        figure
                        subplot(121)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepArrayLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 1),allSpikesGaussSepArrayLocked.eventLockedGAKS(UE.cueLoc == 1,:))
                        colormap(gray)
                        title('Array locked attend away')
                        subplot(122)
                        imagesc(linspace(-1,2,size(allSpikesGaussSepArrayLocked.eventLockedGAKS,2)),1:sum(UE.cueLoc == 3),allSpikesGaussSepArrayLocked.eventLockedGAKS(UE.cueLoc == 3,:))
                        colormap(gray)
                        title('Array locked attend in')
                        saveas(gcf,'GAKSArrayLockedAttLocs.png')

                    
                        % plot fixation-onset aligned GAKS
                        allSpikesGaussFixLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussFixLocked = createEventLockedGAKS(allSpikesGauss,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue- spikeTimes(1),D.directFs,allSpikesGaussFixLocked.window);
                        
                        figure
                        subplot(211)
                        plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussFixLocked.eventLockedGAKS,1),'LineWidth',2)
                        hold on
                        for loci = 1:4
                            allSpikesGaussFixLockedLoc = [];
                            allSpikesGaussFixLockedLoc.window = [1 2]; % seconds before, after
                            allSpikesGaussFixLockedLoc = createEventLockedGAKS(allSpikesGauss,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(UE.cueLoc==(loci))- spikeTimes(1),D.directFs,allSpikesGaussFixLocked.window);
                            plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussFixLockedLoc.eventLockedGAKS,1))
                        end
                        title('Fixation locked')
                        subplot(212)
                        plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussCueLocked.eventLockedGAKS,1),'LineWidth',2)
                        hold on
                        for loci = 1:4
                            allSpikesGaussCueLockedLoc = [];
                            allSpikesGaussCueLockedLoc.window = [1 2]; % seconds before, after
                            allSpikesGaussCueLockedLoc = createEventLockedGAKS(allSpikesGauss,UE.cueOnset(UE.cueLoc==(loci))- spikeTimes(1),D.directFs,allSpikesGaussCueLocked.window);
                            plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussCueLockedLoc.eventLockedGAKS,1))
                        end
                        title('Cue locked')
                        saveas(gcf,'GAKSlineplot.png')

                        % plot histogram of gaks
                        figure
                        subplot(131)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 700000])
                        subplot(132)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 1000])
                        subplot(133)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 100])
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSzoomingIn.png')
                        
                        figure
                        subplot(121)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 700000])
                        subplot(122)
                        histogram(shuffledSpikesGauss,linspace(0,0.12,10000))
                        title('shuffled data')
                        ylim([0 700000])
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSrealShuff.png')
                        
                        


                        
                        
                        figure
                        subplot(231)
                        histogram(datatmpAttIn,linspace(0,600,100))
                        ylim([0 90])
                        title('real data Att In')
                        subplot(234)
                        histogram(datatmpAttInShuffled,linspace(0,600,100))
                        ylim([0 90])
                        title('shuffled data Att In')
                        subplot(232)
                        histogram(datatmpAttOut,linspace(0,600,100))
                        ylim([0 130])
                        title('real data Att Out')
                        subplot(235)
                        histogram(datatmpAttOutShuffled,linspace(0,600,100))
                        ylim([0 130])
                        title('shuffled data Att Out')
                        subplot(233)
                        histogram(datatmpGAKS,linspace(0,.04,10000))
                        ylim([0 5000])
                        title('real data all GAKS')
                        subplot(236)
                        histogram(datatmpShuffledGAKS,linspace(0,.04,10000))
                        ylim([0 5000])
                        title('shuffled data all GAKS')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramISIrealShuffCueTargetDim.png')
                        
                        allSpikesGaussFixLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussFixLocked = createEventLockedGAKS(allSpikesGauss,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue- spikeTimes(1),D.directFs,allSpikesGaussFixLocked.window);
                        
                        figure
                        subplot(211)
                        plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussFixLocked.eventLockedGAKS,1),'LineWidth',2)
                        hold on
                        for loci = 1:4
                            allSpikesGaussFixLockedLoc = [];
                            allSpikesGaussFixLockedLoc.window = [1 2]; % seconds before, after
                            allSpikesGaussFixLockedLoc = createEventLockedGAKS(allSpikesGauss,UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(UE.cueLoc==(loci))- spikeTimes(1),D.directFs,allSpikesGaussFixLocked.window);
                            plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussFixLockedLoc.eventLockedGAKS,1))
                        end
                        title('Fixation locked')
                        subplot(212)
                        plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussCueLocked.eventLockedGAKS,1),'LineWidth',2)
                        hold on
                        for loci = 1:4
                            allSpikesGaussCueLockedLoc = [];
                            allSpikesGaussCueLockedLoc.window = [1 2]; % seconds before, after
                            allSpikesGaussCueLockedLoc = createEventLockedGAKS(allSpikesGauss,UE.cueOnset(UE.cueLoc==(loci))- spikeTimes(1),D.directFs,allSpikesGaussCueLocked.window);
                            plot(linspace(-1,2,size(allSpikesGaussFixLocked.eventLockedGAKS,2)),mean(allSpikesGaussCueLockedLoc.eventLockedGAKS,1))
                        end
                        title('Cue locked')
                        saveas(gcf,'GAKSlineplot.png')

                        % plot histogram of gaks
                        figure
                        subplot(131)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 700000])
                        subplot(132)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 1000])
                        subplot(133)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 100])
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSzoomingIn.png')
                        
                        figure
                        subplot(121)
                        histogram(allSpikesGauss,linspace(0,0.12,10000))
                        title('real data')
                        ylim([0 700000])
                        subplot(122)
                        histogram(shuffledSpikesGauss,linspace(0,0.12,10000))
                        title('shuffled data')
                        ylim([0 700000])
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,'HistogramGAKSrealShuff.png')

% note that histogram lineplot with same edges can only be used when
% binwidth are specified and equal for both histograms
%                         [N,EDGES,BINS] = histcounts(allSpikesGauss,linspace(0,0.12,10000));
%                         shuffledBinarySpikeTrain = binarySpikeTrain(randperm(length(binarySpikeTrain)));
%                         shuffledSpikesGauss = conv(shuffledBinarySpikeTrain, gaussian);
%                         [Ns,EDGESs,BINSs] = histcounts(shuffledSpikesGauss,linspace(0,0.12,10000));
%                         figure
%                         plot(EDGES(1:end-1), log10(N))
%                         hold on
%                         plot(EDGES(1:end-1), log10(Ns))
%                         legend('real data','shuffled1','shuffled2','shuffled3','shuffled4','shuffled5')
%                         saveas(gcf,'Log histogram GAKS same edges.png')
                        
%                         figure
%                         subplot(211)
%                         plot(allSpikesGauss)
%                         hold on
%                         plot(shuffledSpikesGauss)
%                         xlim([0 1000000])
%                         ylim([0 0.04])
%                         subplot(212)
%                         plot(EDGES(1:end-1), N)
%                         hold on
%                         plot(EDGES(1:end-1), Ns)
%                         ylim([0 100000])
%                         xlim([0 0.05])
%                         saveas(gcf,'GaksSnippetAndGaksHist.png')
                        
%                         figure
%                         subplot(211)
%                         plot(EDGES(1:end-1), N)
%                         xlim([0.04 0.12])
%                         subplot(212)
%                         plot(EDGES(1:end-1), N)
%                         xlim([0.08 0.12])

                        % You can do plot(x, gaussian) to make sure this is a normal distribution with the right mean and sd. x goes from -5 sigma to 5 sigma because values should be negligible outside of that range. You can also use -4 sigma to 4 sigma. I multiply by 1000 so that x is in steps of milliseconds, which is what your binarySpikeTrain is in.
                        % You can test this with fake data:
%                         t = 0.001:0.001:1; % one second of data
%                         binarySpikeTrainfk = zeros(numel(t), 1); % one second of data
%                         binarySpikeTrainfk(100) = 1; % set a spike at t = 0.1
%                         binarySpikeTrainfk(101) = 1; % set a spike at t = 0.9
%                         allSpikesGaussfk = conv(binarySpikeTrainfk, gaussian, 'same');  % the 'same' argument restricts the result to the span of binarySpikeTrain
%                         figure;plot(t, allSpikesGaussfk)
                        
                        % create heatmap of convolved data
                        allSpikesGauss = [allSpikesGauss zeros(1,623*10000-size(allSpikesGauss,2))];
                        spikesGaussFakeTrials = reshape(allSpikesGauss,10000,size(allSpikesGauss,2)/10000);
                        figure
                        imagesc(spikesGaussFakeTrials',[-0.05 0.05])
                        colorbar
                        title('real data')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,['Fig1_burstInfo_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        spikesGaussFT4plotWsepTmp = ((spikesGaussFakeTrials<sepLowHigh) * 0.5);
                        spikesGaussFT4plotWsep = spikesGaussFT4plotWsepTmp + (spikesGaussFakeTrials==0);
                        
                        figure
                        imagesc(spikesGaussFT4plotWsep')
                        colormap(gray)
                        saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity 
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % create similar heatmap for cue and array onset
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
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
                        spikeTimesTask = D.allUnitStructs{uniti}.ts;
                        gaussCueOnset.window = [1 2]; % seconds before, after
                        gaussCueOnset.spdfWindowOffset = [-0.9 1.9]; % tighter window for spdf to avoid edge effects
                        kernelSigma = 0.100;
                        gaussCueOnset = createTimeLockedSpdf(spikeTimesTask, UE.cueOnset, UE.cueOnsetByLoc, gaussCueOnset, kernelSigma, startTime, endTime);
                        
                        spikesGaussCueLockedFT4plotWsepTmp = ((gaussCueOnset.singleTrialSpdf<sepLowHigh) * 0.5);
                        spikesGaussCueLockedFT4plotWsep = spikesGaussCueLockedFT4plotWsepTmp + (gaussCueOnset.singleTrialSpdf==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedFT4plotWsep)
                        caxis([0 1])
                        colormap(gray)
                        title('cue locked t=91')
                        
                        gaussArrayOnset.window = [1 2]; % seconds before, after
                        gaussArrayOnset.spdfWindowOffset = [-0.9 1.9]; % tighter window for spdf to avoid edge effects
                        kernelSigma = 0.100;
                        gaussArrayOnset = createTimeLockedSpdf(spikeTimesTask, UE.arrayOnset, UE.arrayOnsetByLoc, gaussArrayOnset, kernelSigma, startTime, endTime);
                        
                        spikesGaussArrayLockedFT4plotWsepTmp = ((gaussArrayOnset.singleTrialSpdf<sepLowHigh) * 0.5);
                        spikesGaussArrayLockedFT4plotWsep = spikesGaussArrayLockedFT4plotWsepTmp + (gaussArrayOnset.singleTrialSpdf==0);
                        
                        figure
                        contour(linspace(-0.9,1.9,length(gaussCueOnset.t)),...
                            1:size(spikesGaussArrayLockedFT4plotWsep,1),...
                            spikesGaussArrayLockedFT4plotWsep,'LineWidth',0.1)
                        caxis([0 1])
                        colormap(gray)
                        title('array locked t=91')
                        
                        
                        
                        
                        % for cue locked data
                        spikeTimesTask = D.allUnitStructs{uniti}.ts;
                        windowOI = [1 2];
                        [spikeTimesCueLocked,spikeTimesCueLockedInd] = createnonemptydatamatpt(spikeTimesTask, UE.cueOnset, windowOI);
                        [spikeTimesCueLockedAttOut,spikeTimesCueLockedAttOutInd] = createnonemptydatamatpt(spikeTimesTask, UE.cueOnset(UE.cueLoc == 1), windowOI);
                        [spikeTimesCueLockedAttIn,spikeTimesCueLockedAttInInd] = createnonemptydatamatpt(spikeTimesTask, UE.cueOnset(UE.cueLoc == 3), windowOI);
                        timeOI = 0:.001:3;
                        
                        % for all cue locked trials
                        binaryCueLocked = zeros(size(spikeTimesCueLocked,2),length(timeOI));
                        cueLockedSpikes = 0;
                        gaussCueLocked = zeros(size(spikeTimesCueLocked,2),(length(gaussian) +length(timeOI))-1);
                        for triali = 1:size(spikeTimesCueLocked,2)
                            binaryCueLocked(triali,ismembertol(timeOI,round(spikeTimesCueLocked(triali).times,3),.00000001)) = 1; 
                            cueLockedSpikes = cueLockedSpikes + length(spikeTimesCueLocked(triali).times);
                            gaussCueLocked(triali,:) = conv(binaryCueLocked(triali,:), gaussian);
                        end
                        spikesGaussCueLockedFT4plotWsepTmp = ((gaussCueLocked<sepLowHigh) * 0.5);
                        spikesGaussCueLockedFT4plotWsep = spikesGaussCueLockedFT4plotWsepTmp + (gaussCueLocked==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedFT4plotWsep)
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        % for att away cue locked trials
                        binaryCueLockedAttOut = zeros(size(spikeTimesCueLockedAttOut,2),length(timeOI));
                        cueLockedAttOutSpikes = 0;
                        gaussCueLockedAttOut = zeros(size(spikeTimesCueLockedAttOut,2),(length(gaussian) * 2)-1);
                        for triali = 1:size(spikeTimesCueLockedAttOut,2)
                            binaryCueLockedAttOut(triali,ismembertol(timeOI,round(spikeTimesCueLockedAttOut(triali).times,3),.00000001)) = 1; 
                            cueLockedAttOutSpikes = cueLockedAttOutSpikes + length(spikeTimesCueLockedAttOut(triali).times);
                            gaussCueLockedAttOut(triali,:) = conv(binaryCueLockedAttOut(triali,:), gaussian);
                        end
                        spikesGaussCueLockedAttOutFT4plotWsepTmp = ((gaussCueLockedAttOut<sepLowHigh) * 0.5);
                        spikesGaussCueLockedAttOutFT4plotWsep = spikesGaussCueLockedAttOutFT4plotWsepTmp + (gaussCueLockedAttOut==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedAttOutFT4plotWsep)
                        caxis([0 1])
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        % for att in cue locked trials
                        binaryCueLockedAttIn = zeros(size(spikeTimesCueLockedAttIn,2),length(timeOI));
                        cueLockedAttInSpikes = 0;
                        gaussCueLockedAttIn = zeros(size(spikeTimesCueLockedAttIn,2),(length(gaussian) * 2)-1);
                        for triali = 1:size(spikeTimesCueLockedAttIn,2)
                            binaryCueLockedAttIn(triali,ismembertol(timeOI,round(spikeTimesCueLockedAttIn(triali).times,3),.00000001)) = 1; 
                            cueLockedAttInSpikes = cueLockedAttInSpikes + length(spikeTimesCueLockedAttIn(triali).times);
                            gaussCueLockedAttIn(triali,:) = conv(binaryCueLockedAttIn(triali,:), gaussian);
                        end
                        spikesGaussCueLockedAttInFT4plotWsepTmp = ((gaussCueLockedAttIn<sepLowHigh) * 0.5);
                        spikesGaussCueLockedAttInFT4plotWsep = spikesGaussCueLockedAttInFT4plotWsepTmp + (gaussCueLockedAttIn==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedAttInFT4plotWsep)
                        caxis([0 1])
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        % for array locked data
                        [spikeTimesArrayLocked,spikeTimesArrayLockedInd] = createnonemptydatamatpt(spikeTimesTask, UE.arrayOnset, windowOI);
                        [spikeTimesArrayLockedAttOut,spikeTimesArrayLockedAttOutInd] = createnonemptydatamatpt(spikeTimesTask, UE.arrayOnset(UE.cueLoc == 1), windowOI);
                        [spikeTimesArrayLockedAttIn,spikeTimesArrayLockedAttInInd] = createnonemptydatamatpt(spikeTimesTask, UE.arrayOnset(UE.cueLoc == 3), windowOI);
                        
                        % for all array locked trials
                        binaryArrayLocked = zeros(size(spikeTimesArrayLocked,2),length(timeOI));
                        ArrayLockedSpikes = 0;
                        gaussArrayLocked = zeros(size(spikeTimesArrayLocked,2),(length(gaussian) +length(timeOI))-1);
                        for triali = 1:size(spikeTimesArrayLocked,2)
                            binaryArrayLocked(triali,ismembertol(timeOI,round(spikeTimesArrayLocked(triali).times,3),.00000001)) = 1; 
                            ArrayLockedSpikes = ArrayLockedSpikes + length(spikeTimesArrayLocked(triali).times);
                            gaussArrayLocked(triali,:) = conv(binaryArrayLocked(triali,:), gaussian);
                        end
                        spikesGaussArrayLockedFT4plotWsepTmp = ((gaussArrayLocked<sepLowHigh) * 0.5);
                        spikesGaussArrayLockedFT4plotWsep = spikesGaussArrayLockedFT4plotWsepTmp + (gaussArrayLocked==0);
                        
                        figure
                        imagesc(spikesGaussArrayLockedFT4plotWsep)
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        % for att away cue locked trials
                        binaryCueLockedAttOut = zeros(size(spikeTimesCueLockedAttOut,2),length(timeOI));
                        cueLockedAttOutSpikes = 0;
                        gaussCueLockedAttOut = zeros(size(spikeTimesCueLockedAttOut,2),(length(gaussian) * 2)-1);
                        for triali = 1:size(spikeTimesCueLockedAttOut,2)
                            binaryCueLockedAttOut(triali,ismembertol(timeOI,round(spikeTimesCueLockedAttOut(triali).times,3),.00000001)) = 1; 
                            cueLockedAttOutSpikes = cueLockedAttOutSpikes + length(spikeTimesCueLockedAttOut(triali).times);
                            gaussCueLockedAttOut(triali,:) = conv(binaryCueLockedAttOut(triali,:), gaussian);
                        end
                        spikesGaussCueLockedAttOutFT4plotWsepTmp = ((gaussCueLockedAttOut<sepLowHigh) * 0.5);
                        spikesGaussCueLockedAttOutFT4plotWsep = spikesGaussCueLockedAttOutFT4plotWsepTmp + (gaussCueLockedAttOut==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedAttOutFT4plotWsep)
                        caxis([0 1])
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        % for att away cue locked trials
                        binaryCueLockedAttIn = zeros(size(spikeTimesCueLockedAttIn,2),length(timeOI));
                        cueLockedAttInSpikes = 0;
                        gaussCueLockedAttIn = zeros(size(spikeTimesCueLockedAttIn,2),(length(gaussian) * 2)-1);
                        for triali = 1:size(spikeTimesCueLockedAttIn,2)
                            binaryCueLockedAttIn(triali,ismembertol(timeOI,round(spikeTimesCueLockedAttIn(triali).times,3),.00000001)) = 1; 
                            cueLockedAttInSpikes = cueLockedAttInSpikes + length(spikeTimesCueLockedAttIn(triali).times);
                            gaussCueLockedAttIn(triali,:) = conv(binaryCueLockedAttIn(triali,:), gaussian);
                        end
                        spikesGaussCueLockedAttInFT4plotWsepTmp = ((gaussCueLockedAttIn<sepLowHigh) * 0.5);
                        spikesGaussCueLockedAttInFT4plotWsep = spikesGaussCueLockedAttInFT4plotWsepTmp + (gaussCueLockedAttIn==0);
                        
                        figure
                        imagesc(spikesGaussCueLockedAttInFT4plotWsep)
                        caxis([0 1])
                        colormap(gray)
%                         saveas(gcf,'sepOnGAKShist.png')
                        % white is silence; black is bursty; gray is normal activity
                        
                        
                        
                        
                        

                        
                        [~,idx] = sort(mean(spikesGaussFakeTrials,1));
                        figure
                        imagesc(spikesGaussFakeTrials(:,idx)',[-0.05 0.05])
                        colorbar
                        title('sorted real data')
                        saveas(gcf,['Fig2_sortedBurstInfo_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        spikesGauss2plot = nan(100,10000);
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,1:50)';
                        clim = [-0.05 0.05];
                        figure;
                        subplot(171)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,51:100)';
                        subplot(172)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,101:150)';
                        subplot(173)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,151:200)';
                        subplot(174)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,201:250)';
                        subplot(175)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,251:300)';
                        subplot(176)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,301:350)';
                        subplot(177)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,351:400)';
                        title('real data with bands')
                        saveas(gcf,['Fig3_BIwithBands_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        figure;
                        subplot(171)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,401:450)';
                        subplot(172)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,451:500)';
                        subplot(173)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,501:550)';
                        subplot(174)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:99,:) = spikesGaussFakeTrials(:,551:600)';
                        subplot(175)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        spikesGauss2plot(1:2:45,:) = spikesGaussFakeTrials(:,601:end)';
                        subplot(176)
                        imagesc(spikesGauss2plot(1:100,:),clim)
                        colormap( [1 1 1; parula(256)] )
                        title('real data with bands II')
                        saveas(gcf,['Fig3b_BIwithBands_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        % create shuffled data to compare
                        shuffledBinarySpikeTrain = binarySpikeTrain(randperm(length(binarySpikeTrain)));
                        shuffledSpikesGauss = conv(shuffledBinarySpikeTrain, gaussian);
                        
                        figure
                        plot(allSpikesGauss)
                        hold on
                        plot(shuffledSpikesGauss)
                        legend('real','shuffled')
                        title('real and shuffled "PSTH"')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,['Fig4_shuffledAndBurst_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        shuffledSpikesGauss = [shuffledSpikesGauss zeros(1,623*10000-size(shuffledSpikesGauss,2))];
                        shuffledSpikesGaussFakeTrials = reshape(shuffledSpikesGauss,10000,size(shuffledSpikesGauss,2)/10000);
                        figure
                        imagesc(shuffledSpikesGaussFakeTrials',[-0.05 0.05])
                        colorbar
                        title('shuffled data')
                        cd('/Users/labmanager/Documents/MATLAB/BurstSep4')
                        saveas(gcf,['Fig1_shuffled_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        [~,idx] = sort(mean(shuffledSpikesGaussFakeTrials,1));
                        figure
                        imagesc(shuffledSpikesGaussFakeTrials(:,idx)',[-0.05 0.05])
                        colorbar
                        title('sorted shuffled data')
                        saveas(gcf,['Fig2_shuffled_' unitStruct.name '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti) '.png'])
                        close all
                        
                        
                        % ADD PPC
                        % move to fieldtrip format
                        
                        
                        
                        % save and close figures
                        
                        
                        
                        % collect info
                        AllBurstEvents = [length(burstEventArrayRelated) length(burstEventCueRelated) ...
                            length(burstEventEnterFix) length(burstEventExitFix) length(burstEventLeverPress) ...
                            length(burstEventLeverRelease) length(burstEventTargetDim) length(burstEventPostTargetDim) ...
                            length(burstEventUnrelated) length(burstEventArrayRelatedAttAway) length(burstEventArrayRelatedAttIn) ...
                            length(burstEventCueRelatedAttAway) length(burstEventCueRelatedAttIn)];
                        cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                        save(['burstInfo' sessionInfo{2}{sessioni}(2:end-1) '_sessioni' num2str(sessioni) ...
                            '_uniti' num2str(uniti)],'burstVect','burstEventOnset','spikeTimes','UE','AllBurstEvents')
                        
                        clear burstVect burstEventOnset spikeTimes AllBurstEvents
                    else
                        counterFR(sessioni) = 1;
                    end
                end
            end
        end
    end
end

