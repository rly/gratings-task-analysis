close all

clear
outputDir = '/Users/labmanager/Documents/MATLAB/BurstSep4/';
v = 14;

fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigMcCartney.csv.bak');
% fid = fopen('pulvinarRecordingInfoByProbeAndArrayConfigFerdy.csv.bak');
sessionInfo = textscan(fid, '%d8%s%s', 'Delimiter', ',', 'HeaderLines' ,1-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fid);
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
        
        for uniti = 2%1:nUnits
            [~,iidx] = ismember(sessioni,burstsession);
            if ismember(uniti,burstunit(iidx,:))
                
                unitStruct = D.allUnitStructs{uniti};
                unitName = unitStruct.name;
                fprintf('-----------------------------\n');
                fprintf('\n');
                if ~isempty(unitStruct.ts)
                    
                    spikeTimes = D.allUnitStructs{uniti}.ts;
                    totalTimeOverall = spikeTimes(end) - spikeTimes(1);
                    firingRateOverall = numel(spikeTimes) / totalTimeOverall;
                    fprintf('SPNA %s (%d/%d = %d%%)... \n', unitName, uniti, ...
                        nUnits, round(uniti/nUnits*100));
                    
                    nUnitsPerSession(sessioni) = nUnits;
                    
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
                        
                        % calculate the boundary with shuffled spikes
                        for shuffi = 1:10;
                            shuffledBinarySpikeTrain = binarySpikeTrain(randperm(length(binarySpikeTrain)));
                            shuffledSpikesGauss = conv(shuffledBinarySpikeTrain, gaussian,'same');
                            highestShuf(shuffi) = max(shuffledSpikesGauss);
                        end
                        sepLowHigh = mean(highestShuf);

                        % threshold GAKS based on shuffled data and cut in trials
                        allSpikesGaussSepTmp = ((allSpikesGauss<sepLowHigh) * 0.5);
                        allSpikesGaussSep = allSpikesGaussSepTmp + (allSpikesGauss==0);
                        
                        % cue locked
                        allSpikesGaussSepCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepCueLocked = createEventLockedGAKS(allSpikesGaussSep,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussSepCueLocked.window);
                        figure
                        imagesc(allSpikesGaussSepCueLocked.eventLockedGAKS)
                        caxis([0 1])
                        colormap(gray)
                        title('cue locked t=1000')
                        saveas(gcf,'cueLockedGAKS.png')

                        % array locked
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

                        % cue locked
                        allSpikesGaussSepCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussSepCueLocked = createEventLockedGAKS(allSpikesGaussSep,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussSepCueLocked.window);
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

                        % array locked
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
                        allSpikesGaussCueLocked.window = [1 2]; % seconds before, after
                        allSpikesGaussCueLocked = createEventLockedGAKS(allSpikesGauss,UE.cueOnset- spikeTimes(1),D.directFs,allSpikesGaussCueLocked.window);

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
                    end
                end
            end
        end
    end
end

