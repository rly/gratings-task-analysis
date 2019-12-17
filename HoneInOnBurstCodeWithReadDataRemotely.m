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
        
        for uniti = 2%1:nUnits
            [~,iidx] = ismember(sessioni,burstsession);
            if ismember(uniti,burstunit(iidx,:))
                
                unitStruct = D.allUnitStructs{uniti};
                unitName = unitStruct.name;
                fprintf('-----------------------------\n');
                fprintf('\n');
                
                spikeTimes = D.allUnitStructs{uniti}.ts;
                totalTimeOverall = spikeTimes(end) - spikeTimes(1);
                
                %                     spikeTimes(50000:end) = spikeTimes(50000:end) + 10000;
                
                firingRateOverall = numel(spikeTimes) / totalTimeOverall;
                fprintf('SPNA %s (%d/%d = %d%%)... \n', unitName, uniti, ...
                    nUnits, round(uniti/nUnits*100));
                
                nUnitsPerSession(sessioni) = nUnits;
                
                %     saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                %             processedDataDir, unitName, blockName, v);
                fprintf('\tOverall firing rate = %0.2f Hz > minimum firing rate = %0.2f Hz in these blocks.\n', ...
                    firingRateOverall, minFiringRateOverall);
                %         fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);
                
                allDiffSpikeTimes = diff(spikeTimes); % get the ISIs
                allDiffSpikeTimesPre = [NaN; allDiffSpikeTimes]; % miss align ISIs by adding NaN to create pre and postISI
                allDiffSpikeTimesPost = [allDiffSpikeTimes; NaN];
                allDiffSpikeTimesPre = log10(allDiffSpikeTimesPre*1000); % logISI in ms
                allDiffSpikeTimesPost = log10(allDiffSpikeTimesPost*1000);
                
                silencethreshold = [100 50];
                isiburstthreshold = [5 20];
                
                for threshi = 1:length(silencethreshold)
                    figure
                    % return plot
                    plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','LineWidth',1)
                    preTimesOverThreshLong = find(allDiffSpikeTimesPre > log10(silencethreshold(threshi)));
                    postTimesUnderThreshShort = find(allDiffSpikeTimesPost < log10(isiburstthreshold(threshi)));
                    hold on
                    y1 = allDiffSpikeTimesPost(preTimesOverThreshLong);
                    x2 = allDiffSpikeTimesPre(postTimesUnderThreshShort);
                    bothTimesOverThreshPre = find(allDiffSpikeTimesPre(postTimesUnderThreshShort) > log10(silencethreshold(threshi)));
                    bothTimesOverThreshPost = find(allDiffSpikeTimesPost(preTimesOverThreshLong) < log10(isiburstthreshold(threshi)));
                    x3 = x2(bothTimesOverThreshPre);
                    y3 = y1(bothTimesOverThreshPost);
                    plot(x3,y3,'k.','LineWidth', 2)
                    xlim([-1 3])
                    ylim([-1 3])
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig1_burstInfo' unitStruct.name '_s' num2str(sessioni) ...
                        '_u' num2str(uniti) 'silT' num2str(silencethreshold(threshi)) ...
                        '_isiT' num2str(isiburstthreshold(threshi)) '.png'])
                    close all
                    
                    silence_period = log10(silencethreshold(threshi));
                    min_spikes = 3;
                    min_isi = log10(isiburstthreshold(threshi));
                    i = 1;
                    burstVect = zeros(1, length(spikeTimes));
                    while i < length(spikeTimes) - 1
                        if allDiffSpikeTimesPre(i) > silence_period
                            possible_burst = true; n_sib = 0;
                            while possible_burst
                                if allDiffSpikeTimesPost(i+n_sib) < min_isi
                                else
                                    possible_burst = false;
                                end
                                n_sib = n_sib + 1;
                            end
                            
                            if n_sib >= min_spikes
                                burstVect(i) = 2;                   % 2: start of burst
                                burstVect(i+n_sib-1) = 3;           % 3: end of burst
                                if min_spikes > 2
                                    burstVect(i+1:i+n_sib-2) = 1;     % 1: within burst
                                else
                                    if n_sib-1 > 2
                                        burstVect(i+1:i+n_sib-2) = 1; % 1: within burst
                                    end
                                end
                            else
                                burstVect(i:i+n_sib) = 0;         % 0: not in burst
                            end
                            i = i+n_sib;
                        end
                        i = i + 1;
                    end
                    
                    burstEventOnset = spikeTimes(burstVect==2);
                    burstEventArrayRelated = []; burstEventCueRelated = []; burstEventEnterFix = [];
                    burstEventExitFix = []; burstEventLeverPress =[]; burstEventLeverRelease =[];
                    burstEventTargetDim = []; burstEventPostTargetDim = [];
                    burstEventUnrelated = [];
                    burstEventArrayRelatedAttAway = []; burstEventArrayRelatedAttIn = [];
                    burstEventCueRelatedAttAway = []; burstEventCueRelatedAttIn = [];
                    for bursti = 1:length(burstEventOnset)
                        for triali = 1:length(UE.arrayOnset)
                            % check whether burst event falls within 200ms
                            % of arrayOnset
                            if burstEventOnset(bursti) >= UE.arrayOnset(triali) && burstEventOnset(bursti) < ( UE.arrayOnset(triali) + 0.2 )
                                burstEventArrayRelated = [burstEventArrayRelated bursti];
                                if UE.cueLoc(triali) == 1
                                    burstEventArrayRelatedAttAway = [burstEventArrayRelatedAttAway bursti];
                                elseif UE.cueLoc(triali) == 3
                                    burstEventArrayRelatedAttIn = [burstEventArrayRelatedAttIn bursti];
                                end
                            elseif burstEventOnset(bursti) >= UE.cueOnset(triali) && burstEventOnset(bursti) < ( UE.cueOnset(triali) + 0.2 )
                                burstEventCueRelated = [burstEventCueRelated bursti];
                                if UE.cueLoc(triali) == 1
                                    burstEventCueRelatedAttAway = [burstEventCueRelatedAttAway bursti];
                                elseif UE.cueLoc(triali) == 3
                                    burstEventCueRelatedAttIn = [burstEventCueRelatedAttIn bursti];
                                end
                            elseif burstEventOnset(bursti) >= UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(triali) && burstEventOnset(bursti) < ( UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(triali) + 0.2 )
                                burstEventEnterFix = [burstEventEnterFix bursti];
                            elseif burstEventOnset(bursti) >= UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(triali) && burstEventOnset(bursti) < ( UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(triali) + 0.2 )
                                burstEventExitFix = [burstEventExitFix bursti];
                            elseif burstEventOnset(bursti) >= UE.fixationAndLeverTimes.firstLeverPressTimesPreCue(triali) && burstEventOnset(bursti) < ( UE.fixationAndLeverTimes.firstLeverPressTimesPreCue(triali) + 0.2 )
                                burstEventLeverPress = [burstEventLeverPress bursti];
                            elseif burstEventOnset(bursti) >= UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(triali) && burstEventOnset(bursti) < ( UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(triali) + 0.2 )
                                burstEventLeverRelease = [burstEventLeverRelease bursti];
                            elseif triali < length(UE.targetDim)
                                if burstEventOnset(bursti) >= (UE.targetDim(triali) - 0.2) && burstEventOnset(bursti) < UE.targetDim(triali)
                                    burstEventTargetDim = [burstEventTargetDim bursti];
                                elseif burstEventOnset(bursti) >= UE.targetDim(triali) && burstEventOnset(bursti) < (UE.targetDim(triali) + 0.2)
                                    burstEventPostTargetDim = [burstEventPostTargetDim bursti];
                                end
                            end
                        end
                    end
                    
                    sumlength = sum([length(burstEventArrayRelated) length(burstEventCueRelated) ...
                        length(burstEventEnterFix) length(burstEventExitFix) length(burstEventLeverPress) ...
                        length(burstEventLeverRelease) length(burstEventTargetDim) length(burstEventPostTargetDim)] );
                    figure
                    bar([1 2 3 4 5 6 7 8 9],[length(burstEventArrayRelated) length(burstEventCueRelated) ...
                        length(burstEventEnterFix) length(burstEventExitFix) length(burstEventLeverPress) ...
                        length(burstEventLeverRelease) length(burstEventTargetDim) length(burstEventPostTargetDim)...
                        (length(burstEventOnset)-sumlength)])
                    xticklabels({'Array','Cue','EnterFix','ExitFix','LeverPress','LeverRel','TargetDim','PostTarget','Unrelated'})
                    title('Burst events related to task events')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig2_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    figure
                    bar([1 2 3 4],[length(burstEventCueRelatedAttAway) length(burstEventCueRelatedAttIn)...
                        length(burstEventArrayRelatedAttAway) length(burstEventArrayRelatedAttIn)])
                    xticklabels({'Cue Att Away','Cue Att In','Array Att Away','Array Att In'})
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig3_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    arrayOnset.window = [0.2 0.8]; % seconds before, after
                    arrayOnset.spdfWindowOffset = [-0.1 0.7]; % tighter window for spdf to avoid edge effects
                    kernelSigma = 0.01;
                    arrayOnset = createTimeLockedSpdf(spikeTimes, UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma, spikeTimes(1), spikeTimes(end));
                    
                    cueOnset.window = [0.2 0.8]; % seconds before, after
                    cueOnset.spdfWindowOffset = [-0.1 0.7]; % tighter window for spdf to avoid edge effects
                    kernelSigma = 0.01;
                    cueOnset = createTimeLockedSpdf(spikeTimes, UE.cueOnset, UE.cueOnsetByLoc, cueOnset, kernelSigma, spikeTimes(1), spikeTimes(end));
                    
                    
                    
                    %             spikeTimesTask = spikeTimesTask(spikeTimesTask >= startTime & spikeTimesTask <= endTime);
                    [spikeTimesCueLocked,spikeTimesCueLockedInd] = createnonemptydatamatpt(spikeTimes, UE.cueOnset, arrayOnset.window);
                    [spikeTimesCueLockedAttOut,spikeTimesCueLockedAttOutInd] = createnonemptydatamatpt(spikeTimes, UE.cueOnset(UE.cueLoc == 1), arrayOnset.window);
                    [spikeTimesCueLockedAttIn,spikeTimesCueLockedAttInInd] = createnonemptydatamatpt(spikeTimes, UE.cueOnset(UE.cueLoc == 3), arrayOnset.window);
                    [spikeTimesArrayLocked,spikeTimesArrayLockedInd] = createnonemptydatamatpt(spikeTimes, UE.arrayOnset, arrayOnset.window);
                    [spikeTimesArrayLockedAttOut,spikeTimesArrayLockedAttOutInd] = createnonemptydatamatpt(spikeTimes, UE.arrayOnset(UE.cueLoc == 1), arrayOnset.window);
                    [spikeTimesArrayLockedAttIn,spikeTimesArrayLockedAttInInd] = createnonemptydatamatpt(spikeTimes, UE.arrayOnset(UE.cueLoc == 3), arrayOnset.window);
                    
                    figure;
                    subplot(211);
                    plot(arrayOnset.t,cueOnset.spdfByLoc(1,:),'b')
                    hold on
                    plot(arrayOnset.t,cueOnset.spdfByLoc(3,:),'r')
                    subplot(212);
                    plot(arrayOnset.t,arrayOnset.spdfByLoc(1,:),'b')
                    hold on
                    plot(arrayOnset.t,arrayOnset.spdfByLoc(3,:),'r')
                    legend('Att Away','Att In')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig4_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    % raster plots with bursts color coded
                    % all trials
                    % cue onset
                    figure
                    hold on
                    spikeIdx = 0;
                    for triali = 1:length(cueOnset.spikeCount)
                        for spikei = 1:length(cueOnset.spikeTimes(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(cueOnset.spikeIndices(spikeIdx)) == 0
                                plot([cueOnset.spikeTimes(triali).times(spikei) cueOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([cueOnset.spikeTimes(triali).times(spikei) cueOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([cueOnset.spikeTimes(triali).times(spikei) cueOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([cueOnset.spikeTimes(triali).times(spikei) cueOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    %                         xlim([0 0.3])
                    title('cueOnset all')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig5_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    % all trials
                    % array onset
                    figure
                    hold on
                    spikeIdx = 0;
                    %                         AOidx = 0;
                    for triali = 1:length(arrayOnset.spikeCount)
                        for spikei = 1:length(arrayOnset.spikeTimes(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(arrayOnset.spikeIndices(spikeIdx)) == 0
                                plot([arrayOnset.spikeTimes(triali).times(spikei) arrayOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([arrayOnset.spikeTimes(triali).times(spikei) arrayOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                                %                                     AOidx = AOidx + 1;
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([arrayOnset.spikeTimes(triali).times(spikei) arrayOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([arrayOnset.spikeTimes(triali).times(spikei) arrayOnset.spikeTimes(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    title('arrayOnset all')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig6_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    
                    % split on attend in and attend away conditions
                    
                    % attend in trials
                    % cue onset
                    figure
                    hold on
                    spikeIdx = 0;
                    for triali = 1:size(cueOnset.spikeTimesByLoc{3},2)
                        % Att In
                        for spikei = 1:length(cueOnset.spikeTimesByLoc{3}(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(cueOnset.spikeIndices(spikeIdx)) == 0
                                plot([cueOnset.spikeTimesByLoc{3}(triali).times(spikei) cueOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([cueOnset.spikeTimesByLoc{3}(triali).times(spikei) cueOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([cueOnset.spikeTimesByLoc{3}(triali).times(spikei) cueOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([cueOnset.spikeTimesByLoc{3}(triali).times(spikei) cueOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    %                         xlim([0 0.3])
                    title('cueOnset Att In')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig7_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    % attend away trials
                    % cue onset
                    figure
                    hold on
                    spikeIdx = 0;
                    for triali = 1:size(cueOnset.spikeTimesByLoc{1},2)
                        % Att Away
                        for spikei = 1:length(cueOnset.spikeTimesByLoc{1}(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(cueOnset.spikeIndices(spikeIdx)) == 0
                                plot([cueOnset.spikeTimesByLoc{1}(triali).times(spikei) cueOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([cueOnset.spikeTimesByLoc{1}(triali).times(spikei) cueOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([cueOnset.spikeTimesByLoc{1}(triali).times(spikei) cueOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(cueOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([cueOnset.spikeTimesByLoc{1}(triali).times(spikei) cueOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    %                         xlim([0 0.3])
                    title('cueOnset Att Away')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig8_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    figure
                    hold on
                    spikeIdx = 0;
                    %                         AOidx = 0;
                    for triali = 1:size(arrayOnset.spikeTimesByLoc{3},2)
                        % Att In
                        for spikei = 1:length(arrayOnset.spikeTimesByLoc{3}(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(arrayOnset.spikeIndices(spikeIdx)) == 0
                                plot([arrayOnset.spikeTimesByLoc{3}(triali).times(spikei) arrayOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([arrayOnset.spikeTimesByLoc{3}(triali).times(spikei) arrayOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                                %                                     AOidx = AOidx + 1;
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([arrayOnset.spikeTimesByLoc{3}(triali).times(spikei) arrayOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([arrayOnset.spikeTimesByLoc{3}(triali).times(spikei) arrayOnset.spikeTimesByLoc{3}(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    title('arrayOnset Att In')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig9_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
                        '_uniti' num2str(uniti) '.png'])
                    close all
                    
                    figure
                    hold on
                    spikeIdx = 0;
                    %                         AOidx = 0;
                    for triali = 1:size(arrayOnset.spikeTimesByLoc{1},2)
                        % Att Away
                        for spikei = 1:length(arrayOnset.spikeTimesByLoc{1}(triali).times)
                            spikeIdx = spikeIdx + 1;
                            if burstVect(arrayOnset.spikeIndices(spikeIdx)) == 0
                                plot([arrayOnset.spikeTimesByLoc{1}(triali).times(spikei) arrayOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'k')
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 2 % start of burst
                                plot([arrayOnset.spikeTimesByLoc{1}(triali).times(spikei) arrayOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'r','LineWidth',4)
                                %                                     AOidx = AOidx + 1;
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 3 % end of burst
                                plot([arrayOnset.spikeTimesByLoc{1}(triali).times(spikei) arrayOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'b','LineWidth',4)
                            elseif burstVect(arrayOnset.spikeIndices(spikeIdx)) == 1 % within burst
                                plot([arrayOnset.spikeTimesByLoc{1}(triali).times(spikei) arrayOnset.spikeTimesByLoc{1}(triali).times(spikei)],...
                                    [triali triali+1],'m','LineWidth',4)
                            end
                        end
                    end
                    title('arrayOnset Att Away')
                    
                    cd('/Users/labmanager/Documents/MATLAB/burstEvents')
                    saveas(gcf,['Fig10_burstInfo' unitStruct.name '_sessioni' num2str(sessioni) ...
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
                end
            end
        end
    end
end






% stimes = spikeTimes.';
%
%                 %note the firing rate for every neuron
%                 tmax = max(stimes) - stimes(1);
%                 rate(uniti) = size(stimes,2)/(tmax);
%                 recTime(uniti) = tmax/(60);
%                 spikes = size(stimes,2);
%                 numofSpikes(uniti)= spikes;
%
%                 isi = diff(stimes*1000);
%                 nonISI = find(isi(1,:)<0);
%                 isi(nonISI) = [];
%                 log_isi = log10(isi);
%                 sdISI= std(isi);
%                 aveISI = mean(isi);
%                 cv(uniti) = sdISI/aveISI;
%                 meanISI(uniti) = aveISI;
%                 stdISI(uniti) = std(isi);
%
%                 isi = isi;
%                 ISIaxis1 = stimes(1:end-1);
%                 ISIaxis2 = stimes(2:end);
%                 ISIaxis = (ISIaxis1+ISIaxis2)/2;
%
%                 % calc stats
%                 mean_logisi = mean(log_isi);
%                 std_logisi = std(log_isi);
%                 skew_logisi = skewness(log_isi);
%                 kurt_logisi = kurtosis(log_isi);
%
%
%                 fprintf(1,'nspikes = %d\n', size(stimes,2));
%                 fprintf(1,'rate = %.1f Hz\n', size(stimes,2)/(tmax));
%
%
%
%                 %%%%%%%%%%%%%%%%%%%%%%%%
%                 % SHORT-LONG CLUSTERING
%                 %%%%%%%%%%%%%%%%%%%%%%%%
%
%                 % cutoff is the logisi value that separates short from long
%                 % list of cutoff values to test
%                 cutoffs = 0.05:0.05:10;
%                 sumdsq = zeros(size(cutoffs)); % summed distance-squared
%
%                 log_isi = log_isi.';
%                 X = [log_isi(1:end-1) log_isi(2:end)];
%
%                 for j = 1:length(cutoffs)
%                     cut = cutoffs(j);
%
%                     XSS = X(X(:,1) <= cut & X(:,2) <= cut, :);
%                     deltaSS = XSS - repmat(mean(XSS),size(XSS,1),1);
%                     dsqSS = sum(deltaSS.*deltaSS, 2);
%
%                     XSL = X(X(:,1) <= cut & X(:,2) > cut, :);
%                     deltaSL = XSL - repmat(mean(XSL),size(XSL,1),1);
%                     dsqSL = sum(deltaSL.*deltaSL, 2);
%
%                     XLS = X(X(:,1) > cut & X(:,2) <= cut, :);
%                     deltaLS = XLS - repmat(mean(XLS),size(XLS,1),1);
%                     dsqLS = sum(deltaLS.*deltaLS, 2);
%
%                     XLL = X(X(:,1) > cut & X(:,2) > cut, :);
%                     deltaLL = XLL - repmat(mean(XLL),size(XLL,1),1);
%                     dsqLL = sum(deltaLL.*deltaLL, 2);
%
%                     sumdsq(j) = sum(dsqSS) + sum(dsqSL) + sum(dsqLS) + sum(dsqLL);
%                 end
%
%                 %find the optimal cutoff value
%                 cut = cutoffs(sumdsq == min(sumdsq));
%                 cut = cut(1); %pick the last cut if two or more cutoffs are equally best
%                 XSS = X(X(:,1) <= cut & X(:,2) <= cut, :);
%                 XSL = X(X(:,1) <= cut & X(:,2) > cut, :);
%                 XLS = X(X(:,1) > cut & X(:,2) <= cut, :);
%                 XLL = X(X(:,1) > cut & X(:,2) > cut, :);
%
%                 % note the cut for every neuron
%                 noteCut(uniti) = (10^cut);
%
%
%                 %centroids
%                 cSS = mean(XSS);
%                 cSL = mean(XSL);
%                 cLS = mean(XLS);
%                 cLL = mean(XLL);
%
%                 % differece between the inter & intra burst spike interval
%                 difference(uniti) = 0.001*(10^(cLS(1)) - 10^(cSS(1)));
%
%
%                 %point number ratio between different clusters
%                 if size(XSS,1)>=1
%                     SStoLS = size(XSS,1)/size(XLS,1);
%                     SStoLL = size(XSS,1)/size(XLL,1);
%                     SStoSL = size(XSS,1)/size(XSL,1);
%                 else
%                     SStoLS = 0;
%                     SStoLL = 0;
%                     SStoSL = 0;
%                 end
%
%                     ratioSStoLS =[];
%                     ratioSStoLL =[];
%                     ratioSStoSL=[];
%                 ratioSStoLS =[ratioSStoLS SStoLS];
%                 ratioSStoLL =[ratioSStoLL SStoLL];
%                 ratioSStoSL =[ratioSStoSL SStoSL];
%
%
%
%                 %%%%%%%%%%%
%                 %GRAPHICS
%                 %%%%%%%%%%
%                 figure
%                 %spike train
%                 subplot(3,2,[1:2])
%                 plot([stimes', stimes'] ./ 60000, [0, 1], 'k-');
%                 ylim([0, 1]);
%                 xlabel ('time (minutes)');
% %                 title(source_files(i).name);
%
%                 subplot(3,2,3)
%                 bins = 1.0:0.1:6.0;
%                 hist(log_isi, bins);
%                 xlim([1,5]);
%                 str = sprintf('log isi (mean, std = %.2f, %.2f)', ...
%                     mean_logisi, ...
%                     std_logisi);
%                 xlabel(str)
%                 title('log isi')
%                 grid on
%
%                 subplot(3,2,4);
%                 plot(log_isi(1:end-1), log_isi(2:end), '.');
%                 xlabel('log isi (n)');
%                 ylabel('los isi (n+1)');
%                 axis square
% %                 axis([1 5 1 5]);
%                 grid on;
%
%                 %plot summed dsq
%                 subplot(3,2,5);
%                 plot(cutoffs, sumdsq, '*-');
%                 xlabel('logisi cutoff');
%                 ylabel('summed dist-squared');
%                 grid on;
%
%                 %plot clusters
%                 subplot(3,2,6)
%                 plot(XSS(:,1),XSS(:,2),'m.')
%                 hold on
%                 plot(XSL(:,1),XSL(:,2),'r.')
%                 plot(XLS(:,1),XLS(:,2),'g.')
%                 plot(XLL(:,1),XLL(:,2),'c.')
%
%                 % plot centroids
%                 plot(cSS(1), cSS(2),'ko');
%                 plot(cSL(1), cSL(2),'ko');
%                 plot(cLS(1), cLS(2),'ko');
%                 plot(cLL(1), cLL(2),'ko');
%
%                 title ('Cluster Assignments and Centroids');
%                 hold off
%                 axis square
% %                 axis([1 5 1 5]);
%                 grid on
%
%                 %generate bursts and intervals for the spike data based on cut
%                 cd('/Users/labmanager/Documents/MATLAB/burstAnalysiseNeuroPaper')
%                 [bursts] = findbursts(stimes, cut);
%                 allBursts = [];
%                 allBursts = [allBursts; bursts];
%
%                 %give burst infomation for the spike data
%                 BurstDuration(uniti) = mean(bursts(:,1));
%                 BurstSpike(uniti) = mean(bursts(:,2));
%                 BurstRate(uniti) = mean(bursts(:,3));
%                 FRinBursts(uniti) = mean(bursts(:,4));
%                 threshold(uniti) = mean(bursts(:,5));
%
%
%
%                 % try GMM with E-M
%                 X = [log_isi(1:end-1), log_isi(2:end)];
%                 k = 1:9;
%                 nK = numel(k);
%                 Sigma = {'diagonal','full'};
%                 nSigma = numel(Sigma);
%                 SharedCovariance = {false};
%                 SCtext = {'false'};
%                 nSC = numel(SharedCovariance);
%                 RegularizationValue = 0.01;
%                 options = statset('MaxIter',10000);
%
%                 % Preallocation
%                 gm = cell(nK,nSigma,nSC);
%                 aic = zeros(nK,nSigma,nSC);
%                 bic = zeros(nK,nSigma,nSC);
%                 converged = false(nK,nSigma,nSC);
%
%                 % Fit all models
%                 for m = 1:nSC;
%                     for j = 1:nSigma;
%                         for i = 1:nK;
%                             gm{i,j,m} = fitgmdist(X,k(i),...
%                                 'CovarianceType',Sigma{j},...
%                                 'SharedCovariance',SharedCovariance{m},...
%                                 'RegularizationValue',RegularizationValue,...
%                                 'Options',options);
%                             aic(i,j,m) = gm{i,j,m}.AIC;
%                             bic(i,j,m) = gm{i,j,m}.BIC;
%                             converged(i,j,m) = gm{i,j,m}.Converged;
%                         end
%                     end
%                 end
%
%                 allConverge = (sum(converged(:)) == nK*nSigma*nSC)
%
%                 figure;
%                 subplot(121)
%                 plot(aic(:,1))
%                 hold on; plot(aic(:,2))
% %                 bar(reshape(aic,nK,nSigma*nSC));
%                 title('AIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
%                 xlabel('$k$','Interpreter','Latex');
%                 ylabel('AIC');
%                 legend({'Diagonal-shared','Full-shared'});
%
%                 subplot(122)
%                 plot(bic(:,1))
%                 hold on; plot(bic(:,2))
%                 xlim([0 10])
% %                 bar(reshape(bic,nK,nSigma*nSC));
%                 title('BIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
%                 xlabel('$c$','Interpreter','Latex');
%                 ylabel('BIC');
%                 legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
%                     'Full-unshared'});
%
%                 [~,kclust] = min(bic(:,2));
%                 gmBest = gm{kclust,2};
%                 clusterX = cluster(gmBest,X);
%                 kGMM = gmBest.NumComponents;
%                 d = 500;
%                 x1 = linspace(min(X(:,1)) - 2,max(X(:,1)) + 2,d);
%                 x2 = linspace(min(X(:,2)) - 2,max(X(:,2)) + 2,d);
%                 [x1grid,x2grid] = meshgrid(x1,x2);
%                 X0 = [x1grid(:) x2grid(:)];
%                 mahalDist = mahal(gmBest,X0);
%                 threshold = sqrt(chi2inv(0.99,2));
%
%                 figure
%     %             subplot(232)
%                 h1 = gscatter(X(:,1),X(:,2),clusterX);
%                 hold on;
%                 for j = 1:kGMM;
%                     idx = mahalDist(:,j)<=threshold;
%                     Color = h1(j).Color*0.75 + -0.5*(h1(j).Color - 1);
% %                     h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
% %                     uistack(h2,'bottom');
%                 end
%                 h3 = plot(gmBest.mu(:,1),gmBest.mu(:,2),'kx','LineWidth',2,'MarkerSize',10);
% %                 set(gca,'xscale','log')
% %                 set(gca,'yscale','log')
%                 xlabel('Pre ISI')
%                 ylabel('Post ISI')
%                 title('Return plot')
%                 legend(h1,'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Location','NorthWest');
%                 hold off
%
%
% %
% %
% %
% %
% %
% %
% %
%                 FigH = figure;
%                 % return plot
%                 subplot(221)
%                 plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','MarkerSize', 0.05)
%                 preTimesOverThreshLong = find(allDiffSpikeTimesPre > log10(100));
%                 postTimesOverThreshLong = find(allDiffSpikeTimesPost > log10(100));
%                 %             loglog(allDiffSpikeTimes(:,1),allDiffSpikeTimes(:,2),'.')
%                 hold on
%                 x1 = allDiffSpikeTimesPre(preTimesOverThreshLong);
%                 y1 = allDiffSpikeTimesPost(preTimesOverThreshLong);
%                 %             plot(x1,y1,'r.')
%                 x2 = allDiffSpikeTimesPre(postTimesOverThreshLong);
%                 y2 = allDiffSpikeTimesPost(postTimesOverThreshLong);
%                 %             plot(x2,y2,'k.')
%                 bothTimesOverThresh = find(allDiffSpikeTimesPre(postTimesOverThreshLong) > log10(100));
%                 x3 = x2(bothTimesOverThresh);
%                 y3 = y2(bothTimesOverThresh);
%                 plot(x3,y3,'k.','MarkerSize', 0.05)
%                 preTimesOverThreshShort = find(allDiffSpikeTimesPre < log10(4));
%                 postTimesOverThreshShort = find(allDiffSpikeTimesPost < log10(4));
%                 x4 = allDiffSpikeTimesPre(preTimesOverThreshShort);
%                 y4 = allDiffSpikeTimesPost(preTimesOverThreshShort);
%                 %             plot(x4,y4,'b.')
%                 x5 = allDiffSpikeTimesPre(postTimesOverThreshShort);
%                 y5 = allDiffSpikeTimesPost(postTimesOverThreshShort);
%                 %             plot(x5,y5,'m.')
%                 bothTimesOverThreshShort = find(allDiffSpikeTimesPre(postTimesOverThreshShort) < log10(4));
%                 x6 = x5(bothTimesOverThreshShort);
%                 y6 = y5(bothTimesOverThreshShort);
%                 plot(x6,y6,'m.','MarkerSize', 0.05)
%                 xlabel('Pre ISI')
%                 ylabel('Post ISI')
%                 title('Return plot: All data')
%                 hold on
%                 plot([2 2],[log10(1) log10(4)],'--k')
%                 plot([2 4],[log10(4) log10(4)],'--k')
%                 %             plot([3 3],[log10(1) log10(4)],'--k')
%                 plot([log10(4) log10(4)],[0 10],'--k')
%                 xlim([0 4])
%                 ylim([0 4])
%                 xticks(0:4)
%                 xticklabels({'1', '10', '100', '1000'})
%                 xlabel('PreISI in ms on log scale')
%                 yticks(0:4)
%                 yticklabels({'1', '10', '100', '1000'})
%                 ylabel('PostISI in ms on log scale')
%
%                 % calculate % of bursting activity
%                 %             figure
%                 %             plot(allDiffSpikeTimesPre,allDiffSpikeTimesPost,'.','MarkerSize', 1)
%                 %             hold on
%                 %             plot(allDiffSpikeTimesPre(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4)),allDiffSpikeTimesPost(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4)),'m*','MarkerSize', 1)
%                 firstSpikeBurst = sum(allDiffSpikeTimesPre > log10(100) & allDiffSpikeTimesPost < log10(4),1);
%                 %             hold on
%                 %             plot(allDiffSpikeTimesPre(allDiffSpikeTimesPre < log10(4)),allDiffSpikeTimesPost(allDiffSpikeTimesPre < log10(4)),'k*','MarkerSize', 1)
%                 overFirstSpikeBurst = sum(allDiffSpikeTimesPre < log10(4),1);
%
%                 allISIs = size(allDiffSpikeTimesPre,1);
%                 burstISIs = overFirstSpikeBurst + firstSpikeBurst;
%                 percBurst(uniti) = (burstISIs*100)/allISIs;
%                 fprintf('\tBurst activity = %0.2f %% of %d spikes.\n', ...
%                     percBurst(uniti), allISIs);
%
%                 set(FigH, 'NumberTitle', 'off', ...
%                     'Name', sprintf('Burst activity = %0.2f %% of %d spikes.', percBurst(uniti), allISIs));
%
%                 timeXaxis = linspace(min(allDiffSpikeTimesPre), max(allDiffSpikeTimesPre),100);
%                 % plot logISI histogram
%                 subplot(222)
%                 %histogram(log10(allDiffSpikeTimesPre(2:end)),timeXaxis)
%                 histogram(allDiffSpikeTimesPre,timeXaxis)
%                 xticks(-1:3)
%                 xticklabels({'0.1','1', '10', '100', '1000'})
%                 xlabel('ISI in ms on log scale')
%                 title('Histogram logISI')
%
%                 % smoothing procedure to detect peaks
%                 %outputHistogram = histcounts(log10(allDiffSpikeTimesPre(2:end)),timeXaxis);
%                 outputHistogram = histcounts(allDiffSpikeTimesPre,timeXaxis);
%                 smoothISIhist = smooth(outputHistogram,'lowess');
%                 [pks,locs,widths] = findpeaks(smoothISIhist,'MinPeakHeight',2);
%                 hold on; plot(timeXaxis(1:end-1),smoothISIhist','m'); plot(timeXaxis(locs), pks', 'm*')
%
%                 subplot(234)
%                 plot(D.allUnitStructs{uniti}.meanWf,'k','LineWidth',3)
%                 hold on
%                 plot(D.allUnitStructs{uniti}.meanWf - (std(D.allUnitStructs{uniti}.wf,1)*3),'--k');
%                 plot(D.allUnitStructs{uniti}.meanWf + (std(D.allUnitStructs{uniti}.wf,1)*3),'--k');
%                 title('Spike waveform')