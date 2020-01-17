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
idx = 1;
iidx = 1;
for sessioni = 1:numel(sessionInfo{1})
    %     if ismember(sessioni,burstsession)
    clearvars -except sessionInfo sessioni outputDir counterFR nUnitsPerSession v burstsession burstunit ...
        pBurstQual percentageBurst idx iidx
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
%         [~,iidx] = ismember(sessioni,burstsession);
        %             if ismember(uniti,burstunit(iidx,:))
        
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
%                 if isnan(unitStruct.unitStartTime)
                    startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>spikeTimes(1),1,'first'));
%                 else
%                     startTime = UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue(find(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue>unitStruct.unitStartTime,1,'first'));
%                 end
                
                % Match spikeTimes with last lever release
%                 if isnan(unitStruct.unitEndTime)
                    endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<spikeTimes(end),1,'last'));
%                 else
%                     endTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(find(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice<unitStruct.unitEndTime,1,'last'));
%                 end
                
                spikeTimes = spikeTimes(spikeTimes >= startTime & spikeTimes <= endTime);
                spikeTimes2use = (spikeTimes - spikeTimes(1))';
                
                % calculate percentage of spikes >= 200 Hz firing rate
                percentageBurst(idx) =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);

                if percentageBurst(idx) > 1
                    pBurstQual{iidx} = unitStruct.name;
                    iidx = iidx + 1;
                else
                end
                idx = idx + 1;

                clear binarySpikeTrain spikeTimes spikeTimes2use endTime startTime
                
            else
                counterFR(sessioni) = 1;
            end
        end
    end
end
%     end
% end


save('womelsdorfInfo','percentageBurst','pBurstQual')

% load('womelsdorfInfo','percentageBurst')
% sum(percentageBurst >= 1)






% eof