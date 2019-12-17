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
                        
                        % convolve with gaussian
                        kernelSigma = 0.1;
                        x = -5 * kernelSigma * 1000 : 5 * kernelSigma * 1000;
                        gaussian = normpdf(x, 0, kernelSigma * 1000);  % does the same thing as computing a gaussian with an equation as you had done
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;
                        allSpikesGauss = conv(binarySpikeTrain, gaussian);

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
                        binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;
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

