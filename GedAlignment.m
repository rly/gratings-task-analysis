% CAR LFP data for VEP data

%'/Volumes/scratch/rly/gratings-task-analysis/processed_data//M20170324/LFP_VEPM/M20170324-ind11-PUL-ch1-ch32-vepm1-vepm2-vepm3-vepm4-vepm5-CAR-responses-v12.mat'
cd('/Volumes/scratch/rly/gratings-task-analysis/processed_data//M20170324/LFP_VEPM/')
load('M20170324-ind12-PUL-ch33-ch64-vepm1-vepm2-vepm3-vepm4-vepm5-CAR-responses-v12.mat')

%%
clear all; close all
processedDataRootDir = '/Volumes/scratch/rly/gratings-task-analysis/processed_data/';
dataDirRoot = '/Volumes/kastner/ryanly/McCartney/merged';
muaDataDirRoot = '/Volumes/scratch/rly/simple-mua-detection/processed_data/';
recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
sessionInd = 15; %12 18
channelsToLoad = 33:64;
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, 'VEPM', 'LFP_VEPM', 0, 0, 1, 0);
sessionName = R.sessionName;
areaName = R.areaName;

if strcmp(sessionName, 'M20170127') || strcmp(sessionName, 'M20170130') || strcmp(sessionName, 'M20170201')
    origFlashEvents = D.events{6};
    preFlashesEvents = D.events{5};
else
    origFlashEvents = D.events{3};
    preFlashesEvents = D.events{2};
end

%[fixationAndLeverTimes,isMissingData] = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc)

%%
% flash related activity 200ms windows
postFlashWindowOffset = [0.025 0.225]; % seconds after flash
baselineWindowOffset = [-.175 .025]; % seconds around preflashevent
saccadeWindowOffset = [-.425 -.225]; % fixation onset at 325ms before preflashevent

% preprocess LFPs
Fs = D.lfpFs;
D.adjLfpsClean = D.adjLfps; %interpolateLfpOverSpikeTimes(D.adjLfps, channelsToLoad, Fs, D.allMUAStructs);
hiCutoffFreq = 100;
[channelDataCARNorm,channelDataNorm,commonAverageNorm,isNoisyChannel] = preprocessLfps(...
        D.adjLfpsClean, Fs, D.lfpNames, processedDataDir, [], hiCutoffFreq, 1, []);
D.adjLfpsClean = [];
flashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, origFlashEvents, ...
        postFlashWindowOffset, processedDataDir, [], 1, []);
preFlashEventsCleanWindowOffset = [min([baselineWindowOffset saccadeWindowOffset]) max([baselineWindowOffset saccadeWindowOffset])];
preFlashEventsClean = detectOutlierLfpEvents(channelDataCARNorm, Fs, D.lfpNames, preFlashesEvents, ...
        preFlashEventsCleanWindowOffset, processedDataDir, [], 1, []);
   
nFlashes = numel(flashEventsClean);
nTrials = numel(preFlashEventsClean);
startIndicesStim = round((flashEventsClean + postFlashWindowOffset(1)) * Fs); % time to index conversion
endIndicesStim = startIndicesStim + round(diff(postFlashWindowOffset) * Fs) - 1;
startIndicesBl = round((preFlashEventsClean + baselineWindowOffset(1)) * Fs); % time to index conversion
endIndicesBl = startIndicesBl + round(diff(baselineWindowOffset) * Fs) - 1;
startIndicesSac = round((preFlashEventsClean + saccadeWindowOffset(1)) * Fs); % time to index conversion
endIndicesSac = startIndicesSac + round(diff(saccadeWindowOffset) * Fs) - 1;
t = postFlashWindowOffset(1):1/Fs:postFlashWindowOffset(2)-1/Fs;
nTime = numel(t);
% assert(all(nTime == (endIndicesStim - startIndicesStim + 1)));

nChannels = length(channelsToLoad);
responses = nan(nChannels, nTime, nFlashes);
baseline = nan(nChannels, nTime, nTrials);
saccade = nan(nChannels, nTime, nTrials);
for j = 1:nChannels
    % can vectorize???
    for i = 1:nFlashes
        responses(j,:,i) = channelDataCARNorm(j,startIndicesStim(i):endIndicesStim(i));
    end
    for k = 1:nTrials
        baseline(j,:,k) = channelDataCARNorm(j,startIndicesBl(k):endIndicesBl(k));
        saccade(j,:,k) = channelDataCARNorm(j,startIndicesSac(k):endIndicesSac(k));
    end
end
allResponses = cat(3,responses,saccade);

% find first flahses, interflash time is typically 400ms (of which 100ms
% stim presentation)
firstFlashesTmp = [1; find(diff(startIndicesStim)>500)]; 
firstFlashes = firstFlashesTmp(boolean([diff(firstFlashesTmp)==4; 1]));
secondFlashesTmp = firstFlashesTmp + 1; 
secondFlashes = secondFlashesTmp(boolean([diff(secondFlashesTmp)==4; 1]));
thirdFlashesTmp = firstFlashesTmp + 2; 
thirdFlashes = thirdFlashesTmp(boolean([diff(thirdFlashesTmp)==4; 1]));
fourthFlashesTmp = firstFlashesTmp + 3; 
fourthFlashes = fourthFlashesTmp(boolean([diff(fourthFlashesTmp)==4; 1]));
firstFlashResponse = mean(responses(:,:,firstFlashes),3);
secondFlashResponse = mean(responses(:,:,secondFlashes),3);
thirdFlashResponse = mean(responses(:,:,thirdFlashes),3);
fourthFlashResponse = mean(responses(:,:,fourthFlashes),3);

% plot the LFPs to the different flashes with baseline correction
figure;
subplot(151)
imagesc(mean(bsxfun(@minus,responses(:,:,firstFlashes),squeeze(mean(mean(baseline,3),2))),3))
caxis1 = caxis;
subplot(152)
imagesc(mean(bsxfun(@minus,responses(:,:,secondFlashes),squeeze(mean(mean(baseline,3),2))),3))
caxis2 = caxis;
subplot(153)
imagesc(mean(bsxfun(@minus,responses(:,:,thirdFlashes),squeeze(mean(mean(baseline,3),2))),3))
caxis3 = caxis;
subplot(154)
imagesc(mean(bsxfun(@minus,responses(:,:,fourthFlashes),squeeze(mean(mean(baseline,3),2))),3))
caxis4 = caxis;
subplot(155)
imagesc(mean(baseline,3))
caxis5 = caxis;
caxisMin = min([caxis1 caxis2 caxis3 caxis4 caxis5]); caxisMax = max([caxis1 caxis2 caxis3 caxis4 caxis5]);
colorbar
colormap(getCoolWarmMap());
subplot(151)
caxis([caxisMin caxisMax])
title('First Flash')
subplot(152)
caxis([caxisMin caxisMax])
title('Second Flash')
subplot(153)
caxis([caxisMin caxisMax])
title('Third Flash')
subplot(154)
caxis([caxisMin caxisMax])
title('Fourth Flash')
subplot(155)
caxis([caxisMin caxisMax])
title('Baseline')

% plot the LFPs to the different flashes without baseline correction
figure;
subplot(151)
imagesc(firstFlashResponse)
caxis1 = caxis;
subplot(152)
imagesc(secondFlashResponse)
caxis2 = caxis;
subplot(153)
imagesc(thirdFlashResponse)
caxis3 = caxis;
subplot(154)
imagesc(fourthFlashResponse)
caxis4 = caxis;
subplot(155)
imagesc(mean(responses,3))
caxis5 = caxis;
caxisMin = min([caxis1 caxis2 caxis3 caxis4 caxis5]); caxisMax = max([caxis1 caxis2 caxis3 caxis4 caxis5]);
colorbar
colormap(getCoolWarmMap());
subplot(151)
caxis([caxisMin caxisMax])
title('First Flash')
subplot(152)
caxis([caxisMin caxisMax])
title('Second Flash')
subplot(153)
caxis([caxisMin caxisMax])
title('Third Flash')
subplot(154)
caxis([caxisMin caxisMax])
title('Fourth Flash')
subplot(155)
caxis([caxisMin caxisMax])
title('Averaged across Flashes')

% plot flash response and peri-saccade LFP
figure;
subplot(121)
imagesc(mean(bsxfun(@minus,responses, squeeze(mean(mean(baseline,2),3))),3))
caxis1 = caxis;
subplot(122)
imagesc(mean(bsxfun(@minus,saccade, squeeze(mean(mean(baseline,2),3))),3))
caxis2 = caxis;
caxisMin = min([caxis1 caxis2]); caxisMax = max([caxis1 caxis2]);
colorbar
colormap(getCoolWarmMap());
subplot(121)
caxis([caxisMin caxisMax])
title('Flash response')
subplot(122)
caxis([caxisMin caxisMax])
title('Saccade response')

figure;
subplot(121)
imagesc(mean(responses,3))
caxis1 = caxis;
subplot(122)
imagesc(mean(saccade,3))
caxis2 = caxis;
caxisMin = min([caxis1 caxis2]); caxisMax = max([caxis1 caxis2]);
colorbar
colormap(getCoolWarmMap());
subplot(121)
caxis([caxisMin caxisMax])
title('Flash response')
subplot(122)
caxis([caxisMin caxisMax])
title('Saccade response')

% create covariance matrix separately for saccade and flash
% create cov matrix as in Mike's code (method #1)
flashesCont = reshape(responses(:,:,firstFlashes),nChannels,[]);
flashesCont = bsxfun(@minus,flashesCont,mean(flashesCont,2));
flashesCov = (flashesCont*flashesCont') / (size(flashesCont,2)-1);
% create cov matrix as in paper (method #2)
for n = 1:size(responses(:,:,firstFlashes),3)
    flashesCov2all(n,:,:) = responses(:,:,firstFlashes(n))*responses(:,:,firstFlashes(n))' / (size(responses(:,:,firstFlashes(n)),2)-1);
end
flashesCov2 = squeeze(mean(flashesCov2all,1));
figure; subplot(121); imagesc(flashesCov); title('On concatenated trials')
subplot(122); imagesc(flashesCov2); title('Per trial and then averaged')
% concl. methods are not identical but yield very similar results and
% pattern is the same

% apply GED separately for saccade and flash evoked activity
% first calculate GED for saccade activity
saccadeCovAll = nan(size(saccade,3),nChannels,nChannels);
for n = 1:size(saccade,3)
    saccadeCovAll(n,:,:) = saccade(:,:,(n))*saccade(:,:,(n))' / (size(saccade(:,:,(n)),2)-1);
end
saccadeCov = squeeze(mean(saccadeCovAll,1));
% next calculate GED for baseline activity
baselineCovAll = nan(size(baseline,3),nChannels,nChannels);
for n = 1:size(baseline,3)
    baselineCovAll(n,:,:) = baseline(:,:,(n))*baseline(:,:,(n))' / (size(baseline(:,:,(n)),2)-1);
end
baselineCov = squeeze(mean(baselineCovAll,1));
% get GED
[evecs, evals] = eig(flashesCov2,baselineCov);
[~,eigidx] = sort(diag(evals));
evecs = evecs(:,eigidx);
figure; 
subplot(141); imagesc(flashesCov2); title('covariance matrix S')
%caxis1 = caxis;
caxis([-0.5 0.5])
subplot(142); imagesc(evecs); title('eigen vectors');
caxis([-1 1])
subplot(143); imagesc(evals(eigidx,eigidx)); title(['eigen values' sprintf('\\lambda')])
subplot(144); imagesc(baselineCov); title('covariance matrix R')
%caxis2 = caxis;
caxis([-0.5 0.5])
%caxisMin = min([caxis1 caxis2]); caxisMax = max([caxis1 caxis2]);
colorbar
colormap(getCoolWarmMap());
% subplot(141)
% caxis([caxisMin caxisMax])
% subplot(144)
% caxis([caxisMin caxisMax])

firstflashGed = reshape( (reshape(responses(:,:,firstFlashes),nChannels,size(responses(:,:,firstFlashes),2)*size(responses(:,:,firstFlashes),3))' * evecs)',...
    nChannels,nTime,length(firstFlashes));
firstflashGedCompMaps = reshape( (flashesCov2' * evecs(:,end))',1,nChannels,1);
figure; 
plot(firstflashGedCompMaps,1:32); set(gca, 'YDir', 'reverse');
title(['Eigen value: ' num2str(evals(end,end))])
findchangepts(firstflashGedCompMaps,'Statistic','mean','MinThreshold',.05)
figure; plot(mean(firstflashGed(end,:,:),3)) % CHECK WITH RYAN!!!!!!!!!!!!

[evecsSacc, evalsSacc] = eig(saccadeCov,baselineCov);
[~,eigidxSacc] = sort(diag(evalsSacc));
evecsSacc = evecsSacc(:,eigidxSacc);
figure; subplot(121); imagesc(evecsSacc); subplot(122); imagesc(evalsSacc(eigidxSacc,eigidxSacc))

figure; 
subplot(141); imagesc(saccadeCov); title('covariance matrix S')
%caxis1 = caxis;
caxis([-0.5 0.5])
subplot(142); imagesc(evecsSacc); title('eigen vectors');
caxis([-1 1])
subplot(143); imagesc(evalsSacc(eigidxSacc,eigidxSacc)); title(['eigen values' sprintf('\\lambda')])
subplot(144); imagesc(baselineCov); title('covariance matrix R')
%caxis2 = caxis;
caxis([-0.5 0.5])
%caxisMin = min([caxis1 caxis2]); caxisMax = max([caxis1 caxis2]);
colorbar
colormap(getCoolWarmMap());

saccadeGed = reshape( (reshape(saccade,nChannels,size(saccade,2)*size(saccade,3))' * evecsSacc)',...
    nChannels,nTime,size(saccade,3));
saccadeGedCompMaps = reshape( (saccadeCov' * evecsSacc(:,end))',1,nChannels,1);
figure; plot(saccadeGedCompMaps,1:32); set(gca, 'YDir', 'reverse');
title(['Eigen value: ' num2str(evalsSacc(end,end))])
findchangepts(saccadeGedCompMaps,'Statistic','mean','MinThreshold',.05)
figure; plot(mean(saccadeGed(end,:,:),3))

%% collapse across stimulus conditions
for n = 1:size(responses,3)
    allFlashesCovAll(n,:,:) = responses(:,:,n)*responses(:,:,n)' / (size(responses(:,:,n),2)-1);
end
stimulusCovAll = cat(1,allFlashesCovAll,saccadeCovAll);
stimulusCov = squeeze(mean(stimulusCovAll,1));
% get GED
[evecs, evals] = eig(stimulusCov,baselineCov);
[~,eigidx] = sort(diag(evals));
evecs = evecs(:,eigidx);
figure; 
subplot(141); imagesc(stimulusCov); title('covariance matrix S')
%caxis1 = caxis;
caxis([-0.5 0.5])
subplot(142); imagesc(evecs); title('eigen vectors');
caxis([-1 1])
subplot(143); imagesc(evals(eigidx,eigidx)); title(['eigen values' sprintf('\\lambda')])
subplot(144); imagesc(baselineCov); title('covariance matrix R')
%caxis2 = caxis;
caxis([-0.5 0.5])
%caxisMin = min([caxis1 caxis2]); caxisMax = max([caxis1 caxis2]);
colorbar
colormap(getCoolWarmMap());
% subplot(141)
% caxis([caxisMin caxisMax])
% subplot(144)
% caxis([caxisMin caxisMax])

stimulusGed = reshape( (reshape(responses,nChannels,size(responses,2)*size(responses,3))' * evecs)',...
    nChannels,nTime,size(responses,3));
stimulusGedCompMaps = reshape( (stimulusCov' * evecs(:,end))',1,nChannels,1);
figure; 
plot(stimulusGedCompMaps,1:32); set(gca, 'YDir', 'reverse');
title(['Eigen value: ' num2str(evals(eigidx(end),eigidx(end)))])
findchangepts(stimulusGedCompMaps,'Statistic','mean','MinThreshold',.05)
figure; plot(mean(stimulusGed(end,:,:),3)) % CHECK WITH RYAN!!!!!!!!!!!!

figure
bar(sort(diag(evals),'descend'))

% permutation test
nPerm = 500; permCov = cat(1,stimulusCovAll,baselineCovAll); nBaseline = size(baselineCovAll,1);
permEvals = nan(nPerm,nChannels);
for permi = 1:nPerm
    permCovShuffled = permCov(randperm(size(permCov,1)),:,:); 
    baselineCovShuffled = squeeze(mean(permCovShuffled(1:nBaseline,:,:),1));
    stimulusCovShuffled = squeeze(mean(permCovShuffled(nBaseline+1:end,:,:),1));
    [evecsShuffled, evalsShuffled] = eig(stimulusCovShuffled,baselineCovShuffled);
    permEvals(permi,:) = sort(diag(evalsShuffled),'descend')';
end
figure;
plot(sort(real(permEvals(:,1))))
figure
histogram(real(permEvals(:,1)),'BinWidth',0.1,'Normalization','probability')
figure
bar(sorterPermEvals)

% compare GED applied on 1st and 4th flash
for n = 1:size(responses(:,:,fourthFlashes),3)
    fourthFlashesCovAll(n,:,:) = responses(:,:,fourthFlashes(n))*responses(:,:,fourthFlashes(n))' / (size(responses(:,:,fourthFlashes(n)),2)-1);
end
fourthFlashesCov = squeeze(mean(fourthFlashesCovAll,1));
[evecsFourth, evalsFourth] = eig(fourthFlashesCov,baselineCov);
[~,eigidxFourth] = sort(diag(evalsFourth));
evecsFourth = evecsFourth(:,eigidxFourth);
figure; subplot(121); imagesc(evecsFourth); subplot(122); imagesc(evalsFourth(eigidxFourth,eigidxFourth))
fourthflashGed = reshape( (reshape(responses(:,:,fourthFlashes),nChannels,size(responses(:,:,fourthFlashes),2)*size(responses(:,:,fourthFlashes),3))' * evecsFourth)',...
    nChannels,nTime,length(fourthFlashes));
fourthflashGedCompMaps = reshape( (fourthFlashesCov' * evecsFourth(:,end))',1,nChannels,1);
figure; plot(fourthflashGedCompMaps,1:32); set(gca, 'YDir', 'reverse');
title(['Eigen value: ' num2str(evalsFourth(end,end))])
findchangepts(fourthflashGedCompMaps,'Statistic','mean','MinThreshold',.05)
figure; plot(mean(fourthflashGed(end,:,:),3))

% create covariance matrices
responsesCont = reshape(allResponses, nChannels, []);
responsesCont = bsxfun(@minus, responsesCont, mean(responsesCont,2));
responsesCov = (responsesCont*responsesCont') / size(responsesCont,2);
baselineCont = reshape(baseline, nChannels, []);
baselineCont = bsxfun(@minus, baselineCont, mean(baselineCont,2));
baselineCov = (baselineCont*baselineCont') / size(baselineCont,2);
[evecs, evals] = eig(responsesCov,baselineCov);
[~,eigidx] = sort(diag(evals));
evecs = evecs(:,eigidx);

responsesGed = reshape( (responsesCont' * evecs)',nChannels,nTime,nFlashes+nTrials);
responsesGedCompMaps = reshape( (responsesCov' * evecs(:,end))',1,nChannels,1);
responsesGedCompMaps2 = reshape( (responsesCov' * evecs(:,2))',1,nChannels,1);

findchangepts(responsesGedCompMaps,'Statistic','mean','MinThreshold',.05)
findchangepts(responsesGedCompMaps2,'Statistic','mean','MinThreshold',.05)

figure;
plot(responsesGedCompMaps,1:32)
set(gca, 'YDir', 'reverse');
title(['1st component eval = ' num2str(evecs(1))])

figure;
plot(responsesGedCompMaps2,1:32)
set(gca, 'YDir', 'reverse');
title(['2nd component eval = ' num2str(evecs(2))])


baselineGed = reshape( (baselineCont' * evecs)',nChannels,nTime,nTrials);

% plot it
averageResponse = mean(allResponses, 3); % average across flashes
averageResponseGED = mean(responsesGed, 3); % average across flashes
responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(averageResponse)));
figure;
subplot(221)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
%     plot(t, responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
subplot(222)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:2%nChannels
%     plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
    plot(t, responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponseGED(j,:) - j*responsePlotYOffset) - 1]);
    end
end
legend({'comp1' 'comp2'})
baselineResponse = mean(baseline, 3); % average across flashes
baselineResponseGED = mean(baselineGed, 3); % average across flashes
subplot(223)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*baselineResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
subplot(224)
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:2%nChannels
    plot(t, responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*baselineResponseGED(j,:) - j*responsePlotYOffset) - 1]);
    end
end

figure;
subplot(131)
imagesc(baselineResponse)
caxis1 = caxis;
subplot(132)
imagesc(averageResponse)
caxis2 = caxis;
subplot(133)
imagesc(averageSaccade)
caxis3 = caxis;
caxisMin = min([caxis1 caxis2 caxis3]); caxisMax = max([caxis1 caxis2 caxis3]);
subplot(131)
caxis([caxisMin caxisMax])
title(sprintf('%s %s - Baseline (N=%d) (%s)', sessionName, areaName, nTrials, ref));
subplot(132)
caxis([caxisMin caxisMax])
title(sprintf('Full-Field Flash (N=%d)', nFlashes));
subplot(133)
caxis([caxisMin caxisMax])
title(sprintf('Saccade (N=%d)', nTrials));
colorbar
colormap(getCoolWarmMap());

figure;
subplot(131)
imagesc(baselineResponse-baselineResponse)
caxis1 = caxis;
subplot(132)
imagesc(averageResponse-baselineResponse)
caxis2 = caxis;
subplot(133)
imagesc(averageSaccade-baselineResponse)
caxis3 = caxis;
caxisMin = min([caxis1 caxis2 caxis3]); caxisMax = max([caxis1 caxis2 caxis3]);
subplot(131)
caxis([caxisMin caxisMax])
title(sprintf('%s %s - Baseline (N=%d) (%s)', sessionName, areaName, nTrials, ref));
subplot(132)
caxis([caxisMin caxisMax])
title(sprintf('Full-Field Flash (N=%d)', nFlashes));
subplot(133)
caxis([caxisMin caxisMax])
title(sprintf('Saccade (N=%d)', nTrials));
colorbar
colormap(getCoolWarmMap());


xBounds = [-0.1 0.3]; ref = 'CAR';
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageResponse - baselineResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);



averageSaccade = mean(saccade,3);

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageSaccade - baselineResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Saccade (N=%d) (%s)', sessionName, areaName, nTrials, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);






[tBaseline,averageBaselineResponses] = computeResponsesInWindow(channelDataCARNorm, preFlashesEvents, ...
        baselineWindowOffset, Fs);
[tBaselineExp,averageBaselineResponsesExpCAR] = computeResponsesInWindow(channelDataCARNorm, preFlashesEvents, ...
        expandedPlotWindowOffset, Fs);

    
%% save responses to mat file
saveFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDir, plotFileNamePrefix, ref, v);
save(saveFileName, 'D', 'R', 'responses', 't', 'periFlashWindowOffset', 'isNoisyChannel');

%% subtract out baseline
% units are standard deviations from baselined mean
averageResponse = mean(responses, 3); % average across flashes

% subtract mean baseline activity (-0.1, 0] seconds before flash
flashBaselineWindowOffset = [-0.1 0];

% TODO make so that baseline can have arbitrary end time
indexFlashTime = -round(postFlashWindowOffset(1) * Fs);
indexStartBaselineFlashTime = -round((postFlashWindowOffset(1) - flashBaselineWindowOffset(1)) * Fs) + 1;
for j = 1:nChannels
    averageResponse(j,:) = averageResponse(j,:) - mean(averageResponse(j,indexStartBaselineFlashTime:indexFlashTime));
end

%% staggered channel line plot of average visually evoked LFP
xBounds = [-0.1 0.3];

responsePlotYOffset = 1;
responsePlotYScale = 5/max(max(abs(averageResponse)));

figure_tr_inch(8, 10);
hold on;
maxYLim = 3;
minYLim = -(nChannels + 2);
for j = 1:nChannels
    plot(t, responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset);
    if j == 1
        maxYLim = max([maxYLim max(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) + 1]);
    elseif j == nChannels
        minYLim = min([minYLim min(responsePlotYScale*averageResponse(j,:) - j*responsePlotYOffset) - 1]);
    end
end
xlim(postFlashWindowOffset);
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);
set(gca, 'FontSize', 16);
ax = ancestor(gca, 'axes');
ax.YAxis.FontSize = 8; % change y tick font size without changing x tick
ylabel('Channel Number (1 = topmost)', 'FontSize', 16);
xlabel('Time from Flash Onset (s)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));

origYLim = [minYLim maxYLim]; % ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim(origYLim);

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, minYLim+1, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpLines-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);



averageSaccade = mean(saccade,3);

figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageSaccade);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Saccade (N=%d) (%s)', sessionName, areaName, nTrials, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, nChannels, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpColor-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

%% plot boxes around areas abs() > thresh over last plot and re-save
boxAbsThresh = 0.25;
strongResponseGroups = bwlabel(abs(averageResponse) > boxAbsThresh, 4);
strongResponseGroupBoundaries = bwboundaries(strongResponseGroups);
for k = 1:length(strongResponseGroupBoundaries)
   boundary = strongResponseGroupBoundaries{k};
   plot(t(boundary(:,2)), boundary(:,1), 'Color', 0.3*ones(3, 1), 'LineWidth', 2)

   minBY = min(boundary(:,1));
   maxBY = max(boundary(:,1));
   text(-0.01, minBY, sprintf('%d', minBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
   text(-0.01, maxBY, sprintf('%d', maxBY), 'Color', 0.3*ones(3, 1), 'HorizontalAlignment', 'right');
end

plotFileName = sprintf('%s/%s-%s-lfpColorBounds-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');

%% plot vertical line plot of mean activity between 35 and 45 ms after flash
% TODO mark which ones are more than 2 SDs from baseline for that channel

earlyActivityWindowOffset = [0.035 0.045];
earlyActivityTLogical = getTimeLogicalWithTolerance(t, earlyActivityWindowOffset);
meanEarlyActivity = mean(averageResponse(:,earlyActivityTLogical), 2);
channelIndices = 1:nChannels;

figure_tr_inch(7, 8);
subaxis(1, 1, 1, 'ML', 0.12, 'MB', 0.11, 'MR', 0.05);
hold on;
plot(meanEarlyActivity, channelIndices, '.--', 'MarkerSize', 25, 'LineWidth', 2, 'Color', lines(1));
plot([0 0], channelIndices([1 end]) + [-1 1], '-', 'Color', 0.3*ones(3, 1));
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel(sprintf('Mean SDs Early Response to Flash (%d-%d ms)', earlyActivityWindowOffset * 1000));
ylabel('Channel Number');
grid on;
ylim(channelIndices([1 end]) + [-1 1]);
xlim([-1 1] * max(abs(xlim())));
title(sprintf('%s %s - Mean Early Response (%s)', sessionName, areaName, ref));

plotFileName = sprintf('%s/%s-%s-meanEarlyResponse-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');





%% staggered channel color plot of average visually evoked LFP
figure_tr_inch(8, 10);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, 1:nChannels, averageResponse);
set(gca, 'YDir', 'reverse');
plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
xlabel('Time from Flash Onset (s)');
ylabel('Channel Number (1 = topmost)');
title(sprintf('%s %s - Response to Full-Field Flash (N=%d) (%s)', sessionName, areaName, nFlashes, ref));
maxCAxis = max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
colormap(getCoolWarmMap());
colorbar;
set(gca, 'FontSize', 16);

if (nChannels > 31 && strcmp(ref, 'BIP')) || (nChannels > 32 && strcmp(ref, 'CAR'))
    plot(xlim(), [nChannels/2+0.5 nChannels/2+0.5], 'Color', 0.3*ones(3, 1));
end

% plot early latency line at 35 ms, 45 ms
plot([0.035 0.035], [-1000 1000], 'm-');
plot([0.045 0.045], [-1000 1000], 'm-');
text(0.050, nChannels, '35-45 ms', 'Color', 'm');

plotFileName = sprintf('%s/%s-%s-lfpColor-v%d.png', processedDataDir, plotFileNamePrefix, ref, v);
%export_fig(plotFileName, '-nocrop');





% eof