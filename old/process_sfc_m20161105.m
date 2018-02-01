function process_sfc_m20161105(sessionName, spikeVar, spikeVarName, ...
        adjLfp, lfpVarName, nLoc, params, ...
        cueOnset, cueOnsetByLoc, ...
        arrayOnset, arrayOnsetByLoc, ...
        arrayOnsetRel, arrayOnsetHold, ...
        arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
        arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
        targetDim, targetDimByLoc, ...
        targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc)

%%
params.pad = 2;
Fs = params.Fs;

%% align spikes to CUE
periCueOnsetWindow = [1 1];
cueOnsetSpikeTimes = createnonemptydatamatpt(spikeVar, cueOnset, periCueOnsetWindow);

cueOnsetSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    cueOnsetSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{i}, periCueOnsetWindow);
end

%% align spikes to ARRAY
periArrayOnsetWindow = [1 1]; % seconds before, after

arrayOnsetSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnset, periCueOnsetWindow);
arrayOnsetRelSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetRel, periArrayOnsetWindow);
arrayOnsetHoldSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetHold, periArrayOnsetWindow);

arrayOnsetSpikeTimesByLoc = cell(nLoc,1);
arrayOnsetRelSpikeTimesByLoc = cell(nLoc,1);
arrayOnsetHoldSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    arrayOnsetSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{i}, periArrayOnsetWindow);
    arrayOnsetRelSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetRelByLoc{i}, periArrayOnsetWindow);
    arrayOnsetHoldSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{i}, periArrayOnsetWindow);
end

%% align spikes to TARGET DIM
periTargetDimWindow = [1 1];
targetDimSpikeTimes = createnonemptydatamatpt(spikeVar, targetDim, periTargetDimWindow);

targetDimSpikeTimesByLoc = cell(nLoc,1);
targetDimShortHoldDurSpikeTimesByLoc = cell(nLoc,1);
targetDimLongHoldDurSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    targetDimSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimByLoc{i}, periTargetDimWindow);
    targetDimShortHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimShortHoldDurByLoc{i}, periTargetDimWindow);
    targetDimLongHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimLongHoldDurByLoc{i}, periTargetDimWindow);
end

%% align lfp to events
lfpAroundCueOnsetByLoc = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundCueOnsetByLoc{i} = createdatamatc(adjLfp, cueOnsetByLoc{i}, Fs, periCueOnsetWindow);
end

lfpAroundArrayOnsetByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLoc = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetByLoc{i} = createdatamatc(adjLfp, arrayOnsetByLoc{i}, Fs, periArrayOnsetWindow);
    lfpAroundArrayOnsetHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, periArrayOnsetWindow);
    lfpAroundArrayOnsetShortHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetShortHoldByLoc{i}, Fs, periArrayOnsetWindow);
    assert(~any(any(isnan(lfpAroundArrayOnsetShortHoldByLoc{i}))));
    lfpAroundArrayOnsetLongHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetLongHoldByLoc{i}, Fs, periArrayOnsetWindow);
    assert(~any(any(isnan(lfpAroundArrayOnsetLongHoldByLoc{i}))));
end

lfpAroundTargetDimByLoc = cell(nLoc, 1);
lfpAroundTargetDimShortHoldByLoc = cell(nLoc, 1);
lfpAroundTargetDimLongHoldByLoc = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundTargetDimByLoc{i} = createdatamatc(adjLfp, targetDimByLoc{i}, Fs, periTargetDimWindow);
    lfpAroundTargetDimShortHoldByLoc{i} = createdatamatc(adjLfp, targetDimShortHoldDurByLoc{i}, Fs, periTargetDimWindow);
    lfpAroundTargetDimLongHoldByLoc{i} = createdatamatc(adjLfp, targetDimLongHoldDurByLoc{i}, Fs, periTargetDimWindow);
end

%% coherence parameters
cohMovingWindow = [0.3 0.03];  % 300ms width, slide every 30ms
fscorr = 0;

nTimeCueSFC = 57;
nTimeArrSFC = 57;
nFreqSFC = 205;

inRFLoc = 3;
exRFLoc = 1;

%% spike-field coherence
[cueInRFC,~,~,~,~,tCue,f] = cohgramcpt(...
        lfpAroundCueOnsetByLoc{inRFLoc}, cueOnsetSpikeTimesByLoc{inRFLoc}, ...
        cohMovingWindow, params, fscorr);
    
[cueExRFC,~,~,~,~,tCue2,f2] = cohgramcpt(...
        lfpAroundCueOnsetByLoc{exRFLoc}, cueOnsetSpikeTimesByLoc{exRFLoc}, ...
        cohMovingWindow, params, fscorr);

[arrInRFC,~,~,~,~,tArr,f3] = cohgramcpt(...
        lfpAroundArrayOnsetByLoc{inRFLoc}, arrayOnsetSpikeTimesByLoc{inRFLoc}, ...
        cohMovingWindow, params, fscorr);

[arrExRFC,~,~,~,~,tArr2,f4] = cohgramcpt(...
        lfpAroundArrayOnsetByLoc{exRFLoc}, arrayOnsetSpikeTimesByLoc{exRFLoc}, ...
        cohMovingWindow, params, fscorr);
    
[dimInRFC,~,~,~,~,tDim,f5] = cohgramcpt(...
        lfpAroundTargetDimByLoc{inRFLoc}, targetDimSpikeTimesByLoc{inRFLoc}, ...
        cohMovingWindow, params, fscorr);

[dimExRFC,~,~,~,~,tDim2,f6] = cohgramcpt(...
        lfpAroundTargetDimByLoc{exRFLoc}, targetDimSpikeTimesByLoc{exRFLoc}, ...
        cohMovingWindow, params, fscorr);   
    
numInRFTrials = numel(cueOnsetSpikeTimesByLoc{inRFLoc});
numExRFTrials = numel(cueOnsetSpikeTimesByLoc{exRFLoc});

conditionsInfo = var2struct(inRFLoc, exRFLoc, periCueOnsetWindow, periArrayOnsetWindow, periTargetDimWindow);
conditionsInfo.numTrials.inRF = numInRFTrials;
conditionsInfo.numTrials.exRF = numExRFTrials;

% sanity checks
assert(isequal(size(cueInRFC), size(cueExRFC), [nTimeCueSFC nFreqSFC]));
assert(isequal(size(arrInRFC), size(arrExRFC), [nTimeArrSFC nFreqSFC]));
assert(all(tCue == tCue2));
assert(all(tArr == tArr2));
assert(isequal(f, f2, f3, f4, f5, f6));
assert(numel(cueOnsetSpikeTimesByLoc{inRFLoc}) == numel(arrayOnsetSpikeTimesByLoc{inRFLoc}));
assert(numel(cueOnsetSpikeTimesByLoc{exRFLoc}) == numel(arrayOnsetSpikeTimesByLoc{exRFLoc}));

%% plot
plotFileName = sprintf('%s-%s-%s-sfc.png', sessionName, spikeVarName, lfpVarName);
plotSessionCueArrDimSpikeFieldCoherence(sessionName, ...,
        sprintf('%s-%s', spikeVarName, lfpVarName), ...
        cueInRFC, cueExRFC, arrInRFC, arrExRFC, dimInRFC, dimExRFC, ...
        tCue, tArr, tDim, f, conditionsInfo, ...
        0, plotFileName);
    

