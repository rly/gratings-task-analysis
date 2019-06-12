function quickSpdfAllEventsRowPopMeanRunner(name, isIncluded, spdfInfo, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        yBounds, isShowLabels, plotFileName)

fprintf('\n');
fprintf('Plotting normalized mean SPDFs for group %s...\n', name);

cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(isIncluded,:));
cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(isIncluded,:));
arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(isIncluded,:));
arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(isIncluded,:));
targetDimSpdfInRFNormSub = (spdfInfo.targetDimSpdfInRFNorm(isIncluded,:));
targetDimSpdfExRFNormSub = (spdfInfo.targetDimSpdfExRFNorm(isIncluded,:));
exitFixationSpdfInRFNormSub = (spdfInfo.exitFixationSpdfInRFNorm(isIncluded,:));
exitFixationSpdfExRFNormSub = (spdfInfo.exitFixationSpdfExRFNorm(isIncluded,:));

fprintf('\t%s: %d units\n', name, sum(isIncluded));

fprintf('Saving to %s...\n', plotFileName);

quickSpdfAllEvents4InARowPopMean(...
        cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
        arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
        targetDimSpdfExRFNormSub, exitFixationSpdfInRFNormSub, exitFixationSpdfExRFNormSub, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, yBounds, isShowLabels, plotFileName);    

