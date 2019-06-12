function quickImagePlotAllEventsRowRunner(groupName, isIncluded, spdfInfo, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        diffPlotFileName, inRFPlotFileName)

fprintf('\n');
fprintf('Plotting normalized mean SPDFs for group %s...\n', groupName);

cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(isIncluded,:));
cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(isIncluded,:));
arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(isIncluded,:));
arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(isIncluded,:));
targetDimSpdfInRFNormSub = (spdfInfo.targetDimSpdfInRFNorm(isIncluded,:));
targetDimSpdfExRFNormSub = (spdfInfo.targetDimSpdfExRFNorm(isIncluded,:));
exitFixationSpdfInRFNormSub = (spdfInfo.exitFixationSpdfInRFNorm(isIncluded,:));
exitFixationSpdfExRFNormSub = (spdfInfo.exitFixationSpdfExRFNorm(isIncluded,:));

fprintf('\t%s: %d units\n', groupName, sum(isIncluded));

isDiff = 1;
quickImagePlotAllEvents4InARow(...
        cueOnsetSpdfInRFNormSub - cueOnsetSpdfExRFNormSub, ...
        arrayOnsetHoldSpdfInRFNormSub - arrayOnsetHoldSpdfExRFNormSub, ...
        targetDimSpdfInRFNormSub - targetDimSpdfExRFNormSub, ...
        exitFixationSpdfInRFNormSub - exitFixationSpdfExRFNormSub, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        isDiff, diffPlotFileName);

isDiff = 0;
quickImagePlotAllEvents4InARow(...
        cueOnsetSpdfInRFNormSub, ...
        arrayOnsetHoldSpdfInRFNormSub, ...
        targetDimSpdfInRFNormSub, ...
        exitFixationSpdfInRFNormSub, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        isDiff, inRFPlotFileName);
