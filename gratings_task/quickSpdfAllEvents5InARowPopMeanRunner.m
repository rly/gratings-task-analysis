function quickSpdfAllEvents5InARowPopMeanRunner(name, isIncluded, spdfInfo, ...
        enterFixationT, cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        summaryDataDir, v)

fprintf('\n');
fprintf('Plotting normalized mean SPDFs for group %s...\n', name);

enterFixationSpdfInRFNormSub = (spdfInfo.enterFixationSpdfInRFNorm(isIncluded,:));
enterFixationSpdfExRFNormSub = (spdfInfo.enterFixationSpdfExRFNorm(isIncluded,:));
cueOnsetSpdfInRFNormSub = (spdfInfo.cueOnsetSpdfInRFNorm(isIncluded,:));
cueOnsetSpdfExRFNormSub = (spdfInfo.cueOnsetSpdfExRFNorm(isIncluded,:));
arrayOnsetHoldSpdfInRFNormSub = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(isIncluded,:));
arrayOnsetHoldSpdfExRFNormSub = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(isIncluded,:));
targetDimSpdfInRFNormSub = (spdfInfo.targetDimSpdfInRFNorm(isIncluded,:));
targetDimSpdfExRFNormSub = (spdfInfo.targetDimSpdfExRFNorm(isIncluded,:));
exitFixationSpdfInRFNormSub = (spdfInfo.exitFixationSpdfInRFNorm(isIncluded,:));
exitFixationSpdfExRFNormSub = (spdfInfo.exitFixationSpdfExRFNorm(isIncluded,:));

fprintf('\t%s: %d units\n', name, sum(isIncluded));

%     plotFileName = sprintf('%s/allSessions-%s-meanSpdfs3-v%d.png', summaryDataDir, subdivision, v);
%     fprintf('Saving to %s...\n', plotFileName);
%     
%     quickSpdfAllEvents3InARowPopMean(cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
%             arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
%             targetDimSpdfExRFNormSub, cueOnsetT, arrayOnsetT, targetDimT, plotFileName);

plotFileName = sprintf('%s/allSessions-%s-meanSpdfs5-v%d.png', summaryDataDir, name, v);
fprintf('Saving to %s...\n', plotFileName);

quickSpdfAllEvents5InARowPopMean(enterFixationSpdfInRFNormSub, enterFixationSpdfExRFNormSub, ...
        cueOnsetSpdfInRFNormSub, cueOnsetSpdfExRFNormSub, ...
        arrayOnsetHoldSpdfInRFNormSub, arrayOnsetHoldSpdfExRFNormSub, targetDimSpdfInRFNormSub, ...
        targetDimSpdfExRFNormSub, exitFixationSpdfInRFNormSub, exitFixationSpdfExRFNormSub, ...
        enterFixationT, cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, plotFileName);    
