%% pul units w/ sig pos cue resp at contralateral (P3) 
precondition = isInPulvinar & isSignificantCueResponseInc & isInRFP3;
fprintf('Units in the pulvinar with significantly increased cue response compared to baseline and InRF P3: %d\n', sum(precondition));

suaAllPath = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\PUL_SUA_GRATINGS_ALL';

suaFiles = dir(suaAllPath);
esFilePathSub = esFileNames(precondition);
outputPdfs = cell(numel(esFilePathSub), 1);
for i = 1:numel(esFilePathSub)
    fprintf('Processing file %d/%d = %d%%\n', i, numel(esFilePathSub), round(100*i/numel(esFilePathSub)));
    [~,esFileNameNoExt] = fileparts(esFilePathSub{i});
    suaPath = dir(sprintf('%s/%s.mat', suaAllPath, esFileNameNoExt));
    assert(numel(suaPath) == 1);
    match = regexp(suaPath(1).name, sprintf('(.*)-evokedSpiking-v%d.mat', v), 'tokens');
    assert(numel(match{1}) == 1);
    pngPath = sprintf('%s/%s-visual-v%d.png', suaPath(1).folder, match{1}{1}, v);
    assert(exist(pngPath, 'file') > 0);
    outputPdfs{i} = quickPngToPdf(pngPath);
end

outputFileName = sprintf('%s/pulCueIncP3-v%d.pdf', suaAllPath, v);
fprintf('Printing to PDF: %s...\n', outputFileName);
deleteAndAppendPdfs(outputFileName, outputPdfs);
