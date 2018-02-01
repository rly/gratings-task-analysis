ENV = getEnv();
dataDir = ENV.dataDir;

dataDirContents = dir(dataDir);

for k = length(dataDirContents):-1:1
    % remove non-directories and hidden directories and directories that do
    % not start with C or L
    if ~dataDirContents(k).isdir || dataDirContents(k).name(1) == '.' || ...
            ~(dataDirContents(k).name(1) == 'C' || dataDirContents(k).name(1) == 'L')
        dataDirContents(k) = [];
        continue;
    end
end

hasAD = nan(length(dataDirContents), 4);
for i = 1:length(dataDirContents)
    fprintf('Inspecting %s...\n', dataDirContents(i).name);
    origLfpFileName = [dataDir filesep dataDirContents(i).name filesep ...
            'lfps' filesep dataDirContents(i).name '_lfps.mat'];
    if exist(origLfpFileName, 'file')
        dataDirContents(i).lfpFileContents = who('-file', origLfpFileName);
        % TODO check for LFPs and events
        hasAD(i,1) = any(strcmp('AD01', dataDirContents(i).lfpFileContents));
        hasAD(i,2) = any(strcmp('AD02', dataDirContents(i).lfpFileContents));
        hasAD(i,3) = any(strcmp('AD03', dataDirContents(i).lfpFileContents));
        hasAD(i,4) = any(strcmp('AD04', dataDirContents(i).lfpFileContents));
    else
        dataDirContents(i).lfpFileContents = {};
    end
    
    spikeDirContents = dir([dataDir filesep dataDirContents(i).name filesep ...
            'spikes']);
    for k = length(spikeDirContents):-1:1
        % remove hidden files
        if spikeDirContents(k).name(1) == '.'
            spikeDirContents(k) = [];
            continue;
        end
        % TODO check for area names and a sessionName_spikes.mat file within
    end
    dataDirContents(i).spikeDirContents = spikeDirContents;
end