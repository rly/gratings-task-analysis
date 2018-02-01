function [dataDir, sessionNames, goodSessionNames, hasAD, sessionsAndAD, hasLfpSessionNames, hasSpikesSessionNames] = getUsableSessionNames()
% this function should yield the 60 sessions with a non-empty LFP data
% directory and a non-empty spike data directory, plus L110810 which should
% be handled specially

ENV = getEnv();
dataDir = ENV.dataDir;

dataDirContents = dir(dataDir);

skipSessions = {...
    'C110519'; % no events
    'C110520'; % no events
    'C110523'; % no events
    'C110525'; % no events
    'C110526'; % no events
    'C110606'; % no events
    'C110621'; % no events
    'C110628'; % no events
    'C110706'; % no events
    'C110707'; % no events
    'C110719'; % unusual data
    'C110729'; % no events
    'L100913'; % missing AD01, presentation logs
    'L101004'; % no events
    'L101005'; % no events
    'L101006'; % no events
%     'L110504'; % missing AD01
};

for k = length(dataDirContents):-1:1
    % remove non-directories and hidden directories and directories that do
    % not start with C or L
    if ~dataDirContents(k).isdir || ...
            dataDirContents(k).name(1) == '.' || ...
            ~(dataDirContents(k).name(1) == 'C' || ...
              dataDirContents(k).name(1) == 'L')
        dataDirContents(k) = [];
    end
end
% these are all the sessions (102) based on an ls of klab_data
sessionNames = {dataDirContents.name};
hasLfpSessionNames = {};
hasSpikesSessionNames = {};
% 'L110810' has spike data but no LFP data

numSessions = numel(sessionNames);
goodSessionNames = {};
sessionsAndAD = cell(0,2);
hasAD = nan(numSessions, 4);

% process each session one by one, and remove bad ones
for j = 1:numSessions
	sessionName = sessionNames{j};
    
    % make sure the relevant dirs and files exist
    % L110810 is special case -- no lfp data but want to process spike data
    lfpDir = [dataDir sessionName filesep 'lfps'];
    spikeDir = [dataDir sessionName filesep 'spikes'];
    
    skipReason = '';
    hasLfp = 1;
    hasSpikes = 1;
    if find(strcmp(skipSessions,sessionName));
        fprintf('Session %s (%d/%d): %s\n', sessionName, ...
                j, numSessions, 'Skipped by exclusion');
        continue;
    end
    if ~isdir(lfpDir)
        skipReason = [skipReason 'no LFP dir, '];
        hasLfp = 0;
    end
    if numel(dir(lfpDir)) == 2
        skipReason = [skipReason 'LFP dir is empty, '];
        hasLfp = 0;
    end
    if ~isdir(spikeDir)
        skipReason = [skipReason 'no spike dir, '];
        hasSpikes = 0;
    end
    if numel(dir(spikeDir)) == 2
        skipReason = [skipReason 'spike dir is empty, '];
        hasSpikes = 0;
    end
    if strcmp(sessionName, 'L110810')
        skipReason = [];
        hasLfp = 0;
        fprintf('Including L110810, but note that there is no LFP data\n');
    end
    if ~isempty(skipReason)
        fprintf('Session %s (%d/%d): %s\n', sessionName, ...
                j, numSessions, skipReason(1:end-2));
    else
        goodSessionNames = [goodSessionNames sessionName];
    end
    if hasLfp
        hasLfpSessionNames = [hasLfpSessionNames sessionName];
        origLfpFileName = [lfpDir filesep sessionName '_lfps.mat'];
        if exist(origLfpFileName, 'file')
            dataDirContents(j).lfpFileContents = who('-file', origLfpFileName);
            hasAD(j,1) = any(strcmp('AD01', dataDirContents(j).lfpFileContents));
            hasAD(j,2) = any(strcmp('AD02', dataDirContents(j).lfpFileContents));
            hasAD(j,3) = any(strcmp('AD03', dataDirContents(j).lfpFileContents));
            hasAD(j,4) = any(strcmp('AD04', dataDirContents(j).lfpFileContents));
            sessionsAndAD(size(sessionsAndAD,1)+1,:) = {sessionName hasAD(j,:)};
        end
    end
    if hasSpikes
        hasSpikesSessionNames = [hasSpikesSessionNames sessionName];
    end
end
