% Remove all generated decoded stimulus events
% You may also want to those files from the example data directory sessions
% as well 

clear
fprintf('\n--------------------------------------------------------\n')

origDir = pwd;

[dataDir, sessionNames] = getUsableSessionNames();

% clean session dirs one by one
for j = 1:numel(sessionNames)
    sessionName = sessionNames{j};
    sessionDirPath = [dataDir filesep sessionName filesep];

    fprintf('Removing generated .mat files for session %s (%d/%d)...\n',...
        sessionName, j, numel(sessionNames));
    
    % print and delete all .mat files in the session directory
    allMatFiles = dir([sessionDirPath '*.mat']);
    for k = 1:numel(allMatFiles)
        matFileToDelete = [sessionDirPath allMatFiles(k).name];
        fprintf('\tDeleting %s\n', matFileToDelete);
        delete(matFileToDelete);
    end
    
    fprintf(' done\n');
end

cd(origDir);