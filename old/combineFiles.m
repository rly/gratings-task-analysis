% combine files

% 1105 g3-g6
% 1104 g3-g6
% 1102 g2-g9
% 1101 g2-g3
% 0127 g2-g5

clear;
sessionName = '20170127';
startFile = 2;
nonFirstFileRange = 3:5;
cd(['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName '\']);
% gStart = load(sprintf('%s-g%d-fp.mat', sessionName, startFile));
gStart = load(sprintf('%s-g%d.mat', sessionName, startFile));

fns = fieldnames(gStart);
gMerge = gStart;
for j = nonFirstFileRange
    fprintf('Loading g%d...\n', j);
%     gj = load(sprintf('%s-g%d-fp.mat', sessionName, j));
    gj = load(sprintf('%s-g%d.mat', sessionName, j));
    offset = gMerge.Stop(end) + 1000; % add 1000 seconds
    
    for i = 1:numel(fns)
        if strcmp(fns{i}(1:2), 'FP')
            if ~isempty(strfind(fns{i}, '_ts_step'))
                assert(gStart.(fns{i}) == gj.(fns{i}));
                gNew.(fns{i}) = gj.(fns{i});
            elseif ~isempty(strfind(fns{i}, '_ts'))
                gNew.(fns{i}) = [gMerge.(fns{i}); (gj.(fns{i}) + offset)];
            elseif ~isempty(strfind(fns{i}, '_ind'))
                % add index
                gNew.(fns{i}) = [gMerge.(fns{i}); (gj.(fns{i}) + numel(gMerge.(fns{i}(1:5))))];  
            else
                % append
                gNew.(fns{i}) = [gMerge.(fns{i}); gj.(fns{i})];
            end
        elseif strcmp(fns{i}(1:3), 'SPK') && isfield(gj, fns{i})
            % ignore any cells that appear later
            % comment out below to skip all cells
%             gNew.(fns{i}) = [gMerge.(fns{i}); (gj.(fns{i}) + offset)];
        elseif strcmp(fns{i}(1:3), 'EVT')
            % add g3 stop and append
            gNew.(fns{i}) = [gMerge.(fns{i}); (gj.(fns{i}) + offset)];
        end
    end
    gNew.Start = [gMerge.Start (gj.Start + offset)];
    gNew.Stop = [gMerge.Stop (gj.Stop + offset)];
    gMerge = gNew;
end

save([sessionName '-gmerged-fp.mat'], '-struct', 'gMerge');
