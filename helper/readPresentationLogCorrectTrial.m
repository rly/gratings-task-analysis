function trialParams = readPresentationLogCorrectTrial(filepath, filter, logIndices)
%
% filepath: main file path for those log files
% filter: the string part of log files, like 'edit6_flanker_task_27_with_eye_tracker_12w2_and_eyeErr_'
% startnumber: the start point for the first log file (those number must be continuous), like 5.
% The file name will be [filter startnumber], like 'edit6_flanker_task_27_with_eye_tracker_12w2_and_eyeErr_5'
%
% array: the code for correct responses, corresponding to Event010 in the
% same order
%

P = pwd;
cd(filepath);
% logfile = dir('*.log');
% fileName = [];
% count = 0;
% while count < length(logfile)*5 % NOTE: clean me up
%     if exist([filter num2str(logIndices(1)+count) '.log'], 'file')
%         sortfile = char(sortfile,[filter num2str(logIndices(1)+count) '.log']);
%     end
%     count = count + 1;
% end
% sortfile(1,:) = [];

trialParams = [];

% the important bits are after two blank lines after line 5
for ifile = 1:numel(logIndices)%size(sortfile,1)
    fileName = deblank([filter num2str(logIndices(ifile)) '.log']);
    if ~exist([filter num2str(logIndices(ifile)) '.log'], 'file')
        error('Cannot find file %s', fileName);
    end
    fid = fopen(fileName, 'r+');
    % read the first 5 rows to change the point location
    L = 5;
    for irow = 1:L, fgetl(fid); end
    while 1
        tline = fgetl(fid);
        if isempty(tline)
            break;
        end
    end
    tline = fgetl(fid); % read empty line
    assert(isempty(tline));
    tline = fgetl(fid); % read header line
    assert(~isempty(tline));
    tline = fgetl(fid); % read empty line
    assert(isempty(tline));
    data = textscan(fid, '%s %s %s %d %d %d %d %d %d %d %d %s', 'delimiter', '\t');
    fclose(fid);
    
    nrow = size(data{1},1);
    for irow = 1:nrow
        if ~isempty(strfind(data{1,2}{irow}, 'correct')) && isempty(strfind(data{1,2}{irow}, 'incorrect'))
            if ~isempty(strfind(data{1,2}{irow-5}, 'begin')) || ~isempty(strfind(data{1,2}{irow-4}, 'begin'))
                if ~isempty(strfind(data{1,2}{irow-5}, 'begin'))
                    beginRow = irow - 5;
                else
                    beginRow = irow - 4;
                end
                % note to self: use regexp
                beginCell = textscan(data{1,2}{beginRow}, 'begin %d %d  new%d P%d theta%f delay%d %c%d');
                prevCorrectTrialCount = beginCell{1};
                prevTrialCount = beginCell{2};
                isNewTrial = beginCell{3};
                cueLocation = beginCell{4};
                thetaRad = beginCell{5};
                thetaPos = round(thetaRad*4/pi);
                delayTime = beginCell{6};
                trialType = beginCell{7};
                if trialType == 'h'
                    holdTime = beginCell{8};
                else
                    assert(trialType == 'r');
                    assert(isempty(beginCell{8}));
                    holdTime = -1;
                end
            else
                warning('a correct trial without a begin event: %d, %s', ...
                        irow, data{1,2}{irow-5});
            end
            
            if ~isempty(strfind(data{1,2}{beginRow+1}, 'cue on'))
                cueOnCell = textscan(data{1,2}{beginRow+1}, 'cue on %d %d  new%d P%d theta%f delay%d %c%d');
                % TODO more asserts
                assert(delayTime == cueOnCell{6});
            else
                warning('a correct trial without a cue on event: %d, %s', ...
                        irow, data{1,2}{beginRow+1});
            end
            assert(~isempty(strfind(data{1,2}{beginRow+2}, 'cue off')));
            
            if trialType == 'h'
                arrayStr = data{1,2}{irow-2};
            else % 'r'
                arrayStr = data{1,2}{irow-1};
            end
            
            if ~isempty(strfind(arrayStr, 'array on'))
                arrayOnCell = textscan(arrayStr, 'array on %d %d  new%d P%d theta%f delay%d %c%d');
                % TODO more asserts
                assert(delayTime == arrayOnCell{6});
                assert(trialType == arrayOnCell{7});
                if trialType == 'h'
                    assert(holdTime == arrayOnCell{8});
                else
                    assert(isempty(arrayOnCell{8}));
                end
            else
                warning('a correct trial without a cue off event: %d, %s, %s', ...
                        irow, data{1,2}{irow-2}, data{1,2}{irow-1});
            end
            
            correctCell = textscan(data{1,2}{irow}, 'correct %dx juice');
%             assert(trialType == correctCell{1});
            rt = double(data{1,7}(irow) - data{1,7}(irow-1)) / 10000;
            
            % append new row to trial param -- correct trials only
            trialParams = [trialParams; ...
                    double(prevCorrectTrialCount) double(prevTrialCount) ...
                    double(isNewTrial) double(cueLocation) ...
                    double(thetaPos) double(delayTime) ...
                    double(holdTime) double(rt) double(data{1,7}(irow))];
        end
    end
end

cd(P)