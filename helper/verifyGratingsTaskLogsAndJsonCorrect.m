function verifyGratingsTaskLogsAndJsonCorrect(trialParamsAllCorrect, trialStructs)

trialResults = cellfun(@(x) x.response.type, trialStructs, 'UniformOutput', false);
isCorrect = strcmp(trialResults, 'correct-response');
trialStructsCorrect = trialStructs(isCorrect);

assert(size(trialParamsAllCorrect, 1) == numel(trialStructsCorrect));
% assert(numel(rt) == numel(trialStructsCorrect));

rtTol = 0.025; % seconds, at least for comparison between RT computed using time of target dim to time of juice onset 
for i = 1:size(trialParamsAllCorrect, 1)
    assert(trialParamsAllCorrect(i,1) == trialStructsCorrect{i}.numTrialsCorrectNoRpt);
    assert(trialParamsAllCorrect(i,2) == trialStructsCorrect{i}.trialNum - 1);
    assert(trialParamsAllCorrect(i,3) == ~trialStructsCorrect{i}.isRedo);
    assert(trialParamsAllCorrect(i,4) == trialStructsCorrect{i}.cueLoc);
    assert(trialParamsAllCorrect(i,5) == round(trialStructsCorrect{i}.cueTheta * 4/pi));
    assert(trialParamsAllCorrect(i,6) == trialStructsCorrect{i}.cueArrayDelayDuration);
    if trialParamsAllCorrect(i,7) == -1
        assert(15000 == trialStructsCorrect{i}.holdShapeDuration);
    else
        assert(trialParamsAllCorrect(i,7) == trialStructsCorrect{i}.holdShapeDuration);
    end
%     if trialStructsCorrect{i}.holdShapeDuration >= 15000
%         assert(abs(rt(i) - trialStructsCorrect{i}.response.responseTime / 1000) < rtTol);
%     else
%         assert(abs(rt(i) - (trialStructsCorrect{i}.response.responseTime - trialStructsCorrect{i}.holdShapeDuration) / 1000) < rtTol);
%     end
end
