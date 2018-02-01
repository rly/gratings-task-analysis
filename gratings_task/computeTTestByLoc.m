function tTestStatsByLoc = computeTTestByLoc(group1, group2ByLoc)
% check assumptions
% test difference in means
% assumes normal distribution with equal variance
% does not test difference in variances (see F-test of equality of
% variances)

nLoc = numel(group2ByLoc);
tTestStatsByLoc = struct();
for i = 1:nLoc
    if ~isempty(group2ByLoc{i})
        [~,tTestStatsByLoc(i).p,tTestStatsByLoc(i).ci,stats] = ttest2(group1, group2ByLoc{i});
        fn = fieldnames(stats);
        for j = 1:numel(fn)
            tTestStatsByLoc(i).(fn{j}) = stats.(fn{j});
        end
    end
    % else the fields are left empty
end