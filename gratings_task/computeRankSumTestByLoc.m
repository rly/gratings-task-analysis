function rankSumStatsByLoc = computeRankSumTestByLoc(group1, group2ByLoc)
% check assumptions
% test difference in medians
% does not test difference in variance (see levene's test)

nLoc = numel(group2ByLoc);
rankSumStatsByLoc = struct();
for i = 1:nLoc
    if ~isempty(group2ByLoc{i})
        [rankSumStatsByLoc(i).p,~,stats] = ranksum(group1, group2ByLoc{i});
        fn = fieldnames(stats);
        for j = 1:numel(fn)
            rankSumStatsByLoc(i).(fn{j}) = stats.(fn{j});
        end
    end
    % else the fields are left empty
end