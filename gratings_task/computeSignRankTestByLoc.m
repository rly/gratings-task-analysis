function signRankStatsByLoc = computeSignRankTestByLoc(group1ByLoc, group2ByLoc)
% check assumptions
% test difference in medians
% does not test difference in variance (see levene's test)

nLoc = numel(group2ByLoc);
assert(nLoc == numel(group1ByLoc));
signRankStatsByLoc = struct();
for i = 1:nLoc
    if ~isempty(group1ByLoc{i}) && ~isempty(group2ByLoc{i})
        assert(numel(group1ByLoc{i}) == numel(group2ByLoc{i}));
        [signRankStatsByLoc(i).p,~,stats] = signrank(group1ByLoc{i} - group2ByLoc{i});
        fn = fieldnames(stats);
        for j = 1:numel(fn)
            signRankStatsByLoc(i).(fn{j}) = stats.(fn{j});
        end
    else
        signRankStatsByLoc(i).p = NaN;
    end
end