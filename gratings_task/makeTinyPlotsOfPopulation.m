function makeTinyPlotsOfPopulation(RInRF, RInRFErr, RExRF, RExRFErr, t, unitNames, titleBase, plotFileBaseName)

cols = lines(4);
inRFCol = [0.9 0 0];
exRFCol = [0 0 0.9];

nCells = size(RInRF, 1);
nFigures = ceil(nCells / 60);

assert(nCells == numel(unitNames));

numSessionsUsed = 0;
for j = 1:nCells
    % create another counter in case you want to skip some cells
    numSessionsUsed = numSessionsUsed + 1;
    if mod(numSessionsUsed, 60) == 1
        if numSessionsUsed > 60
            % save last figure
            plotFileName = sprintf('%s_pg%d.png', plotFileBaseName, ceil(numSessionsUsed / 60) - 1);
            if ~isempty(plotFileName)
                export_fig(plotFileName, '-nocrop');
            end
        end
        figure_tr_inch(13.5, 7); clf;
        set(gcf, 'Color', 'white');

        % make title axis
        axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
        set(get(axBig, 'Title'), 'Visible', 'on')
        % make title
        titleText = sprintf('%s - Page %d of %d', titleBase, ceil(numSessionsUsed / 60), nFigures);
        titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
        title(titleText, 'FontSize', 14, titleParams{:});
    end

    subaxis(6, 10, mod(numSessionsUsed - 1, 60) + 1, 'ML', 0.01, 'MB', 0.01, 'MT', 0.05, 'MR', 0.01, 'SH', 0.02, 'SV', 0.02, 'hold', 1);

    xBounds = [-0.25 0.25];
    hold on;

    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1));
    plot(t([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1));

    jbfill(t, ...
            RExRF(j,:) + RExRFErr(j,:), ...
            RExRF(j,:) - RExRFErr(j,:), ...
            exRFCol, exRFCol, 0.3);
    plot(t, RExRF(j,:), 'Color', exRFCol, 'LineWidth', 3);
    
    jbfill(t, ...
            RInRF(j,:) + RInRFErr(j,:), ...
            RInRF(j,:) - RInRFErr(j,:), ...
            inRFCol, inRFCol, 0.3);
    plot(t, RInRF(j,:), 'Color', inRFCol, 'LineWidth', 3);
    
    text(1, 0.95, unitNames{j}, 'HorizontalAlignment', 'right', ...
            'FontSize', 7, 'Interpreter', 'none', 'Units', 'normalized');

    xlim(xBounds);
    ylim([-1 1]);
%         set(gca, 'XTick', -0.2:0.2:0.2);
%         set(gca, 'TickLength', get(gca, 'TickLength') * 3);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca ,'XTickLabel', []);
    set(gca ,'YTickLabel', []);
    set(gca, 'box', 'off');
end

%% save
plotFileName = sprintf('%s_pg%d.png', plotFileBaseName, ceil(numSessionsUsed / 60));
if ~isempty(plotFileName)
    export_fig(plotFileName, '-nocrop');
end
