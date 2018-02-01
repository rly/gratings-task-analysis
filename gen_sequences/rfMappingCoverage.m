% plot the extent of the RF mapping stimuli

% 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
distsToFix = 150:60:330;
% 7 possible, coded in events 11-14: 0001, 0010, ..., 1011, in the order below
polarAngles = (-5:5)*pi/8;
% 4 possible, coded in events 15-16: 00, 01, 10, 11, in the order below
gratingAngles = (0:3)*pi/4;

nPtsDrawPerCircle = 50;
fixationBoundsRadius = 68;

fixW = 1.8 * 28;
fixH = fixW;

figure_tr_inch(6, 6);
hold on;
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles)
        diameter = distsToFix(k) / 150 * 3 * 28; % pixels
        theta = linspace(0, 2 * pi, nPtsDrawPerCircle);
        rho = ones(1, nPtsDrawPerCircle) * diameter / 2;
        [centerX,centerY] = pol2cart(polarAngles(l), distsToFix(k));
        [X,Y] = pol2cart(theta, rho);
        X = X + centerX;
        Y = Y + centerY;
        plot(X, Y, 'Color', zeros(3,1));
    end
end

rectangle('Position', [-fixW/2 -fixH/2 fixW fixH], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
coords = circle([0 0], fixationBoundsRadius, nPtsDrawPerCircle);
plot(coords(:,1), coords(:,2), '--', 'Color', [1 1 1]);
xlim([-500 500]);
ylim([-500 500]);
axis square;
set(gca, 'Color', 127/255*ones(3,1));
