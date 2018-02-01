% plot the extent of the RF mapping stimuli
% lit says smallest RF is 1.5 degrees in radius (at 5 deg ecc) and
% largest is 5 degrees in radius (at 10 deg ecc). 
% (5,1.5), (10,5)
% m=3.5/5, b=-2. radius = mx+b

% 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
distsToFix = [180 240 320];
% max 21 possible, coded in events 11-16: 000001, 0010, ..., in the order
% below, separately for each distToFix
polarAngles = {[-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14};

nPtsDrawPerCircle = 50;
fixationBoundsRadius = 56; % <-- need to reduce this...

fixW = 1.8 * 28;
fixH = fixW;

figure_tr_inch(9, 9);
hold on;
numCircles = 0;
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles{k})
        % diameter = (mx+b) * 0.9 just scales circles smaller than RF
        diameter = (3.5 / 5 * distsToFix(k) - 2 * 28)*0.85; % pixels
        theta = linspace(0, 2 * pi, nPtsDrawPerCircle);
        rho = ones(1, nPtsDrawPerCircle) * diameter / 2;
        [centerX,centerY] = pol2cart(polarAngles{k}(l), distsToFix(k));
        [X,Y] = pol2cart(theta, rho);
        X = X + centerX;
        Y = Y + centerY;
        plot(X, Y, 'Color', zeros(3,1));
        numCircles = numCircles + 1;
    end
end

rectangle('Position', [-fixW/2 -fixH/2 fixW fixH], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
coords = circle([0 0], fixationBoundsRadius, nPtsDrawPerCircle);
plot(coords(:,1), coords(:,2), '--', 'Color', [1 1 1]);
xlim([-500 500]);
ylim([-500 500]);
axis square;
set(gca, 'Color', 127/255*ones(3,1));
fprintf('%d circles drawn.\n', numCircles);

%%
% plot the extent of the RF mapping stimuli
% lit says smallest RF is 1.5 degrees in radius (at 5 deg ecc) and
% largest is 5 degrees in radius (at 10 deg ecc). 
% (5,1.5), (10,5)
% m=3.5/5, b=-2. radius = mx+b

% 4 possible, coded in events 9-10: 01, 01, 10, 11, in the order below
distsToFix = [90 110 140 180 240 320];
% max 21 possible, coded in events 11-16: 000001, 0010, ..., in the order
% below, separately for each distToFix
polarAngles = {[-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14, [-11 -8:8 11]*pi/14};

nPtsDrawPerCircle = 50;
fixationBoundsRadius = 56; % <-- need to reduce this...

fixW = 1.8 * 28;
fixH = fixW;

figure_tr_inch(9, 9);
hold on;
numCircles = 0;
for k = 1:numel(distsToFix)
    for l = 1:numel(polarAngles{k})
        % diameter = (mx+b) * 0.9 just scales circles smaller than RF
        diameter = (3.5 / 5 * distsToFix(k) - 2 * 28)*0.85; % pixels
        if k < 3
            diameter = (3.5 / 5 * distsToFix(3) - 2 * 28)*0.85; % pixels
        end
        theta = linspace(0, 2 * pi, nPtsDrawPerCircle);
        rho = ones(1, nPtsDrawPerCircle) * diameter / 2;
        [centerX,centerY] = pol2cart(polarAngles{k}(l), distsToFix(k));
        [X,Y] = pol2cart(theta, rho);
        X = X + centerX;
        Y = Y + centerY;
        plot(X, Y, 'Color', zeros(3,1));
        numCircles = numCircles + 1;
    end
end

rectangle('Position', [-fixW/2 -fixH/2 fixW fixH], 'FaceColor', [1 1 1], 'EdgeColor', [1 1 1]);
coords = circle([0 0], fixationBoundsRadius, nPtsDrawPerCircle);
plot(coords(:,1), coords(:,2), '--', 'Color', [1 1 1]);
xlim([-500 500]);
ylim([-500 500]);
axis square;
set(gca, 'Color', 127/255*ones(3,1));
fprintf('%d circles drawn.\n', numCircles);
