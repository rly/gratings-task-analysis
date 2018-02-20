function fig_h = figure_tr_inch(w, h, x, y)
% create a figure of size (w,h) in inches, given a dpi of 96 (which is the
% windows default). default location (x,y) = (0,0) is top-left corner of 
% the screen
if nargin < 3
    x = 0; y = 0;
end

% convert inches to pixels using windows default dpi
dpi = 96;
wp = dpi * w;
hp = dpi * h;
xp = dpi * x;
yp = dpi * y;

% for reference, matlab default on windows is 96 dpi so 7 inch = 672 pixels
% matlab default on mac/linux is 72 dpi so 7 inch = 504 pixels
% use get(0, 'ScreenPixelsPerInch') to check
scrsz = get(0, 'ScreenSize'); % get screensize in pixels

left = xp + 1; % adjust +1 for default windows border
btm = max(0, scrsz(4) - hp - yp - 79); % adjust -79 for windows title bar

fig_h = figure('Position', [left btm wp hp], 'Color', 'w');
subaxis(1, 1, 1, 'MR', 0.2/w, 'ML', 0.8/w, 'MB', 0.6/h, 'MT', 0.5/h);
