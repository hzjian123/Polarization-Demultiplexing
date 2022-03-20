function RDEDP_plot_learningCurve(fig, err_x, err_y)

if ischar(fig)
    fig = figure('Name', fig);
else
    set(0, 'CurrentFigure', fig);
end

for i = 1:size(err_x, 2)
    lgd{i} = ['it. #', num2str(i)];
end

subplot(121);
plot(err_x, '.');
ylim([0 2])
title('X Pol.');
legend(lgd{:})
subplot(122);
plot(err_y, '.');
ylim([0 2])
title('Y Pol.');
legend(lgd{:})
