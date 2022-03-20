function RDEDP_plot_BER(fig,err_x,err_y)
if ischar(fig)
    fig = figure('Name', fig);
else
    set(0, 'CurrentFigure', fig);
end
subplot(121);
plot(err_x(:,1), '.');
title('X Pol.');
subplot(122);
plot(err_y(:,1), '.');
title('Y Pol.');
