function RDEDP_plot_centerTaps(fig, h_xx, h_xy, h_yx, h_yy, trmask)

if ischar(fig)
    fig = figure('Name', fig);
else
    set(0, 'CurrentFigure', fig);
end

T = (size(h_xx, 1) - 1) / 2;
h_xx = h_xx(T+1,:);
h_yx = h_yx(T+1,:);
h_xy = h_xy(T+1,:);
h_yy = h_yy(T+1,:);

lgdoff = @(h) set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

untrind = find(trmask==0);
trind = find(trmask==1);

hold on;
plot(untrind, abs(h_xx(untrind)));
plot(untrind, abs(h_yx(untrind)));
plot(untrind, abs(h_xy(untrind)));
plot(untrind, abs(h_yy(untrind)));
legend('h_{xx}', 'h_{yx}', 'h_{xy}', 'h_{yy}');
if numel(trind) > 0
    p(1) = plot(trind, abs(h_xx(trind)), 'LineWidth', 6);
    p(2) = plot(trind, abs(h_yx(trind)), 'LineWidth', 6);
    p(3) = plot(trind, abs(h_xy(trind)), 'LineWidth', 6);
    p(4) = plot(trind, abs(h_yy(trind)), 'LineWidth', 6);
    arrayfun(lgdoff, p);
end
hold off;
grid on;

