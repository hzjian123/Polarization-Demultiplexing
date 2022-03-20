function RDEDP_plot_freqRes(fig, h_xx, h_xy, h_yx, h_yy)

    if ischar(fig)
        fig = figure('Name', fig);
    else
        set(0, 'CurrentFigure', fig);
    end

    N = size(h_xx, 2) - 1;
    T = (size(h_xx, 1) - 1) / 2;

    sld = uicontrol('Parent',fig,'Style','slider', 'Position',[70,5,400,15],...
                'value', N, 'min', 1, 'max', N);
    txt = uicontrol('Parent',fig, 'Style','text', 'Position',[480 5 100 15], ...
                'String', num2str(N));
    sld.Callback = @(es, ed) plot_(fig, txt, es.Value);
    plot_(fig, txt, N);

    function plot_(fig, txt, k)
        set(0, 'CurrentFigure', fig);

        k = round(k);
        if k == 0
            k = 1;
        end

        if k > N
            k = N;
        end

        ampres = @(x) abs(fftshift(fft([x; zeros(1e4-numel(x),1)])));

        set(txt, 'String', num2str(k));
        subplot(221);
        plot(ampres(h_xx(:,k)));
        title('|H_{xx}|')
        subplot(222);
        plot(ampres(h_xy(:,k)));
        title('|H_{xy}|')
        subplot(223);
        plot(ampres(h_yx(:,k)));
        title('|H_{yx}|')
        subplot(224);
        plot(ampres(h_yy(:,k)));
        title('|H_{yy}|')
    end
end
