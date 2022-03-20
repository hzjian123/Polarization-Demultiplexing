function RDEDP_plot_timeRes(fig, h_xx, h_xy, h_yx, h_yy)

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

        set(txt, 'String', num2str(k));
        subplot(221);
        plot(-T:T, abs(h_xx(:,k)), '-*');
        title('h_{xx}')
        xlim([-T,T]);
        subplot(222);
        plot(-T:T, abs(h_xy(:,k)), '-*');
        title('h_{xy}')
        xlim([-T,T]);
        subplot(223);
        plot(-T:T, abs(h_yx(:,k)), '-*');
        title('h_{yx}')
        xlim([-T,T]);
        subplot(224);
        plot(-T:T, abs(h_yy(:,k)), '-*');
        title('h_{yy}')
        xlim([-T,T]);
    end
end
