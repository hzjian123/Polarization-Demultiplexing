function [Ex_out, Ey_out] = RDEDP_impl(Ex_in, Ey_in, opt)
%Ex_in=rxt; Ey_in=ryt;
% Constant Modulus Algorithm for 2-D modulation format
% Arguments
%     Ex_in (:,1) double
%     Ey_in (:,1) double
%     opt.SaPerSym  (1,1) double {mustBePositive, mustBeInteger} = 2
%     opt.SaOffset  (1,1) double {mustBePositive, mustBeInteger} = 1
%     opt.NumTaps   (1,1) double {mustBePositive, mustBeInteger} = 15
%     opt.StepSize  (1,1) double {mustBePositive} = 1e-3
%     opt.AidedData (:,2) double = []
%     opt.Overheads (2,1) double = [1e3; 1e4]
%     opt.NumIter   (1,1) double {mustBePositive, mustBeInteger} = 1
%     opt.RadGain   (1,1) double {mustBePositive} = 1
%     opt.Const     (:,1) double = abs(qammod(0:15,16)/sqrt(10)).^2
%     opt.InitWt    (:,4) double
% 
%     opt.Name      char = 'RDEDP'
%     opt.Active    (1,1) logical = true
%     opt.Verbose   (1,1) logical = true
%     opt.Plot      (1,1) logical = true

opt_.SaPerSym = 2;
opt_.SaOffset = 1;
opt_.NumTaps  = 9;
opt_.StepSize = 10^-3;
opt_.AidedData = [];
opt_.Overheads = [10^3; 10^4];
opt_.NumIter = 3;
opt_.RadGain = 1;
opt_.Const   = abs(qammod(0:3,4)/sqrt(10)).^2;
opt_.Name    = 'RDEDP';
opt_.Active  = true;
opt_.Verbose = true;
opt_.Plot    = true;

opt = mergestruct(opt, opt_);

if ~isfield(opt, 'InitWt')
    spike = @(l) [zeros(floor(l/2),1); 1; zeros(floor(l/2),1)];
    opt.InitWt = [spike(opt.NumTaps), ...
                  zeros(opt.NumTaps,1), ...
                  zeros(opt.NumTaps,1), ...
                  spike(opt.NumTaps)];
end

% shorten var. names
sps   = opt_.SaPerSym;
so    = opt_.SaOffset;
L     = opt_.NumTaps;
mu    = opt_.StepSize;
NI    = opt_.NumIter;
W     = opt.InitWt;
ref   = opt_.AidedData;
ad    = opt_.AidedData;
oh    = opt_.Overheads;
rg    = opt_.RadGain;
const = opt_.Const;

% helper functions
logger  = @(varargin) opt.Verbose && (exist('Log', 'file') && Log.dbg(varargin{:}) || ~exist('Log', 'file') && fprintf(varargin{:}));
normPwr = @(x) x/sqrt(mean(abs(x).^2));
indref  = @(x, r) subsref(x, struct('type', '()', 'subs', {{r}}));
replike = @(x, y) indref(repmat(x, ceil(numel(y)/numel(x)),1), 1:numel(y));

if ~opt.Active
    Ex_out = Ex_in;
    Ey_out = Ey_in;
    return;
end

N = numel(Ex_in);
T = floor(L/2);

h_xx = zeros(L,N+1);
h_xy = zeros(L,N+1);
h_yx = zeros(L,N+1);
h_yy = zeros(L,N+1);

h_xx(:,end) = W(:,1);
h_xy(:,end) = W(:,2);
h_yx(:,end) = W(:,3);
h_yy(:,end) = W(:,4);

% reference radius
R = unique(const) / mean(const) * rg;

% normalize input
Ex_in  = normPwr(Ex_in);
Ey_in  = normPwr(Ey_in);
Ex_out = zeros(size(Ex_in));
Ey_out = zeros(size(Ey_in));

% normalize aided data to unit power

if ~isempty(ref)
    trmask = zeros(numel(ref(:,1)) * sps, 1);
    logger('data aided\n');
    % derive training mask
    if numel(oh) == 2
        % periodic training
        trmask = replike([ones(oh(1),1); zeros(oh(2)-oh(1),1)], ref(:,1));
    else
        % one-time training
        trmask = [ones(oh, 1); zeros(size(ref,1)-oh, 1)];
    end
    % zero-order interp.
    refx   = normPwr(ref(:,1));
    refy   = normPwr(ref(:,2));
    R      = unique(abs(refx).^2);
    refx   = reshape(repmat(refx.', sps, 1), [], 1);
    refy   = reshape(repmat(refy.', sps, 1), [], 1);
    trmask = reshape(repmat(trmask.', sps, 1), [], 1);
end

% wrap input with zeros to 'keep' convolution
Ex_in = [zeros(T,1); Ex_in; zeros(T,1)];%T: half length of tap
Ey_in = [zeros(T,1); Ey_in; zeros(T,1)];

eps_x = zeros(N, NI);%N:Input length, NI: Iteration number
eps_y = zeros(N, NI);

for it = 1:NI
    % inherits the lastest taps
    h_xx(:,1) = h_xx(:,end);
    h_xy(:,1) = h_xy(:,end);
    h_yx(:,1) = h_yx(:,end);
    h_yy(:,1) = h_yy(:,end);

    % adaptive learning taps
    for n = 1:N
        x_in  = Ex_in(n+L-1:-1:n); % shape: [L, 1] L: Tap number
        y_in  = Ey_in(n+L-1:-1:n); % shape: [L, 1]
        x_out = h_xx(:,n).' * x_in + h_xy(:,n).' * y_in; % shape: [1, 1]
        y_out = h_yy(:,n).' * y_in + h_yx(:,n).' * x_in; % shape: [1, 1]
        % radius decision
       if isempty(ref) || ~isempty(ref) && ~trmask(n)
        % unsupervised
            a1 = sqrt(R)*x_out/abs(x_out)-x_out;
            [~, idx1] = min(abs(a1));
            rx = R(idx1);
            a2 = sqrt(R)*y_out/abs(y_out)-y_out;
            [~, idx2] = min(abs(a2));
            ry = R(idx2);
        else
        % supervised
            rx = abs(refx(n))^2;
            ry = abs(refy(n))^2;
        end
%rx=1/2;ry=1/2;
        eps_x(n, it) = rx - abs(x_out)^2;%Error tern
        eps_y(n, it) = ry - abs(y_out)^2;

        if mod(n, sps) == mod(so, sps) % update filter taps every 'sps' samples at so+sps*x position
            h_xx(:,n+1) = h_xx(:,n) + mu * eps_x(n,it) * x_out * conj(x_in);
            h_xy(:,n+1) = h_xy(:,n) + mu * eps_x(n,it) * x_out * conj(y_in);
            h_yx(:,n+1) = h_yx(:,n) + mu * eps_y(n,it) * y_out * conj(x_in);
            h_yy(:,n+1) = h_yy(:,n) + mu * eps_y(n,it) * y_out * conj(y_in);
        else
            h_xx(:,n+1) = h_xx(:,n);
            h_xy(:,n+1) = h_xy(:,n);
            h_yx(:,n+1) = h_yx(:,n);
            h_yy(:,n+1) = h_yy(:,n);
        end
        Ex_out(n) = x_out;
        Ey_out(n) = y_out;
    end
end

err_x = eps_x(so:sps:end, :).^2;
err_y = eps_y(so:sps:end, :).^2;

if 0%exist('Viewer')
    viewer = Viewer(opt.Name);
    pltfun(1).hfun = @RDEDP_plot_learningCurve;
    pltfun(1).args = {'RDEDP - learning curve', err_x, err_y};
    pltfun(2).hfun = @RDEDP_plot_centerTaps;
    pltfun(2).args = {'RDEDP - |h|', h_xx, h_xy, h_yx, h_yy, trmask};
    pltfun(3).hfun = @RDEDP_plot_timeRes;
    pltfun(3).args = {'RDEDP - time response', h_xx, h_xy, h_yx, h_yy};
    pltfun(4).hfun = @RDEDP_plot_freqRes;
    pltfun(4).args = {'RDEDP - freqency response', h_xx, h_xy, h_yx, h_yy};
    viewer.load(pltfun).register();
elseif opt.Plot
    %RDEDP_plot_learningCurve('RDEDP - learning curve', err_x, err_y);
   %RDEDP_plot_timeRes('RDEDP - time response', h_xx, h_xy, h_yx, h_yy);
    %RDEDP_plot_centerTaps('RDEDP - center tapsa', h_xx, h_xy, h_yx, h_yy);
    %RDEDP_plot_freqRes('RDEDP - freqency response',  h_xx, h_xy, h_yx, h_yy);
    %RDEDP_plot_BER('RDEDP - BER',err_x,err_y)
end

