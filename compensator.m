function[exo,eyo]=compensator(exi,eyi)
refx=1;refy=1;
%exi=fft(exi);eyi=fft(eyi);
%p=sum(x.*conj(x))/length(x)% power in time domain
%pow=sum(exi.*conj(exi))/length(x)^2 %power in frequency domain
mu=0.001;%update gain
N=length(exi);
normPwr = @(x) x/sqrt(mean(abs(x).^2));
NI=30;
L=11;NumTaps=L;%Num of tap
exi  = normPwr(exi);
eyi  = normPwr(eyi);
exo = zeros(size(exi));
eyo = zeros(size(eyi));
spike = @(l) [zeros(floor(l/2),1); 1; zeros(floor(l/2),1)];
W=[spike(NumTaps), ...
                  zeros(NumTaps,1), ...
                  zeros(NumTaps,1), ...
                  spike(NumTaps)];
T = floor(L/2);
h_xx = zeros(L,N+1);
h_xy = zeros(L,N+1);
h_yx = zeros(L,N+1);
h_yy = zeros(L,N+1);
h_xx(:,end) = W(:,1);
h_xy(:,end) = W(:,2);
h_yx(:,end) = W(:,3);
h_yy(:,end) = W(:,4);
eps_x = zeros(N, NI);%N:Input length, NI: Iteration number
eps_y = zeros(N, NI);
exi = [zeros(T,1); exi; zeros(T,1)];%T: half length of tap
eyi = [zeros(T,1); eyi; zeros(T,1)];
for it = 1:NI
    % inherits the lastest taps
    h_xx(:,1) = h_xx(:,end);
    h_xy(:,1) = h_xy(:,end);
    h_yx(:,1) = h_yx(:,end);
    h_yy(:,1) = h_yy(:,end);

    % adaptive learning taps
    for n = 1:N
        xi  = exi(n+L-1:-1:n); % shape: [L, 1] L: Tap number
        yi  = eyi(n+L-1:-1:n); % shape: [L, 1]
        xo = h_xx(:,n).' * xi + h_xy(:,n).' * yi; % shape: [1, 1]
        yo = h_yy(:,n).' * yi + h_yx(:,n).' * xi; % shape: [1, 1]
         rx = abs(refx)^2;
         ry = abs(refy)^2;
        eps_x(n, it) = rx - abs(xo)^2;%Error tern
        eps_y(n, it) = ry - abs(yo)^2;
            h_xx(:,n+1) = h_xx(:,n) + mu * eps_x(n,it) * xo * conj(xi);
            h_xy(:,n+1) = h_xy(:,n) + mu * eps_x(n,it) * xo * conj(yi);
            h_yx(:,n+1) = h_yx(:,n) + mu * eps_y(n,it) * yo * conj(xi);
            h_yy(:,n+1) = h_yy(:,n) + mu * eps_y(n,it) * yo * conj(yi);
            exo(n)=xo;
            eyo(n)=yo;
%exo=ifft(exo);eyo=ifft(eyo);
    end
end
end