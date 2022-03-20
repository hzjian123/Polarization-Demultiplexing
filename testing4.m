% Demonstration of Eb/N0 Vs BER for QPSK (waveform simulation)
clear all;clc;
init
% include optix libraries in path
DEBUG = false; % output figures and extra information
N=5000;%Number of bits for each arm
alpha =0; Nsym =8;L=8;% RC filter alpha and filter span in symbols
rcPulse= raisedCosineFunction(alpha,L,Nsym);%;%Number of bits to tran

EbN0dB = -4:2:10; % Eb/N0 range in dB for simulation
%OF =L; %oversampling factor, sampling frequency will be fs=OF*fc
EbN0lin = 10.^(EbN0dB/10); %converting dB values to linear scale
BER = zeros(length(EbN0dB),1); %For BER values for each Eb/N0
M=dlmread('D:\Hzj\WDM.vtmu_pack\Resources\..\Inputs\data\Channel1.txt');
m=dlmread('D:\Hzj\WDM.vtmu_pack\Resources\..\Outputs\data\Channel1.txt');
IXb = M(:,1);%X polarization Q arm bit stream initialization
QXb = M(:,2);
IYb = M(:,3);
QYb = M(:,4);
a=[IXb;QXb;IYb;QYb];% Combined vector for BER computation
IXt= m(:,1);QXt=m(:,2);IYt=m(:,3);QYt=m(:,4);
sx=IXt+1i*QXt;sy=IYt+1i*QYt;
sx=sx/(mean(abs(sx).^2));
sy=sy/(mean(abs(sy).^2));
[nx,ny]=deal(zeros(length(sx),length(EbN0dB)));
for i=1:8%1:length(EbN0dB)%for each SNR
    spike = @(l) [zeros(floor(l/2),1); 1; zeros(floor(l/2),1)];
    Ebx=2*sum(abs(sx).^2)/(length(sx)); %compute energy per bit
    Eby=2*sum(abs(sy).^2)/(length(sy)); %compute energy per bit
    N0x= Ebx/EbN0lin(i); %required noise spectral density from Eb/N0
    N0y= Eby/EbN0lin(i); %required noise spectral density from Eb/N0
    
    nxr = sqrt(N0x/2)*(randn(1,length(sx)))';nxi = sqrt(N0x/2)*(randn(1,length(sx)))';%computed noise for x
    nyr = sqrt(N0y/2)*(randn(1,length(sy)))';nyi = sqrt(N0y/2)*(randn(1,length(sy)))';%computed noise for y
    nx(:,i)=nxr + 1i*nxi;
    ny(:,i)=nyr + 1i*nyi;
    rxt = sx + nx(:,i);%add noise
    ryt = sy + ny(:,i);%add noise
    %rx=conv(rx,rcPulse,'same');
    %ry=conv(ry,rcPulse,'same');
    sps=2;
    
    sxorigin=IXb+1i*QXb;
    syorigin=IYb+1i*QYb;
    % Adaptive linear equalization
    % equalize the signal power to a reasonable range
    islms=0;
    if islms==1
        % frame sync.
        dbp  = optix.DigitalBackProp('SpanLength', 1226e3/15, 'SpanNumber', 15, 'LaunchedPower', 5, 'xi_1', 0., 'xi_2', 0., 'StPS', 1, 'BlockProcMode', true, 'TimeDomainMode', false,'BlockSize', 512-128);
        %dbp  = optix.DigitalBackProp('SpanLength', 1226e3/15, 'SpanNumber', 15, 'VirtualSpanNumber', 15, 'LaunchedPower', 5, 'xi_1', 0., 'xi_2', 0., 'StPS', 1, 'BlockProcMode', false, 'TimeDomainMode', true,'NumTaps', 350);
        cma  = optix.LinearEqualizerMIMO('Algorithm','CMA','StepSize', 1e-3, 'NumTaps', 21, 'ReferenceTap',11, ...
            'TapWeightsOutputPort', true, 'ErrorOutputPort', true);
        lms  = optix.LinearEqualizerMIMO('Algorithm','LMS','StepSize',1e-3, 'NumTaps',15, 'ReferenceTap',8, ...
            'TapWeightsOutputPort', true, 'ErrorOutputPort', true, 'TrainingFlagInputPort', true,'InputSamplesPerSymbol',2,'AdaptWeights',true, 'Constellation', qammod([0:3], 4)/2);%,'InitialWeights',zeros(5,4));
        foex = optix.FreqOffsetEqualizer('FrequencyOffsetOutputPort', true, 'MagnitudeOutputPort', true);
        foey = optix.FreqOffsetEqualizer('FrequencyOffsetOutputPort', true, 'MagnitudeOutputPort', true);
        cprx = optix.CarrierPhaseSynchronizer('EstimatedPhaseOutputPort', true);
        cpry = optix.CarrierPhaseSynchronizer('EstimatedPhaseOutputPort', true);
        [refx,TX]=deal([((sxorigin)/mean(abs(sxorigin).^2));((sxorigin)/mean(abs(sxorigin).^2))]);
        [refy,TY]=deal([((syorigin)/mean(abs(syorigin).^2));((syorigin)/mean(abs(syorigin).^2))]);
        pdx  = optix.PreambleDetector('PreambleSource', 'Input port', 'DetectionMetricOutput', true, 'Threshold', -inf);
        pdy  = optix.PreambleDetector('PreambleSource', 'Input port', 'DetectionMetricOutput', true, 'Threshold', -inf);
        normpwr = @(x,v) x/sqrt(mean(abs(x).^2))*sqrt(v); % normalize 1D power to given value
        
        xcorlen = 100; % number of symbols used to do cross-correlation, a proper value should be searched to maxmize the corr. magnitude
        [indx, detmetx] = pdx(refx, rxt(1:sps:sps*xcorlen));
        [indy, detmety] = pdy(refy, ryt(1:sps:sps*xcorlen));
        [~,iix] = max(abs(detmetx));
        ix = indx(iix);
        [~,iiy] = max(abs(detmety));
        iy = indy(iiy);
        if 0%DEBUG
            ix
            iy
            figure;
            subplot(211);plot(detmetx);title('magnitude of xcorr. in frame sync. of X pol.');
            subplot(212);plot(detmety);title('magnitude of xcorr. in frame sync. of Y pol.');
        end
        refx = refx(ix-xcorlen+1:ix-xcorlen+length(rxt)/sps);
        refy = refy(iy-xcorlen+1:iy-xcorlen+length(ryt)/sps);
        
        rx_temp = rxt/sqrt(mean(abs(rxt).^2));
        ry_temp = ryt/sqrt(mean(abs(ryt).^2));
        outlen = length(refx);
        [rxt,ryt,ex,ey] = deal(ones(outlen,1));
        [wxx,wxy,wyx,wyy] = deal(ones(lms.NumTaps,outlen));
        ReferenceTap=1;
        train_period = 10;
        train_cnt = 0;
        train_flag = true;
        nff=9;stepSize=0.001;
        [w_xy,w_yx,w_xx,w_yy]=deal(zeros(nff,length(refx)));
        [w_xx(:,1),w_yy(:,1)]=deal(spike(nff));
        [e_x,e_y]=deal(zeros(1,length(refx)));
        
        rxtp=rx_temp(1:2:end);rytp=ry_temp(1:2:end);
        %rxtp = [zeros(ReferenceTap-1,1); rxtp; zeros(nff-ReferenceTap,1)];%T: half length of taprytp = [zeros(ReferenceTap-1,1); rytp; zeros(nff-ReferenceTap,1)];%T: half length of tap
        for ind = 1:length(refx)
            if ind<ceil(nff/sps)
                [u_x,u_y]=deal(zeros(nff,1));
                u_x(1:sps*ind)=rx_temp(sps*ind:-1:1);u_y(1:sps*ind)=ry_temp(sps*ind:-1:1);
            else
                u_x=rx_temp(sps*ind:-1:sps*ind-nff+1);u_y=rx_temp(sps*ind:-1:sps*ind-nff+1);
            end
            
            %if ind>1 [w_xx(:,ind),w_xy(:,ind),w_yx(:,ind),w_yy(:,ind)]=deal(w_xx(:,ind-1),w_xy(:,ind-1),w_yx(:,ind-1),w_yy(:,ind-1));end
            %if ind==1 symbolCounter=0;end
            %trainingSymbolsX=refx(ind);trainingSymbolsY=refy(ind);
            x_i=rx_temp((ind-1)*sps+1:ind*sps);y_i=ry_temp((ind-1)*sps+1:ind*sps);
            % [rxt(ind),ryt(ind),e_x(ind),e_y(ind),w_xx(:,ind),w_xy(:,ind),w_yx(:,ind),w_yy(:,ind)]=LMS(x_i,y_i,u_x,u_y,symbolCounter,...
            %w_xx(:,ind),w_xy(:,ind),w_yx(:,ind),w_yy(:,ind),trainingSymbolsX,trainingSymbolsY,train_flag,stepSize,nff);
             for cn=1:10
            [rxt(ind),ryt(ind),ex(ind),ey(ind),wxx(:,ind),wxy(:,ind),wyx(:,ind),wyy(:,ind)] = ...
                lms(rx_temp((ind-1)*sps+1:ind*sps),ry_temp((ind-1)*sps+1:ind*sps),refx(ind),refy(ind),train_flag);
            %symbolCounter=symbolCounter+1; if symbolCounter==nff symbolCounter=0; end
            train_cnt = train_cnt + 1;
            if train_cnt == train_period
                train_cnt = 0;
                train_flag = true;
            else
                train_flag = false;
            end
             end
        end
        
        if i==8
            DEBUG=1;
            % scatterplot(rxt); title('after LMS X pol.');
            %scatterplot(ryt); title('after LMS Y pol.');
            figure;
            %subplot(411); plot(abs(wxx(:,end)), '-o'); title('converged filter taps of X pol.');
            %subplot(412); plot(abs(wyy(:,end)), '-o'); title('converged filter taps of Y pol.');
            subplot(413); plot(abs(e_x)); title('learning curve of X pol.')
            subplot(414); plot(abs(e_y)); title('learning curve of Y pol.')
            subplot(411); plot(abs(ex)); title('learning curve of X pol.')
            subplot(412); plot(abs(ey)); title('learning curve of Y pol.')
            
        end
    end
    rxt=rxt(3:end);ryt=ryt(3:end);
    % refx=refx(1:end-1);refy=refy(1:end-1);
    rxt=[rxt;zeros(2,1)];ryt=[ryt;zeros(2,1)];
    
    
    
    % pfox = comm.PhaseFrequencyOffset('PhaseOffset',rand(N*2,1)*90-90/2);
    % pfoy = comm.PhaseFrequencyOffset('PhaseOffset',rand(N*2,1)*90-90/2);
    % rxt=pfox(rxt);
    % ryt=pfoy(ryt);
    opt=struct;
    %[H,r] = corrmtx(rxt,length(rxt)-1); rr=eig(r);1/max(rr);
    
    %opt.AidedData=[sxorigin(1:200),syorigin(1:200)];
    
    %if i==8
     [rxt,ryt]=RDEDP_impl(rxt,ryt,opt);%Compensator
    
    % [rxt,ryt]=PU_CMA(rxt,ryt,opt);%Compensator
    % [rxt,ryt]=RDEDP(rxt,ryt,opt);%Compensator
    
    
    
    % [rxt,ryt]=MCMA(rxt,ryt,opt);%Compensator
    
    
    %carrierSync = comm.CarrierSynchronizer;
    %carrierSync.Modulation='QPSK';
    %rxt=carrierSync(rxt);
    %ryt=carrierSync(ryt);
    
    IXrt=real(rxt);%Received IX arm with oversampling of 2
    IYrt=real(ryt);
    QXrt=imag(rxt);
    QYrt=imag(ryt);
    
    if islms==1
        IXr=IXrt;%(1:2:end);%Downsampling to 1 time
        IYr=IYrt;%(1:2:end);
        QXr=QXrt;%(1:2:end);
        QYr=QYrt;%(1:2:end);
    else
        IXr=IXrt(1:2:end);%Downsampling to 1 time
        IYr=IYrt(1:2:end);
        QXr=QXrt(1:2:end);
        QYr=QYrt(1:2:end);
    end
    rxd=IXr+1i*QXr;
    ryd=IYr+1i*QYr;
    
    IXr=real(rxd);%Received IX arm with oversampling of 2
    IYr=real(ryd);
    QXr=imag(rxd);
    QYr=imag(ryd);
    ber=zeros(1,10);
    %IXr=upfirdn(IXr,rcPulse);IXr=IXr(L*Nsym/2+1:end-L*Nsym/2);
    %IYr=upfirdn(IXr,rcPulse);IYr=IYr(L*Nsym/2+1:end-L*Nsym/2);
    %QXr=upfirdn(QXr,rcPulse);QXr=QXr(L*Nsym/2+1:end-L*Nsym/2);
    if 0%i==8
        eyediagram(IXr,2);
    end
    ac=100;%Accuracy of phase compensation
    ph=0;
    if ph==1
        for p=1:ac
            IXrt=(IXr>0)*2-1;%Convert to NRZ
            IYrt=(IYr>0)*2-1;
            QXrt=(QXr>0)*2-1;
            QYrt=(QYr>0)*2-1;
            a_cap=[IXrt;IYrt;QXrt;QYrt];%Combined signal for BER computation
            ber(p)= sum(a~=a_cap)/length(a);%Bit Error Rate Computation
            [thetax,rhox] = cart2pol(IXr,QXr);
            thetax=thetax+2*pi/ac;
            [IXr,QXr]=pol2cart(thetax,rhox);
            [thetay,rhoy] = cart2pol(IYr,QYr);
            thetay=thetay+2*pi/ac;
            [IYr,QYr]=pol2cart(thetay,rhoy);
            %hold on
            %subplot(1,2,1)
            %hold on
            %plot(IXr,QXr,'o');
            %subplot(1,2,2)
            %hold on
            %plot(IYr,QYr,'o')
            close all
            if p==ac
                [mi,pos]=min(ber);
                [thetax,rhox] = cart2pol(IXr,QXr);
                thetax=thetax+2*pi/ac*(pos);
                [IXr,QXr]=pol2cart(thetax,rhox);
                [thetay,rhoy] = cart2pol(IYr,QYr);
                thetay=thetay+2*pi/ac*(pos);
                [IYr,QYr]=pol2cart(thetay,rhoy);
            end
        end
    end
    
    plotwave=0;
    if plotwave==1&&i==8
        xx=1:N;
        plot(xx(1:200),IX(1:200),'-o')
        hold on
        plot(xx(1:2:400),IXrt(1:200),'-gx')
        plot(xx(1:8:1600),IXr(1:200),'-b*')
        plot(xx(1:8:1600),IXb(1:200),'r+')
        
        xlim([0 100])
        title('Example of sending signal');
        xlabel('Sample symbols');ylabel('Amplitude');
        legend('Oversampling of 8', 'Downsample to 2','Downsample to 1','original signal');
    end
    
    
    doPlot=0;
    if i==8
        doPlot=1;%=1 to plot constellation
    end
    if doPlot==1        
        subplot(1,2,1)
        plot(IXr,QXr,'o');
        xlim([-1.5,1.5]);ylim([-1.5,1.5]);
        title('X Polarization')
        subplot(1,2,2)
        plot(IYr,QYr,'o')
        title('Y Polarization')
         xlim([-1.5,1.5]);ylim([-1.5,1.5]);
        
    end
    IXr=(IXr>0)*2-1;%Convert to NRZ
    IYr=(IYr>0)*2-1;
    QXr=(QXr>0)*2-1;
    QYr=(QYr>0)*2-1;
    
    a_cap=[IXr;IYr;QXr;QYr];%Combined signal for BER computation
    BER(i) = sum(a~=a_cap)/length(a);%Bit Error Rate Computation
    if i==8
        for be=1:N
            BERc2(1,be)=sum(a(1:be)~=a_cap(1:be))/be;
        end
    end
end
%BERc=reshape(BERc,1,N*8);

%------Theoretical Bit Error Rate-------------
theoreticalBER = 0.5*erfc(sqrt(EbN0lin));%Theoretical bit error rate
%-------------Plot performance curve------------------------
%figure(2);semilogy(EbN0dB,BER,'k*','LineWidth',1.5); %simulated BER
plo=3;
switch (plo)
    case 1
        plot(1:N,10*log(BERc2),1:N,10*log(BERc3));%,1:N,10*log(BERc2),1:N,10*log(BERc3),'-')
        xlabel('Symbol number');ylabel('BER in dB');
        legend('CMA', 'PU-CMA');
        xlim([0,4096])
        title('Convergent Rate')
    case 2
        snr= 10*log(sum(abs(sx).^2)/sum(abs(nx(:,8).^2)));
        title('Probability of Bit Error for QPSK modulation');
        xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
        legend('Simulated', 'Theoretical');grid on;
    case 3
        figure;semilogy(EbN0dB,BER,'k*','LineWidth',1.5); %simulated BER
        hold on; semilogy(EbN0dB,theoreticalBER,'r-','LineWidth',1.5);
        title('Probability of Bit Error for QPSK modulation');
        xlabel('E_b/N_0 (dB)');ylabel('Probability of Bit Error - P_b');
        legend('Simulated', 'Theoretical');grid on;
end