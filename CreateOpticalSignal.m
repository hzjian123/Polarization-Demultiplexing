% A Matlab cosimulation function that creates an optical signal.
% The function checks if number of samples per bit is a power of two
% and that the vector length (which is the product of TimeWindow and the SampleRate)
% is also a power of two.
% The structures created for the optical ONSL signals MUST be created in the order
% created here, that is first the noise bins, then the channels, then the sampled bands.
% The signal MUST contain all three signal types.
% The contents of the structures for each signal type MUST also be created in the order 
% they are created here.
%function output = CreateOpticalSignal()
N=2^11;%Number of bits for each arm
N1=N*8;
SampleRate =40e9;
BitRate=SampleRate/8;
duration=N1/SampleRate;
dt=0;

y.type='osignal';
y.boundaries = 'Periodic';
y.dt = 1/SampleRate;%Period in pico second
y.df=1/duration;
y.t0=int16(0);
y.T=int16(N1);
y.f0 = floor(188.0e12/y.df);
y.f1 = ceil(198.0e12/y.df);% Set the highest frequency of the structure (grid points)

% The sample rate has to give a power of two number of samples per bit
if ~isequal(ceil(log2(SampleRate/BitRate)),floor(log2(SampleRate/BitRate)))
	error(['SampleRate does not yield a number of samples per bit which is a power of two']);
end 
% The sample rate multiplied by the time window must yield a vector which is a power of two
if ~isequal(floor(log2(SampleRate*duration)),ceil(log2(SampleRate*duration)))
	error(['SampleRate times TimeWindow must be a power of two']);
end

    


%y.fc=int32(4e9);
alpha = 0; Nsym =8;L=8;% RC filter alpha and filter span in symbols
rcPulse = raisedCosineFunction(alpha,L,Nsym);%;%Number of bits to transmit
%t=1:length(rcPulse);plot(t,rcPulse,'-o')%Plot RC waveform
EbN0dB = -4:2:10; % Eb/N0 range in dB for simulation
%OF =L; %oversampling factor, sampling frequency will be fs=OF*fc
EbN0lin = 10.^(EbN0dB/10); %converting dB values to linear scale
BER = zeros(length(EbN0dB),1); %For BER values for each Eb/N0
IXb = (rand(N,1)>0.5)*2-1;%X polarization Q arm bit stream initialization
IYb = (rand(N,1)>0.5)*2-1;
QXb = (rand(N,1)>0.5)*2-1;
QYb = (rand(N,1)>0.5)*2-1;
a=[IXb;IYb;QXb;QYb];% Combined vector for BER computation

IXt= zeros(N*L+length(rcPulse),N);%Modulation or carrier matrix
for n = 1:N%for each bit
    IXt(L*(n-1)+1:length(rcPulse)+L*(n-1),n)=rcPulse;%Adding Rc pulse with time shift
end
[IYt,QXt,QYt]=deal(IXt);%Initialization of other 3 arms
IXt=IXt*IXb;IX=IXt(Nsym*L/2+1:end-Nsym*L/2-1);%Multily with bit streams and truncate
IYt=IYt*IYb;IY=IYt(Nsym*L/2+1:end-Nsym*L/2-1);
QXt=QXt*QXb;QX=QXt(Nsym*L/2+1:end-Nsym*L/2-1);
QYt=QYt*QYb;QY=QYt(Nsym*L/2+1:end-Nsym*L/2-1);

sx=IX+QX*1i;%sending signal
sy=IY+QY*1i;

y.noise=cell(0,0);
%y.noise{1}.type = 'onoise';
y.channels=cell(0,0);
%y.channels{1}.type = 'ochannel';
%y.bands=cell(0,0);
y.bands.type = 'oband';
% Set lower frequency of sampled band
y.bands.f0 = floor((193.0e12 + (1-1) * 500e9 - 0.5*SampleRate)/y.df);
% Set upper frequency of sampled band
y.bands.f1 = ceil(y.bands.f0 + SampleRate/y.df);
% Set polarization to x-polarization
y.bands.azi = pi/4;
y.bands.ell = pi/4;
y.bands.E = [];
y.bands.Ex=sx';
y.bands.Ey=sy';




nt = length(IX);% number of points
T=dt*nt;% time window (time duration)
z =100000;                    % propagation distance (in meter)
nz =16;                   % number of steps
dz = z/nz;                  % step size
betapa = [0.000,0.0007];betapb = [-0.000,-0.0007]; %Dispersion polynomial
output=y;