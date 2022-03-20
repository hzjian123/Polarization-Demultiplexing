% A Matlab cosimulation function that creates an optical signal.
% The function checks if number of samples per bit is a power of two
% and that the vector length (which is the product of TimeWindow and the SampleRate)
% is also a power of two.
% The structures created for the optical ONSL signals MUST be created in the order
% created here, that is first the noise bins, then the channels, then the sampled bands.
% The signal MUST contain all three signal types.
% The contents of the structures for each signal type MUST also be created in the order 
% they are created here.gg

function output = CreateOpticalSignal2()
SampleRate=1024*8;BitRate=1024;Duration=2;

% The sample rate has to give a power of two number of samples per bit
if ~isequal(ceil(log2(SampleRate/BitRate)),floor(log2(SampleRate/BitRate)))
	error(['SampleRate does not yield a number of samples per bit which is a power of two']);
end 
% The sample rate multiplied by the time window must yield a vector which is a power of two
if ~isequal(floor(log2(SampleRate*Duration)),ceil(log2(SampleRate*Duration)))
	error(['SampleRate times TimeWindow must be a power of two']);
end

% Set the type of the signal
y.type = 'osignal';
% Set the boundary conditions of the signal
y.boundaries = 'Periodic';
% Set the time grid spacing
y.dt = 1/SampleRate;
% Set the frequency grid spacing
y.df = 1/Duration;
% Set the time stamp
y.t0 = 0;
% Set the duration of the signal (=Duration/time grid spacing)
y.T = 1/(y.dt * y.df);
% Set the lowest frequency of the structure (grid points)
y.f0 = floor(188.0e12/y.df);
% Set the highest frequency of the structure (grid points)
y.f1 = ceil(198.0e12/y.df);

% Create 10 noise bins with equal spacing describing white noise. 
% The frequency of the first noise bin is 188.0 THz, the bandwidth of
% each noise bin is 1 THz.
y.noise = cell(1,10);
for i = 1 : 10
	% Set the type of the noise bin
	y.noise{i}.type = 'onoise';
	% Set the lower frequency of the noise bin
	y.noise{i}.f0 = floor((188e12 + (i-1) * 1e12) / y.df);
	% Set the upper frequency of the noise bin
	y.noise{i}.f1 = ceil(y.noise{i}.f0 + 1e12/y.df);
	% Initialize the four element stokes vector of the noise bin
	y.noise{i}.S = [1e-10 0 0 0];
end

% Create 5 parameterized channels, the first channel has a carrier
% frequency of 190.0 THz, the channel spacing is 500 GHz.
% The channel power is 1 mW. The bandwidth of the channels is set 
% to 'BitRate'.
y.channels = cell(1,5);
for i = 1 : 5
	% Set the type of the signal
	y.channels{i}.type = 'ochannel';
	% Set the lower frequency of the channel
	y.channels{i}.f0 = floor((190e12 + (i-1) * 500e9 - 0.5*BitRate) /y.df); 
	% Set the upper frequency of the channel
	y.channels{i}.f1 = ceil(y.channels{i}.f0 + BitRate/y.df);
	% Initialize the four element stokes vector of the channel
	y.channels{i}.S = [1e-3 0 0 0];
	% Set the signal statistics attributes to zero
	y.channels{i}.AveragePulsePosition = 0.0;
	y.channels{i}.ExtinctionRatio = 0.0;
	y.channels{i}.PulsePositionVariance = 0.0;
	% Set the default tracking information
	y.channels{i}.Event = [];
	y.channels{i}.Entry = [];
	y.channels{i}.PhysicalPath = [];
	y.channels{i}.TopologicalPosition = [];
	y.channels{i}.Path =[];
	y.channels{i}.AccumulatedGVD = [];
	y.channels{i}.Power = [];
	y.channels{i}.S1 = [];
	y.channels{i}.S2 = [];
	y.channels{i}.S3 = [];
	y.channels{i}.NonlinearCoefficient = [];
	y.channels{i}.AccumulatedDGD = [];
	y.channels{i}.SPM = [];
	y.channels{i}.Bo = [];
	y.channels{i}.Be = [];
	y.channels{i}.NoisePSD = [];
	y.channels{i}.NoisePower = [];
	y.channels{i}.DistCohTotalPower = [];
	y.channels{i}.DistIncohTotalPower = [];
	y.channels{i}.TransitTime = [];
end

% Create 3 sampled bands. 
% The peak power is 1mW, the extinction ratio is ideal, 
% the bit sequence is an alternating (1 0 1 0 1 0 1 ....)
% bit sequence. The channel spacing is 500 GHz.
numberOfBits = BitRate/y.df;
samplesPerBit = SampleRate / BitRate;
y.bands = cell(1,3);
for u = 1 : 3
	% Set the type of the signal
	y.bands{u}.type = 'oband';
	% Set lower frequency of sampled band
	y.bands{u}.f0 = floor((193.0e12 + (u-1) * 500e9 - 0.5*SampleRate)/y.df);
	% Set upper frequency of sampled band
	y.bands{u}.f1 = ceil(y.bands{u}.f0 + SampleRate/y.df);
	% Set polarization to x-polarization
	y.bands{u}.azi = 0;
	y.bands{u}.ell = 0;
	y.bands{u}.E = [];
	% Create an alternating bit stream
	for i = 1 :2: numberOfBits
		for t = 1 : samplesPerBit
			y.bands{u}.E((i-1)*samplesPerBit + t) = complex(sqrt(1e-3),0);
		end
		for t = 1 : samplesPerBit
			y.bands{u}.E(i*samplesPerBit + t) = complex(0,0);
		end
	end 
	% Create empty arrays corresponding to arbitrary polarization
	y.bands{u}.Ex = [];
	y.bands{u}.Ey = [];
end
% Assign the new signal to the output variable
output=y;
