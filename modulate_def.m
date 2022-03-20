function [y] = modulate_def(ss, Option, Sign_mod) 
% 
%modulate  do digital signal modulation/demodulation. 
%  
%   [y] = modulate_def(...) returns modulated/demodulated digital signal. 
% 
%   Input:  
%      ss: received digital signal sequence, vector 
%      Option: modulation type, including DBPSK, DQPSK, M-QAM. String 
%      Sign_mod: 1 -- modulation;  
%               -1 -- demodulation 
% 
%   Output: 
%      y: output signal sequency. vector for modulation, matrix for demoluation. 
 
if nargin < 3 
   	error('MATLAB: modulate_def: NotEnoughInputs',... 
      	'Not enough input arguments. See modulate_def.') 
end 
 
% 
% --------------- Initialization ---------------- 
% 
tmp = strfind(Option, ' '); 
if length(tmp) >= 1 
    Option = Option(1: tmp(1) - 1);  % remove the blank in the last part 
end 
 
if strcmpi(Option, 'DBPSK') 
    N_bit = 1; 
    Map_Mod = [1 -1]; 
     
elseif strcmpi(Option, 'DQPSK') 
    N_bit = 2; 
    % Mapping of dqpsk for modulation 
    Map_Mod = [pi pi/2 pi*3/2 0];       % data: Phase, 00: 0  01: pi/2   11: pi   10: 3*pi/2 
 
elseif strcmpi(Option, 'QPSK') 
    N_bit = 2; 
    % Mapping of QPSK modulation 
    Map_Mod = [1+j]; 
     
elseif strcmpi(Option, '16QAM') 
    N_bit = 4; 
    % Mapping of QAM for modulation 
    Map_16QAM_Modu = [1+1j  1+3j  1-1j  1-3j ... 
                      3+1j  3+3j  3-1j  3-3j ... 
                     -1+1j -1+3j -1-1j -1-3j ... 
                     -3+1j -3+3j -3-1j -3-3j];           %  a + bj : [a0 a1 b0 b1] 
    Map_Mod = Map_16QAM_Modu; 
elseif strcmpi(Option, '32QAM') 
    N_bit = 5; 
    Map_32QAM_Modu = [ -3-5j -1-5j -3+5j -1+5j -5-3j -5-1j -5+3j -5+1j ... 
                       -1-3j -1-1j -1+3j -1+1j -3-3j -3-1j -3+3j -3+1j ... 
                        3-5j  1-5j  3+5j  1+5j  5-3j  5-1j  5+3j  5+1j ... 
                        1-3j  1-1j  1+3j  1+1j  3-3j  3-1j  3+3j  3+1j]; %  a + bj : [a0 a1 a2 b0 b1] 
          % the 32QAM constellation comes from IEEE_2004_"Exact BER Computation for the Cross 32-QAM 
          % Constellation" 
    Map_Mod = Map_32QAM_Modu; 
elseif strcmpi(Option, '64QAM') 
    N_bit = 6; 
    Map_64QAM_Modu = [   3+3j  3+1j  1+3j  1+1j  3+5j  3+7j  1+5j  1+7j ... 
                         5+3j  5+1j  7+3j  7+1j  5+5j  5+7j  7+5j  7+7j ... 
                         3-3j  3-1j  1-3j  1-1j  3-5j  3-7j  1-5j  1-7j ... 
                         5-3j  5-1j  7-3j  7-1j  5-5j  5-7j  7-5j  7-7j ... 
                        -3+3j -3+1j -1+3j -1+1j -3+5j -3+7j -1+5j -1+7j ... 
                        -5+3j -5+1j -7+3j -7+1j -5+5j -5+7j -7+5j -7+7j ... 
                        -3-3j -3-1j -1-3j -1-1j -3-5j -3-7j -1-5j -1-7j ... 
                        -5-3j -5-1j -7-3j -7-1j -5-5j -5-7j -7-5j -7-7j];    %  a + bj : [a0 a1 a2 b0 b1 b2] 
    Map_Mod = Map_64QAM_Modu; 
elseif strcmpi(Option, '128QAM') 
    N_bit = 7; 
    Map_128QAM_Modu = [ -7-7j  -7-5j  -7-1j  -7-3j  -7+7j   -7+5j   -7+1j   -7+3j ...   
                    -5-7j  -5-5j  -5-1j  -5-3j  -5+7j   -5+5j   -5+1j   -5+3j ...   
                    -1-7j  -1-5j  -1-1j  -1-3j  -1+7j   -1+5j   -1+1j   -1+3j  ...  
                    -3-7j  -3-5j  -3-1j  -3-3j  -3+7j   -3+5j   -3+1j   -3+3j  ...  
                     7-7j   7-5j   7-1j   7-3j   7+7j    7+5j    7+1j    7+3j  ...  
                     5-7j   5-5j   5-1j   5-3j   5+7j    5+5j    5+1j    5+3j  ...  
                     1-7j   1-5j   1-1j   1-3j   1+7j    1+5j    1+1j    1+3j  ...  
                     3-7j   3-5j   3-1j   3-3j   3+7j    3+5j    3+1j    3+3j  ...  
                    -9-7j  -9-5j  -9-1j  -9-3j  -9+7j   -9+5j   -9+1j   -9+3j  ...  
                   -11-7j -11-5j -11-1j -11-3j -11+7j  -11+5j  -11+1j  -11+3j  ...  
                    -1-9j -1-11j -7-9j  -7-11j  -1+9j  -1+11j   -7+9j  -7+11j  ...  
                    -3-9j -3-11j -5-9j  -5-11j  -3+9j  -3+11j   -5+9j  -5+11j  ...  
                     9-7j   9-5j  9-1j    9-3j   9+7j    9+5j    9+1j    9+3j  ...  
                    11-7j  11-5j 11-1j   11-3j  11+7j   11+5j   11+1j   11+3j  ...  
                     1-9j  1-11j  7-9j   7-11j   1+9j   1+11j    7+9j   7+11j  ...  
                     3-9j  3-11j  5-9j   5-11j   3+9j   3+11j    5+9j   5+11j];  
                 % the 128QAM constellation comes from 
                 % IEEE_2004_"Implementation Aspects of High-Speed Wireless LAN Systems" 
    Map_Mod = Map_128QAM_Modu; 
elseif strcmpi(Option, '256QAM') 
    N_bit = 8; 
    Map_256QAM_Modu = [ -15+15j -13+15j -9+15j -11+15j -1+15j -3+15j -7+15j -5+15j 15+15j 13+15j 9+15j 11+15j 1+15j 3+15j 7+15j 5+15j ...  
                    -15+13j -13+13j -9+13j -11+13j -1+13j -3+13j -7+13j -5+13j 15+13j 13+13j 9+13j 11+13j 1+13j 3+13j 7+13j 5+13j ...  
                    -15+ 9j -13+ 9j -9+ 9j -11+ 9j -1+ 9j -3+ 9j -7+ 9j -5+ 9j 15+ 9j 13+ 9j 9+ 9j 11+ 9j 1+ 9j 3+ 9j 7+ 9j 5+ 9j ...  
                    -15+11j -13+11j -9+11j -11+11j -1+11j -3+11j -7+11j -5+11j 15+11j 13+11j 9+11j 11+11j 1+11j 3+11j 7+11j 5+11j ...  
                    -15+01j -13+01j -9+01j -11+01j -1+01j -3+01j -7+01j -5+01j 15+01j 13+01j 9+01j 11+01j 1+01j 3+01j 7+01j 5+01j ...  
                    -15+03j -13+03j -9+03j -11+03j -1+03j -3+03j -7+03j -5+03j 15+03j 13+03j 9+03j 11+03j 1+03j 3+03j 7+03j 5+03j ...  
                    -15+07j -13+07j -9+07j -11+07j -1+07j -3+07j -7+07j -5+07j 15+07j 13+07j 9+07j 11+07j 1+07j 3+07j 7+07j 5+07j ...  
                    -15+05j -13+05j -9+05j -11+05j -1+05j -3+05j -7+05j -5+05j 15+05j 13+05j 9+05j 11+05j 1+05j 3+05j 7+05j 5+05j ...  
                    -15+(-15)*j -13+(-15)*j -9+(-15)*j -11+(-15)*j -1+(-15)*j -3+(-15)*j -7+(-15)*j  -5+(-15)*j 15+(-15)*j 13+(-15)*j 9+(-15)*j  11+(-15)*j 1+(-15)*j 3+(-15)*j 7+(-15)*j 5+(-15)*j  ...  
                    -15+(-13)*j -13+(-13)*j -9+(-13)*j -11+(-13)*j -1+(-13)*j -3+(-13)*j -7+(-13)*j  -5+(-13)*j 15+(-13)*j 13+(-13)*j 9+(-13)*j  11+(-13)*j 1+(-13)*j 3+(-13)*j 7+(-13)*j 5+(-13)*j ...  
                     -15+(-9)*j  -13+(-9)*j  -9+(-9)*j  -11+(-9)*j  -1+(-9)*j  -3+(-9)*j  -7+(-9)*j   -5+(-9)*j  15+(-9)*j  13+(-9)*j  9+(-9)*j   11+(-9)*j  1+(-9)*j  3+(-9)*j  7+(-9)*j  5+(-9)*j ...  
                    -15+(-11)*j -13+(-11)*j -9+(-11)*j -11+(-11)*j -1+(-11)*j -3+(-11)*j -7+(-11)*j  -5+(-11)*j 15+(-11)*j 13+(-11)*j 9+(-11)*j  11+(-11)*j 1+(-11)*j 3+(-11)*j 7+(-11)*j 5+(-11)*j ...  
                     -15+(-1)*j  -13+(-1)*j  -9+(-1)*j  -11+(-1)*j  -1+(-1)*j  -3+(-1)*j  -7+(-1)*j   -5+(-1)*j  15+(-1)*j  13+(-1)*j  9+(-1)*j   11+(-1)*j  1+(-1)*j  3+(-1)*j  7+(-1)*j  5+(-1)*j ...  
                     -15+(-3)*j  -13+(-3)*j  -9+(-3)*j  -11+(-3)*j  -1+(-3)*j  -3+(-3)*j  -7+(-3)*j   -5+(-3)*j  15+(-3)*j  13+(-3)*j  9+(-3)*j   11+(-3)*j  1+(-3)*j  3+(-3)*j  7+(-3)*j  5+(-3)*j ...  
                     -15+(-7)*j  -13+(-7)*j  -9+(-7)*j  -11+(-7)*j  -1+(-7)*j  -3+(-7)*j  -7+(-7)*j   -5+(-7)*j  15+(-7)*j  13+(-7)*j  9+(-7)*j   11+(-7)*j  1+(-7)*j  3+(-7)*j  7+(-7)*j  5+(-7)*j ...  
                     -15+(-5)*j  -13+(-5)*j  -9+(-5)*j  -11+(-5)*j  -1+(-5)*j  -3+(-5)*j  -7+(-5)*j   -5+(-5)*j  15+(-5)*j  13+(-5)*j  9+(-5)*j   11+(-5)*j  1+(-5)*j  3+(-5)*j  7+(-5)*j  5+(-5)*j];  
                % IEEE_2003_08__"IEEE P802.15 Wireless Personal Area 
                % Networks" 
 	Map_Mod = Map_256QAM_Modu; 
elseif strcmpi(Option, 'DROP') 
    N_bit = 0; 
else 
    error('Incorrect input parameter, "Option"!') 
end 
P_mean_sqrt = sqrt(mean(abs(Map_Mod).^2));      % Calculate square root of average power 
     
if Sign_mod == 1         
% 
% ------------------------------ Modulation ------------------------------- 
% 
    L_y = length(ss)/N_bit; 
    y = zeros(1, L_y); 
    if N_bit >=4 
        flag = mod(N_bit, 2); 
        tmp_1 = (N_bit - 1: -1: 0)' * ones(1, L_y); 
        tmp_2 = 2.^tmp_1; 
        tmp_3 = reshape(ss, N_bit, []); 
        tmp_4 = sum(tmp_3 .* tmp_2); 
        tmp_4 = tmp_4 + 1;  % make sure that the smallest number is 1. 
       y = (Map_Mod(tmp_4)/P_mean_sqrt);   % Mapping and identify 
%      	tmp_1 = zeros(1, N_O); 
%       	for k_tmp = 1 : N_bit 
%          	tmp_1 = input_mat(Kth_start + k_tmp - 1, : )*2^(N_bit - k_tmp) + tmp_1; 
%        	end 
%         tmp_1 = tmp_1 + 1;                      % make sure that the smallest number is 1. 
%         y = (Map_Mod(tmp_1)/P_mean_sqrt(N_bit));   % Mapping and identify 
         
    elseif strcmpi(Option, 'DQPSK') 
     	% 00: 0  01: pi/2   11: pi   10: 3*pi/2 
        k_tmp = 1; 
    	Phase_tmp = Map_Mod(ss(2*(k_tmp - 1) + 1) + 2*ss(2*(k_tmp - 1) + 2) + 1);    % data location:input(k + 1) input(k) 
      	Phase_last_dqpsk = Phase_tmp; 
        y(k_tmp) = 1*exp(j*Phase_tmp);  
        for k_tmp = 2: L_y 
          	Phase_tmp = mod(Phase_last_dqpsk + Map_Mod(ss(2*(k_tmp - 1) + 1) ... 
                        + 2*ss(2*(k_tmp - 1) + 2) + 1), 2*pi); 
          	Phase_last_dqpsk = Phase_tmp; 
           	y(k_tmp) = 1*exp(j*Phase_tmp);  
        end 
         
    elseif strcmpi(Option, 'QPSK') 
      	temp = (2 * ss(1:2:end) - 1) + j*(2*ss(2:2:end) - 1); 
     	y = temp/P_mean_sqrt;  
         
    elseif strcmpi(Option, 'DBPSK') 
     	tmp_1 = xor(ss(1), 1);           	% iniciate the original phase as 0+j0 
       	tmp_2 = tmp_1; 
     	y(1) = 2*tmp_1 - 1; 
      	for k_tmp = 2 : length(ss) 
          	tmp_1 = xor(ss(k_tmp), tmp_2);      % 1: phase:0 (1)   0: -pi (-1) 
          	tmp_2 = tmp_1; 
         	y(k_tmp) = 2*tmp_1 - 1; 
        end 
 
    else 
        error('Incorrect input parameter, "Option"!') 
    end 
 
 
elseif Sign_mod == -1   % Demodulation 
% 
% ------------------------------ Demodulation ----------------------------- 
%  
    L_y = length(ss); 
    y = zeros(N_bit, L_y); 
    if N_bit >= 4 
        % Initialization 
        if strcmpi(Option, '16QAM')  
            % Mapping of QAM for demodulation 
            Map_16QAM_Demodu = [1 1 1 1; 1 1 1 0; 1 1 0 0; 1 1 0 1; 
                    1 0 1 1; 1 0 1 0; 1 0 0 0; 1 0 0 1; 
                    0 0 1 1; 0 0 1 0; 0 0 0 0; 0 0 0 1;  
                    0 1 1 1; 0 1 1 0; 0 1 0 0; 0 1 0 1]; 
            Map_Demod = Map_16QAM_Demodu; 
         
        elseif strcmpi(Option, '32QAM') 
            Map_32QAM_Demodu = [1 1 1 1 1 ; 0 0 1 0 0 ; 0 0 1 0 1 ; 0 0 1 1 1 ; 0 0 1 1 0 ; 1 1 1 1 1;   1 1 1 1 1;   1 1 1 1 1;  % the last two column group is usefulness   
                    0 0 0 0 0 ; 0 1 1 0 0 ; 0 1 1 0 1 ; 0 1 1 1 1 ; 0 1 1 1 0 ; 0 0 0 1 0;   1 1 1 1 1;   1 1 1 1 1;   
                    0 0 0 0 1 ; 0 1 0 0 0 ; 0 1 0 0 1 ; 0 1 0 1 1 ; 0 1 0 1 0 ; 0 0 0 1 1;   1 1 1 1 1;   1 1 1 1 1;   
                    1 0 0 0 1 ; 1 1 0 0 0 ; 1 1 0 0 1 ; 1 1 0 1 1 ; 1 1 0 1 0 ; 1 0 0 1 1;   1 1 1 1 1;   1 1 1 1 1;   
                    1 0 0 0 0 ; 1 1 1 0 0 ; 1 1 1 0 1 ; 1 1 1 1 1 ; 1 1 1 1 0 ; 1 0 0 1 0;   1 1 1 1 1;   1 1 1 1 1;   
                    1 1 1 1 1 ; 1 0 1 0 0 ; 1 0 1 0 1 ; 1 0 1 1 1 ; 1 0 1 1 0 ; 1 1 1 1 1];  
            Map_Demod = Map_32QAM_Demodu; 
 
        elseif strcmpi(Option, '64QAM') 
            Map_64QAM_Demodu = [1 1 1 1 1 1; 1 1 1 1 1 0; 1 1 1 0 1 0; 1 1 1 0 1 1; 1 0 1 0 1 1; 1 0 1 0 1 0; 1 0 1 1 1 0; 1 0 1 1 1 1; ... 
                    1 1 1 1 0 1; 1 1 1 1 0 0; 1 1 1 0 0 0; 1 1 1 0 0 1; 1 0 1 0 0 1; 1 0 1 0 0 0; 1 0 1 1 0 0; 1 0 1 1 0 1; ... 
                    1 1 0 1 0 1; 1 1 0 1 0 0; 1 1 0 0 0 0; 1 1 0 0 0 1; 1 0 0 0 0 1; 1 0 0 0 0 0; 1 0 0 1 0 0; 1 0 0 1 0 1; ... 
                    1 1 0 1 1 1; 1 1 0 1 1 0; 1 1 0 0 1 0; 1 1 0 0 1 1; 1 0 0 0 1 1; 1 0 0 0 1 0; 1 0 0 1 1 0; 1 0 0 1 1 1; ... 
                    0 1 0 1 1 1; 0 1 0 1 1 0; 0 1 0 0 1 0; 0 1 0 0 1 1; 0 0 0 0 1 1; 0 0 0 0 1 0; 0 0 0 1 1 0; 0 0 0 1 1 1; ... 
                    0 1 0 1 0 1; 0 1 0 1 0 0; 0 1 0 0 0 0; 0 1 0 0 0 1; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 1 0 1; ... 
                    0 1 1 1 0 1; 0 1 1 1 0 0; 0 1 1 0 0 0; 0 1 1 0 0 1; 0 0 1 0 0 1; 0 0 1 0 0 0; 0 0 1 1 0 0; 0 0 1 1 0 1; ... 
                    0 1 1 1 1 1; 0 1 1 1 1 0; 0 1 1 0 1 0; 0 1 1 0 1 1; 0 0 1 0 1 1; 0 0 1 0 1 0; 0 0 1 1 1 0; 0 0 1 1 1 1]; 
            Map_Demod = Map_64QAM_Demodu; 
 
        elseif strcmpi(Option, '128QAM') 
                Map_128QAM_Demodu = [1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 0 0 1 0 0 0 ; 1 0 0 1 0 0 1 ; 1 0 0 1 0 1 1 ; 1 0 0 1 0 1 0 ; 1 0 0 1 1 1 0 ; 1 0 0 1 1 1 1; 
                     1 0 0 1 1 0 1 ; 1 0 0 1 1 0 0 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1;  
                      
                     1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 0 0 0 0 0 0 ; 1 0 0 0 0 0 1 ; 1 0 0 0 0 1 1 ; 1 0 0 0 0 1 0 ; 1 0 0 0 1 1 0 ; 1 0 0 0 1 1 1;  
                     1 0 0 0 1 0 1 ; 1 0 0 0 1 0 0 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1;  
                      
                     1 0 1 0 0 1 1 ; 1 0 1 0 0 1 0 ; 0 0 0 0 0 0 0 ; 0 0 0 0 0 0 1 ; 0 0 0 0 0 1 1 ; 0 0 0 0 0 1 0 ; 0 0 0 0 1 1 0 ; 0 0 0 0 1 1 1; 
                     0 0 0 0 1 0 1 ; 0 0 0 0 1 0 0 ; 1 0 1 0 1 1 0 ; 1 0 1 0 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 0 1 1 0 1 1 ; 1 0 1 1 0 1 0 ; 0 0 0 1 0 0 0 ; 0 0 0 1 0 0 1 ; 0 0 0 1 0 1 1 ; 0 0 0 1 0 1 0 ; 0 0 0 1 1 1 0 ; 0 0 0 1 1 1 1;  
                     0 0 0 1 1 0 1 ; 0 0 0 1 1 0 0 ; 1 0 1 1 1 1 0 ; 1 0 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 0 1 1 0 0 1 ; 1 0 1 1 0 0 0 ; 0 0 1 1 0 0 0 ; 0 0 1 1 0 0 1 ; 0 0 1 1 0 1 1 ; 0 0 1 1 0 1 0 ; 0 0 1 1 1 1 0 ; 0 0 1 1 1 1 1; 
                     0 0 1 1 1 0 1 ; 0 0 1 1 1 0 0 ; 1 0 1 1 1 0 0 ; 1 0 1 1 1 0 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 0 1 0 0 0 1 ; 1 0 1 0 0 0 0 ; 0 0 1 0 0 0 0 ; 0 0 1 0 0 0 1 ; 0 0 1 0 0 1 1 ; 0 0 1 0 0 1 0 ; 0 0 1 0 1 1 0 ; 0 0 1 0 1 1 1;  
                     0 0 1 0 1 0 1 ; 0 0 1 0 1 0 0 ; 1 0 1 0 1 0 0 ; 1 0 1 0 1 0 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 0 0 0 1 ; 1 1 1 0 0 0 0 ; 0 1 1 0 0 0 0 ; 0 1 1 0 0 0 1 ; 0 1 1 0 0 1 1 ; 0 1 1 0 0 1 0 ; 0 1 1 0 1 1 0 ; 0 1 1 0 1 1 1; 
                     0 1 1 0 1 0 1 ; 0 1 1 0 1 0 0 ; 1 1 1 0 1 0 0 ; 1 1 1 0 1 0 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 1 0 0 1 ; 1 1 1 1 0 0 0 ; 0 1 1 1 0 0 0 ; 0 1 1 1 0 0 1 ; 0 1 1 1 0 1 1 ; 0 1 1 1 0 1 0 ; 0 1 1 1 1 1 0 ; 0 1 1 1 1 1 1;  
                     0 1 1 1 1 0 1 ; 0 1 1 1 1 0 0 ; 1 1 1 1 1 0 0 ; 1 1 1 1 1 0 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 1 0 1 1 ; 1 1 1 1 0 1 0 ; 0 1 0 1 0 0 0 ; 0 1 0 1 0 0 1 ; 0 1 0 1 0 1 1 ; 0 1 0 1 0 1 0 ; 0 1 0 1 1 1 0 ; 0 1 0 1 1 1 1;  
                     0 1 0 1 1 0 1 ; 0 1 0 1 1 0 0 ; 1 1 1 1 1 1 0 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 0 0 1 1 ; 1 1 1 0 0 1 0 ; 0 1 0 0 0 0 0 ; 0 1 0 0 0 0 1 ; 0 1 0 0 0 1 1 ; 0 1 0 0 0 1 0 ; 0 1 0 0 1 1 0 ; 0 1 0 0 1 1 1;  
                     0 1 0 0 1 0 1 ; 0 1 0 0 1 0 0 ; 1 1 1 0 1 1 0 ; 1 1 1 0 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 0 0 0 0 0 ; 1 1 0 0 0 0 1 ; 1 1 0 0 0 1 1 ; 1 1 0 0 0 1 0 ; 1 1 0 0 1 1 0 ; 1 1 0 0 1 1 1;  
                     1 1 0 0 1 0 1 ; 1 1 0 0 1 0 0 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1; 
                      
                     1 1 1 1 1 1 1 ; 1 1 1 1 1 1 1 ; 1 1 0 1 0 0 0 ; 1 1 0 1 0 0 1 ; 1 1 0 1 0 1 1 ; 1 1 0 1 0 1 0 ; 1 1 0 1 1 1 0 ; 1 1 0 1 1 1 1;  
                     1 1 0 1 1 0 1 ; 1 1 0 1 1 0 0]; 
            Map_Demod = Map_128QAM_Demodu; 
                  
        elseif strcmpi(Option, '256QAM') 
            Map_256QAM_Demodu = [ 1 0 0 0 0 0 0 0; 1 0 0 1 0 0 0 0; 1 0 1 1 0 0 0 0; 1 0 1 0 0 0 0 0; 1 1 1 0 0 0 0 0; 1 1 1 1 0 0 0 0; 1 1 0 1 0 0 0 0; 1 1 0 0 0 0 0 0; 0 1 0 0 0 0 0 0; 0 1 0 1 0 0 0 0; 0 1 1 1 0 0 0 0; 0 1 1 0 0 0 0 0; 0 0 1 0 0 0 0 0; 0 0 1 1 0 0 0 0; 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0;  
                      1 0 0 0 0 0 0 1; 1 0 0 1 0 0 0 1; 1 0 1 1 0 0 0 1; 1 0 1 0 0 0 0 1; 1 1 1 0 0 0 0 1; 1 1 1 1 0 0 0 1; 1 1 0 1 0 0 0 1; 1 1 0 0 0 0 0 1; 0 1 0 0 0 0 0 1; 0 1 0 1 0 0 0 1; 0 1 1 1 0 0 0 1; 0 1 1 0 0 0 0 1; 0 0 1 0 0 0 0 1; 0 0 1 1 0 0 0 1; 0 0 0 1 0 0 0 1; 0 0 0 0 0 0 0 1;  
                      1 0 0 0 0 0 1 1; 1 0 0 1 0 0 1 1; 1 0 1 1 0 0 1 1; 1 0 1 0 0 0 1 1; 1 1 1 0 0 0 1 1; 1 1 1 1 0 0 1 1; 1 1 0 1 0 0 1 1; 1 1 0 0 0 0 1 1; 0 1 0 0 0 0 1 1; 0 1 0 1 0 0 1 1; 0 1 1 1 0 0 1 1; 0 1 1 0 0 0 1 1; 0 0 1 0 0 0 1 1; 0 0 1 1 0 0 1 1; 0 0 0 1 0 0 1 1; 0 0 0 0 0 0 1 1;  
                      1 0 0 0 0 0 1 0; 1 0 0 1 0 0 1 0; 1 0 1 1 0 0 1 0; 1 0 1 0 0 0 1 0; 1 1 1 0 0 0 1 0; 1 1 1 1 0 0 1 0; 1 1 0 1 0 0 1 0; 1 1 0 0 0 0 1 0; 0 1 0 0 0 0 1 0; 0 1 0 1 0 0 1 0; 0 1 1 1 0 0 1 0; 0 1 1 0 0 0 1 0; 0 0 1 0 0 0 1 0; 0 0 1 1 0 0 1 0; 0 0 0 1 0 0 1 0; 0 0 0 0 0 0 1 0;  
                      1 0 0 0 0 1 1 0; 1 0 0 1 0 1 1 0; 1 0 1 1 0 1 1 0; 1 0 1 0 0 1 1 0; 1 1 1 0 0 1 1 0; 1 1 1 1 0 1 1 0; 1 1 0 1 0 1 1 0; 1 1 0 0 0 1 1 0; 0 1 0 0 0 1 1 0; 0 1 0 1 0 1 1 0; 0 1 1 1 0 1 1 0; 0 1 1 0 0 1 1 0; 0 0 1 0 0 1 1 0; 0 0 1 1 0 1 1 0; 0 0 0 1 0 1 1 0; 0 0 0 0 0 1 1 0;  
                      1 0 0 0 0 1 1 1; 1 0 0 1 0 1 1 1; 1 0 1 1 0 1 1 1; 1 0 1 0 0 1 1 1; 1 1 1 0 0 1 1 1; 1 1 1 1 0 1 1 1; 1 1 0 1 0 1 1 1; 1 1 0 0 0 1 1 1; 0 1 0 0 0 1 1 1; 0 1 0 1 0 1 1 1; 0 1 1 1 0 1 1 1; 0 1 1 0 0 1 1 1; 0 0 1 0 0 1 1 1; 0 0 1 1 0 1 1 1; 0 0 0 1 0 1 1 1; 0 0 0 0 0 1 1 1;  
                      1 0 0 0 0 1 0 1; 1 0 0 1 0 1 0 1; 1 0 1 1 0 1 0 1; 1 0 1 0 0 1 0 1; 1 1 1 0 0 1 0 1; 1 1 1 1 0 1 0 1; 1 1 0 1 0 1 0 1; 1 1 0 0 0 1 0 1; 0 1 0 0 0 1 0 1; 0 1 0 1 0 1 0 1; 0 1 1 1 0 1 0 1; 0 1 1 0 0 1 0 1; 0 0 1 0 0 1 0 1; 0 0 1 1 0 1 0 1; 0 0 0 1 0 1 0 1; 0 0 0 0 0 1 0 1;  
                      1 0 0 0 0 1 0 0; 1 0 0 1 0 1 0 0; 1 0 1 1 0 1 0 0; 1 0 1 0 0 1 0 0; 1 1 1 0 0 1 0 0; 1 1 1 1 0 1 0 0; 1 1 0 1 0 1 0 0; 1 1 0 0 0 1 0 0; 0 1 0 0 0 1 0 0; 0 1 0 1 0 1 0 0; 0 1 1 1 0 1 0 0; 0 1 1 0 0 1 0 0; 0 0 1 0 0 1 0 0; 0 0 1 1 0 1 0 0; 0 0 0 1 0 1 0 0; 0 0 0 0 0 1 0 0;  
                      1 0 0 0 1 1 0 0; 1 0 0 1 1 1 0 0; 1 0 1 1 1 1 0 0; 1 0 1 0 1 1 0 0; 1 1 1 0 1 1 0 0; 1 1 1 1 1 1 0 0; 1 1 0 1 1 1 0 0; 1 1 0 0 1 1 0 0; 0 1 0 0 1 1 0 0; 0 1 0 1 1 1 0 0; 0 1 1 1 1 1 0 0; 0 1 1 0 1 1 0 0; 0 0 1 0 1 1 0 0; 0 0 1 1 1 1 0 0; 0 0 0 1 1 1 0 0; 0 0 0 0 1 1 0 0;  
                      1 0 0 0 1 1 0 1; 1 0 0 1 1 1 0 1; 1 0 1 1 1 1 0 1; 1 0 1 0 1 1 0 1; 1 1 1 0 1 1 0 1; 1 1 1 1 1 1 0 1; 1 1 0 1 1 1 0 1; 1 1 0 0 1 1 0 1; 0 1 0 0 1 1 0 1; 0 1 0 1 1 1 0 1; 0 1 1 1 1 1 0 1; 0 1 1 0 1 1 0 1; 0 0 1 0 1 1 0 1; 0 0 1 1 1 1 0 1; 0 0 0 1 1 1 0 1; 0 0 0 0 1 1 0 1;  
                      1 0 0 0 1 1 1 1; 1 0 0 1 1 1 1 1; 1 0 1 1 1 1 1 1; 1 0 1 0 1 1 1 1; 1 1 1 0 1 1 1 1; 1 1 1 1 1 1 1 1; 1 1 0 1 1 1 1 1; 1 1 0 0 1 1 1 1; 0 1 0 0 1 1 1 1; 0 1 0 1 1 1 1 1; 0 1 1 1 1 1 1 1; 0 1 1 0 1 1 1 1; 0 0 1 0 1 1 1 1; 0 0 1 1 1 1 1 1; 0 0 0 1 1 1 1 1; 0 0 0 0 1 1 1 1;  
                      1 0 0 0 1 1 1 0; 1 0 0 1 1 1 1 0; 1 0 1 1 1 1 1 0; 1 0 1 0 1 1 1 0; 1 1 1 0 1 1 1 0; 1 1 1 1 1 1 1 0; 1 1 0 1 1 1 1 0; 1 1 0 0 1 1 1 0; 0 1 0 0 1 1 1 0; 0 1 0 1 1 1 1 0; 0 1 1 1 1 1 1 0; 0 1 1 0 1 1 1 0; 0 0 1 0 1 1 1 0; 0 0 1 1 1 1 1 0; 0 0 0 1 1 1 1 0; 0 0 0 0 1 1 1 0;  
                      1 0 0 0 1 0 1 0; 1 0 0 1 1 0 1 0; 1 0 1 1 1 0 1 0; 1 0 1 0 1 0 1 0; 1 1 1 0 1 0 1 0; 1 1 1 1 1 0 1 0; 1 1 0 1 1 0 1 0; 1 1 0 0 1 0 1 0; 0 1 0 0 1 0 1 0; 0 1 0 1 1 0 1 0; 0 1 1 1 1 0 1 0; 0 1 1 0 1 0 1 0; 0 0 1 0 1 0 1 0; 0 0 1 1 1 0 1 0; 0 0 0 1 1 0 1 0; 0 0 0 0 1 0 1 0;  
                      1 0 0 0 1 0 1 1; 1 0 0 1 1 0 1 1; 1 0 1 1 1 0 1 1; 1 0 1 0 1 0 1 1; 1 1 1 0 1 0 1 1; 1 1 1 1 1 0 1 1; 1 1 0 1 1 0 1 1; 1 1 0 0 1 0 1 1; 0 1 0 0 1 0 1 1; 0 1 0 1 1 0 1 1; 0 1 1 1 1 0 1 1; 0 1 1 0 1 0 1 1; 0 0 1 0 1 0 1 1; 0 0 1 1 1 0 1 1; 0 0 0 1 1 0 1 1; 0 0 0 0 1 0 1 1;  
                      1 0 0 0 1 0 0 1; 1 0 0 1 1 0 0 1; 1 0 1 1 1 0 0 1; 1 0 1 0 1 0 0 1; 1 1 1 0 1 0 0 1; 1 1 1 1 1 0 0 1; 1 1 0 1 1 0 0 1; 1 1 0 0 1 0 0 1; 0 1 0 0 1 0 0 1; 0 1 0 1 1 0 0 1; 0 1 1 1 1 0 0 1; 0 1 1 0 1 0 0 1; 0 0 1 0 1 0 0 1; 0 0 1 1 1 0 0 1; 0 0 0 1 1 0 0 1; 0 0 0 0 1 0 0 1;  
                      1 0 0 0 1 0 0 0; 1 0 0 1 1 0 0 0; 1 0 1 1 1 0 0 0; 1 0 1 0 1 0 0 0; 1 1 1 0 1 0 0 0; 1 1 1 1 1 0 0 0; 1 1 0 1 1 0 0 0; 1 1 0 0 1 0 0 0; 0 1 0 0 1 0 0 0; 0 1 0 1 1 0 0 0; 0 1 1 1 1 0 0 0; 0 1 1 0 1 0 0 0; 0 0 1 0 1 0 0 0; 0 0 1 1 1 0 0 0; 0 0 0 1 1 0 0 0; 0 0 0 0 1 0 0 0]; 
            Map_Demod = Map_256QAM_Demodu; 
        end 
 
        flag_odd = mod(N_bit,2); 
        sequ_Modu = ss * P_mean_sqrt;           % recover the data identified 
                Val_half = 2^((N_bit - flag_odd)/2)/2;     % half of half an length of 'Zheng fang xing' 
        if flag_odd == 1  	% when the number of bits modulated is odd(ex. 32QAM, 128QAM) 
            if N_bit == 5               % 32 QAM 
                tmp_1 = sequ_Modu; 
                tmp_2 = (abs(real(sequ_Modu)) > (2*Val_half)) & (abs(imag(sequ_Modu)) > (2*Val_half)); 
                tmp_3 = tmp_1.* tmp_2;   % the border points 
                tmp_4 = abs(imag(tmp_3)) >= abs(real(tmp_3)); 
                tmp_5 = (3 + 5j)*(( imag(tmp_3) > 0 & real(tmp_3) > 0 ) & tmp_4) + ... 
                        (5 + 3j)*(( imag(tmp_3) > 0 & real(tmp_3) > 0 ) & not(tmp_4)) + ... 
                       (-3 + 5j)*(( imag(tmp_3) > 0 & real(tmp_3) < 0 ) & tmp_4) + ... 
                       (-5 + 3j)*(( imag(tmp_3) > 0 & real(tmp_3) < 0 ) & not(tmp_4)) + ... 
                       (-3 - 5j)*(( imag(tmp_3) < 0 & real(tmp_3) < 0 ) & tmp_4) + ... 
                       (-5 - 3j)*(( imag(tmp_3) < 0 & real(tmp_3) < 0 ) & not(tmp_4)) + ... 
                        (3 - 5j)*(( imag(tmp_3) < 0 & real(tmp_3) > 0 ) & tmp_4) + ... 
                        (5 - 3j)*(( imag(tmp_3) < 0 & real(tmp_3) > 0 ) & not(tmp_4)); 
                sequ_Modu = tmp_5 + tmp_1.*(not(tmp_2));       
            elseif N_bit == 7           % 128 QAM 
              	tmp_1 = sequ_Modu; 
                tmp_2 = (abs(real(sequ_Modu)) > (2*Val_half)) & (abs(imag(sequ_Modu)) > (2*Val_half));     % outside 'Zheng fang xing' 
                tmp_3 = (abs(real(sequ_Modu)) > (2*Val_half + 2)) & (abs(imag(sequ_Modu)) > (2*Val_half + 2)); % outside 'zheng fang xing + 2' 
                tmp_31 = tmp_2 & ( (abs(real(sequ_Modu)) < (2*Val_half + 2)) & (abs(imag(sequ_Modu)) < (2*Val_half + 2))); 
                tmp_4 = tmp_1.* tmp_31;                 % the inner border points 
                tmp_5 = tmp_1.* tmp_2.*(not(tmp_31));   % the externer border points 
                tmp_6 = abs(imag(tmp_1)) >= abs(real(tmp_1)); 
                tmp_7 = (7 + 9j) *(( imag(tmp_4) > 0 & real(tmp_4) > 0 ) & tmp_6) + ... 
                        (9 + 7j) *(( imag(tmp_4) > 0 & real(tmp_4) > 0 ) & not(tmp_6)) + ... 
                        (7 + 11j)*(( imag(tmp_5) > 0 & real(tmp_5) > 0 ) & tmp_6) + ... 
                        (11 + 7j)*(( imag(tmp_5) > 0 & real(tmp_5) > 0 ) & not(tmp_6)) + ... 
                       (-7 + 9j) *(( imag(tmp_4) > 0 & real(tmp_4) < 0 ) & tmp_6) + ... 
                       (-9 + 7j) *(( imag(tmp_4) > 0 & real(tmp_4) < 0 ) & not(tmp_6)) + ... 
                       (-7 + 11j)*(( imag(tmp_5) > 0 & real(tmp_5) < 0 ) & tmp_6) + ... 
                       (-11 + 7j)*(( imag(tmp_5) > 0 & real(tmp_5) < 0 ) & not(tmp_6)) + ... 
                        (-7 - 9j)*(( imag(tmp_4) < 0 & real(tmp_4) < 0 ) & tmp_6) + ... 
                       (-9 - 7j)*(( imag(tmp_4) < 0 & real(tmp_4) < 0 ) & not(tmp_6)) + ... 
                       (-7 - 11j)*(( imag(tmp_5) < 0 & real(tmp_5) < 0 ) & tmp_6) + ... 
                       (-11 - 7j)*(( imag(tmp_5) < 0 & real(tmp_5) < 0 ) & not(tmp_6)) + ... 
                        (7 - 9j)*(( imag(tmp_4) < 0 & real(tmp_4) > 0 ) & tmp_6) + ... 
                        (9 - 7j)*(( imag(tmp_4) < 0 & real(tmp_4) > 0 ) & not(tmp_6)) + ... 
                        (7 - 11j)*(( imag(tmp_5) < 0 & real(tmp_5) > 0 ) & tmp_6) + ... 
                        (11 - 7j)*(( imag(tmp_5) < 0 & real(tmp_5) > 0 ) & not(tmp_6)); 
                sequ_Modu = tmp_7 + tmp_1.*(not(tmp_2));  
            else 
                disp('There is something wrong!!!') 
                return 
            end 
            flag_tmp = (floor(N_bit/2) - 1);          % 32QAM: 1  128QAM: 2  --02/05/07 
            Val_ss_Re = round( (real(sequ_Modu) + 2*Val_half + 2*flag_tmp - 1)/2 );  % temp_r 
            Val_ss_Im = round( (imag(sequ_Modu) + 2*Val_half + 2*flag_tmp - 1)/2 );  % temp_i 
 
            sequ_Modu_Re = zeros(1, L_y).*(Val_ss_Re<0) + Val_ss_Re .* ((Val_ss_Re >= 0) & (Val_ss_Re <= 2^((N_bit - flag_odd)/2) + 2*flag_tmp - 1) ) ... 
                    + (2^((N_bit - flag_odd)/2) + 2*flag_tmp - 1)*(Val_ss_Re > 2^((N_bit + flag_odd)/2) + 2*flag_tmp - 1); 
            sequ_Modu_Im = zeros(1, L_y).*(Val_ss_Im<0) + Val_ss_Im .* ((Val_ss_Im >= 0) & (Val_ss_Im <= 2^((N_bit - flag_odd)/2) + 2*flag_tmp - 1) ) ... 
                    + (2^((N_bit - flag_odd)/2) + 2*flag_tmp - 1)*(Val_ss_Im > 2^((N_bit - flag_odd)/2) + 2*flag_tmp - 1); 
 
            tmp_1 = zeros(N_bit + flag_odd, L_y); 
            for k_tmp = 1 : (N_bit + flag_odd)/2 
                tmp_1((N_bit + flag_odd)/2 - k_tmp + 1, : ) = mod(floor(sequ_Modu_Re/(2^(k_tmp - 1))),2);                 	% real part 
                tmp_1((N_bit + flag_odd)/2 - k_tmp + 1 + (N_bit + flag_odd)/2, : ) = mod(floor(sequ_Modu_Im/(2^(k_tmp - 1))),2);   	% imaginary part 
            end 
            tmp_2 = zeros(1,L_y); 
            for k_tmp = 1 : (N_bit + flag_odd) 
                tmp_2 = tmp_2 + tmp_1(k_tmp, : )* 2^((N_bit + flag_odd) - k_tmp); 
            end 
            tmp_2 = tmp_2 + 1;           % make sure the smallest number is 1 
            y = (Map_Demod(tmp_2, :))'; 
            
        else 	% when the number of bits modulated is even(ex. 16QAM, 64QAM, 256) 
            Val_ss_Re = round( (real(sequ_Modu) + 2*Val_half - 1)/2 );  % temp_r 
            Val_ss_Im = round( (imag(sequ_Modu) + 2*Val_half - 1)/2 );  % temp_i 
            sequ_Modu_Re = zeros(1, L_y).*(Val_ss_Re<0) + Val_ss_Re .* ((Val_ss_Re >= 0) & (Val_ss_Re <= 2^((N_bit)/2) - 1) ) ... 
                    + (2^((N_bit)/2) - 1)*(Val_ss_Re > 2^((N_bit)/2) - 1); 
            sequ_Modu_Im = zeros(1, L_y).*(Val_ss_Im<0) + Val_ss_Im .* ((Val_ss_Im >= 0) & (Val_ss_Im <= 2^((N_bit)/2) - 1) ) ... 
                    + (2^((N_bit)/2) - 1)*(Val_ss_Im > 2^((N_bit)/2) - 1); 
            tmp_1 = zeros(N_bit, L_y); 
            for k_tmp = 1 : N_bit/2 
                if N_bit == 8 
                    tmp_1(N_bit/2 - k_tmp + 1, : ) = mod(floor(sequ_Modu_Re/(2^(k_tmp - 1))),2);                 	% real part 
                    tmp_1(N_bit/2 - k_tmp + 1 + N_bit/2, : ) = mod(floor(sequ_Modu_Im/(2^(k_tmp - 1))),2);   	% imaginary part 
                else 
                tmp_1(N_bit/2 - k_tmp + 1, : ) = mod(floor(sequ_Modu_Re/(2^(k_tmp - 1))),2);                 	% real part 
                tmp_1(N_bit/2 - k_tmp + 1 + N_bit/2, : ) = mod(floor(sequ_Modu_Im/(2^(k_tmp - 1))),2);   	% imaginary part 
                end 
            end 
            tmp_2 = zeros(1,L_y); 
            for k_tmp = 1 : N_bit 
                tmp_2 = tmp_2 + tmp_1(k_tmp, : )* 2^(N_bit - k_tmp); 
            end 
            tmp_2 = tmp_2 + 1;           % make sure the smallest number is 1 
            y = (Map_Demod(tmp_2, :))'; 
        end 
         
    elseif strcmpi(Option, 'DQPSK') 
        sequ_Modu = ss; 
        k_tmp = 1; 
      	Phase_tmp = angle(sequ_Modu(k_tmp));         
       	Phase_last_dqpsk = angle(sequ_Modu(k_tmp)); 
    	temp = 1*exp(j*(Phase_tmp + pi/4)); 
      	y(1, k_tmp) = imag(temp) > 0;       % data location:output(k + 1) output(k) 
      	y(2, k_tmp) = real(temp) > 0; 
        for k_tmp = 2 : length(ss) 
         	Phase_tmp = mod((angle(sequ_Modu(k_tmp)) - Phase_last_dqpsk + 2*pi), 2*pi); 
        	Phase_last_dqpsk = angle(sequ_Modu(k_tmp)); 
          	temp = 1*exp(j*(Phase_tmp + pi/4)); 
         	y(1, k_tmp) = imag(temp) > 0;       % data location:output(k + 1) output(k) 
           	y(2, k_tmp) = real(temp) > 0; 
        end 
 
    elseif strcmpi(Option, 'QPSK')   
        sequ_Modu = ss; 
      	y(1, :) = real(sequ_Modu * P_mean_sqrt) > 0; 
      	y(2, :) = imag(sequ_Modu * P_mean_sqrt) > 0; 
             
    elseif strcmpi(Option, 'DBPSK')             
      	sequ_Modu = ss > 0; 
        y = [xor(sequ_Modu(1), 1) xor(sequ_Modu(2:end),sequ_Modu(1:end - 1))]; 
 
    end 
     
else 
    error('Incorrect input parameter, "Sign_mod"!') 
end 
% end of the file
