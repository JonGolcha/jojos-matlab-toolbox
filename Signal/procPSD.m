function [PSD, F_AXIS] = procPSD(DATA, WINDOW_TYPE, WINDOW_SIZE, OVERLAP, FS, FREQUENCY_RANGE, MASK)
%  [PSD, F_AXIS] = procPSD(DATA, WINDOW_TYPE, WINDOW_SIZE, OVERLAP, FS, FREQUENCY_RANGE, MASK)
%
% This function computes power spectrum using Welch's method. 
%
% *** Inputs *** 
% - DATA            --> matrix with the data (M samples by N channels)
% - WINDOW_TYPE     --> (optionnal) string with the type of the window
% - WINDOW_SIZE     --> (optionnal) size of the window in number of samples 
% - OVERLAP         --> (optionnal) percentage of the overlapping window (from 0 to 0.99)
% - FS              --> (optionnal) specifies frequency rate for scaling F_AXIS
% - FREQUENCY_RANGE --> (optionnal) range of interesting frequencies to keep, 
%                                   format must be [sample_rate min_frec max_frec]
% - MASK            --> (optionnal) binary mask (1 x DATA_LENGTH) 
%                                   -> epochs with at least one sample set to 0 are not considered
%  *** Outputs *** 
% - S               --> matrix of power spectra (f by N)
% - F_AXIS          --> frequency axis for display (1 by f)
%
%  *** History ***
% First version: Jonas Chatel-Goldman @OIST, 07/29/2014
% jonas.chatel.goldman(at)gmail.com


%% Assess inputs
% default mask
if(nargin < 7)||(isempty(MASK))                 MASK  = ones(1,size(DATA,1));  end
% default frequency range
if(nargin < 6)||(isempty(FREQUENCY_RANGE))      FREQUENCY_RANGE = []; end
% default frequency rate
if(nargin < 5)||(isempty(FS))                   FS = 1;  end
% default overlapping
if(nargin < 4)||(isempty(OVERLAP))              OVERLAP = .75;  end
% default windows size
if(nargin < 3)||(isempty(WINDOW_SIZE))          WINDOW_SIZE = 128;  end
% default windows type
if(nargin < 2)||(isempty(WINDOW_TYPE))          WINDOW_TYPE = 'hanning';  end

% Adjust the window to the next pow of 2 (for the FFT)
WINDOW_SIZE = 2 ^ nextpow2(WINDOW_SIZE);
% Calculate the number of frequencies
number_freqs = (WINDOW_SIZE / 2); % +1;
s_data = size(DATA);


%% Create the window from the inputs and preallocates memory
switch lower(WINDOW_TYPE)
    case 'hamming'
        win = hamming(WINDOW_SIZE);
    case 'hann'
        win = hann(WINDOW_SIZE);
    case 'hanning'
        win = hanning(WINDOW_SIZE);
    case 'blackman'
        win = blackman(WINDOW_SIZE);
    case 'barthannwin'
        win = barthannwin(WINDOW_SIZE);
    case 'blackmanharris'
        win = blackmanharris(WINDOW_SIZE);
    case 'bohmanwin'
        win = bohmanwin(WINDOW_SIZE);
    case 'chebwin'
        win = chebwin(WINDOW_SIZE);
    case 'gausswin'
        win = gausswin(WINDOW_SIZE);
    case 'kaiser'
        win = kaiser(WINDOW_SIZE);
    case 'nuttallwin'
        win = nuttallwin(WINDOW_SIZE);
    case 'parzenwin'
        win = parzenwin(WINDOW_SIZE);
    case 'tukeywin'
        win = tukeywin(WINDOW_SIZE);
    case 'rect'
        win = ones(WINDOW_SIZE,1);
    case 'flattopwin'
        win = flattopwin(WINDOW_SIZE);
    case 'welch'
         %win = welchwin(WINDOW_SIZE, 0);
        Nd2=(WINDOW_SIZE-1)/2;
        for i=1 : WINDOW_SIZE
            win(i)=1-((i-Nd2)/Nd2)^2;
        end
        win = win';
    otherwise
        disp('WARNING - Unknown window type!');
        disp('Switching to default window type : Hamming Window');
        win = hamming(WINDOW_SIZE);
end

if s_data(1) < WINDOW_SIZE
    error('Cospectra computation: data segment too short for this window length, consider changing frequency resolution');
end


% Calculate the number of windows because of the overlapping
if(OVERLAP == 0)
    number_windows = floor(s_data(1) / WINDOW_SIZE);
else
    nbFullWin = floor(s_data(1)/WINDOW_SIZE); 
    number_windows = 1 + (nbFullWin-1)/(1-OVERLAP) + floor((s_data(1)-((nbFullWin)*WINDOW_SIZE))/((1-OVERLAP)*WINDOW_SIZE));
end


% pre-allocation of memory 
S = zeros(WINDOW_SIZE, s_data(2)); 

%% Loop on all time windows
for window_ix = 1 : number_windows
 
    % time markers to select the data
    t1 = floor((window_ix-1) * (1-OVERLAP) * WINDOW_SIZE) +1;   % marker of the beginning of the time window
    t2 = t1 + WINDOW_SIZE -1;                           % marker of the end of the time window
    
    % select current window and apodize it   
    cdata = DATA(t1:t2, :) .* (win*ones(1,s_data(2)));

    % FFT calculation
    fdata = fft(cdata,WINDOW_SIZE) ./ number_windows; % complex data
    
    % Averaging over non-masked windows
    if MASK(t1:t2)  % if at least one mask sample is set to zero, do not consider this period
        S = S + fdata;     
    end
end

%% Keep log magnitude only and make frequency axis
PSD = 20*log10(abs(S)./number_windows);
PSD = PSD(1:WINDOW_SIZE/2,:);
F_AXIS = FS*(1/WINDOW_SIZE:1/WINDOW_SIZE:.5);

%% Adjust Frequency range to specified range (in case it is a parameter)
if ~isempty(FREQUENCY_RANGE)     
    if  FREQUENCY_RANGE(3) < (FREQUENCY_RANGE(1) / 2)        
        Faxis = [0:(1/number_freqs):1-(1/number_freqs)];
        FindexMin = find(Faxis >= FREQUENCY_RANGE(2)/(FREQUENCY_RANGE(1)/2), 1);
        FindexMax = find(Faxis >= FREQUENCY_RANGE(3)/(FREQUENCY_RANGE(1)/2), 1);
        PSD = PSD(FindexMin:FindexMax,:) ;
    else         
        % Assertion on frequency, to comply with Nyquist Theorem
        error('Frequency range exceeds data sampling rate in consideration of Shannon''s limitation.')
    end
end

