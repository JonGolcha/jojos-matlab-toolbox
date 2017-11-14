function [PHASE] = proc_phase(DATA, SRATE, BP_RANGE, BP_ORDER)
% [PHASE] = proc_phase(DATA, SRATE, FILTSPEC)
%
% Computes the instantaneous phase using bandpass filtering and Hilbert
% transform. 
%
% Inputs:
% - DATA    	--> 2D matrix with the data (P channels by N samples)
%                       OR 3D matrix (P channels, N samples, T trials).
% - SRATE    	--> sampling frequency
% - BP_RANGE    --> bandpass filter range (2 entry vector with limits of the frequency band)
%                       for example, BP_RANGE = [35 45] for gamma band.
% - BP_ORDER 	--> filter order (scalar) OPTIONAL
%                   
%
% Outputs:
% - PHASE    	--> 2D or 3D matrix of instantaneous phase values 
%                   (P channels by N samples)
%
%
% History
% Last version:  07/10/2016
% Created by J.Chatel-Goldman, jonas.chatel.goldman(at)gmail.com
% inspired by implementation from Praneeth Namburi

% assess inputs
if(nargin < 2)  
    error('srate and filter BP must be provided'),
end
P = size(DATA, 1);
N = size(DATA, 2);
if ndims(DATA)==3
    T = size(DATA, 3);
end


% if not provided, calculate the order from the parameters using FIRPMORD.
if (nargin < 3) || isempty(BP_ORDER)
    Fstop1 = BP_RANGE(1);
    Fstop2 = BP_RANGE(2);
    Fpass1  = Fstop1+.5;        % First Passband Frequency 
    Fpass2  = Fstop2-.5;        % Second Passband Frequency
    Dstop1  = 1e-005;           % First Stopband Attenuation
    Dpass   = 0.005;            % Passband Ripple
    Dstop2  = 1e-005;           % Second Stopband Attenuation
    dens    = 20;               % Density Factor
    [BP_ORDER, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(SRATE/2), [0 1 0], [Dstop1 Dpass Dstop2]);
	b  = firpm(N, Fo, Ao, W, {dens});
    disp(['estimated filter order for this freq band (' int2str(Fstop1) '-' int2str(Fstop2) 'Hz) is N=' int2str(BP_ORDER)]),
else
    b  = fir1(BP_ORDER, 2/SRATE*BP_RANGE);
end

% Estimate the filter coefficients
a = 1;
% H_bp = dfilt.dffir(b);
% fvtool(H_bp,'Fs',128,'Analysis','freq','Legend','on') % show filter characterist

% zero-phase digital filtering using filtfilt.
% disp('Filtering data...');
DATA = double(DATA);
if ndims(DATA)==3
    filteredData = zeros(P,N,T);
    for trial_ix =1:T
        filteredData(:,:,trial_ix) = filtfilt(b,a,squeeze(DATA(:,:,trial_ix))')';
    end
else
    filteredData = filtfilt(b,a,DATA')';
end
% figure, plot(1:1250,[squeeze(DATA(1,:,2));squeeze(filteredData(1,:,2))])

% compute phase from Hilbert transform
if ndims(DATA)==3
    PHASE = zeros(P,N,T);
else
    PHASE = zeros(P,N);
end
for chan_ix = 1:P
    PHASE(chan_ix, :, :) = angle(hilbert(squeeze(filteredData(chan_ix, :, :))));
end

return;