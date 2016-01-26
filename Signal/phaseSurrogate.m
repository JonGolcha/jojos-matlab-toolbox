function S=phaseSurrogate(X,SHUFF)
% function S=phaseSurrogate(X,SHUFF)
%
% ************************************************************************
% Creates N phase randomized surrogates of multiple time series X.
% FAST AND EFFICIENT IMPLEMENTATION!
%
% ------- Input:
% X  	--> original timeseries (T samples x N series)
% SHUFF	--> number of surrogates
%
% ------- Output:
% S   	--> phase randomized time series (T x N x SHUFF)
%
% ------- Reference:
% Theiler J, Galdrikian B, Longtin A, Eubank S, Farmer D J (1992): Using 
% Surrogate Data to Detect Nonlinearity in Time Series. In Nonlinear Modeling
% and Forecasting, eds. Casdagli M & Eubank S. 163-188. Addison-Wesley
%
% ------- History:
% Last modified 2014/07/14
% inspired from code by Alexandros Leontitsis (phaseran function)
% J.Chatel-Goldman @OIST --> jonas.chatel.goldman(at)gmail.com


% Assess inputs
if nargin<1 | isempty(X)
   error('[phaseSurrogate] you should provide a time serie!');
elseif size(X,1) == 1
      X = X'; % if only 1 serie, then X should be a column vector
end
if nargin<2 | isempty(SHUFF)==1
   error('[phaseSurrogate] bad number of surrogates')
end
if ~isscalar(SHUFF)
  error('[phaseSurrogate] SHUFF must be scalar.');
end

% Various init
[T N]  	= size(X);        
X_fft   = fft(X);           % Fourier transform
X_mag   = abs(X_fft);     	% magnitudes
X_phase = angle(X_fft);   	% phase
i       = sqrt(-1);      	% imaginary unit
T_half  = floor(T/2);       % half of the data points
randPhase = zeros(T,N,SHUFF);

% Randomizing phases
if rem(T,2)==0   % case signal has even length
    dummy = rand(T_half-1,N,SHUFF)*2*pi;
    randPhase(2:T,:,:)= cat(1, dummy, repmat(X_phase(T_half+1,:),[1 1 SHUFF]) , -flipdim(dummy,1));
else            % case signal has odd length
    dummy = rand(T_half,N,SHUFF)*2*pi;
    randPhase(2:T,:,:)= cat(1, dummy, -flipdim(dummy,1));
end

% Back to the complex numbers
dummy=repmat(X_mag,[1 1 SHUFF]).*exp(i*randPhase);

% Back to the time series (phase randomized surrogates)
S=real(ifft(dummy)); 	






