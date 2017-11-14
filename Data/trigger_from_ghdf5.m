function [trigVal trigIndex trigTime] = trigger_from_ghdf5(DATA)
% -------------------------------------------------------------------------
% This function extract trigger times and values from a matlab structure
% imported from a ghdf5 g.tec EEG file. 
% IMPORTANT: triggers must be recorded using rising edges only, resetting triggers with
% long enough inter-trigger delay.
%
% Inputs:
% DATA      --> matlab data structure imported from ghdf5 file
% 
% Outputs
% trigVal   --> trigger values (decimal)
% trigIndex --> trigger sample index
% trigTime  --> trigger time values (in seconds)
%
% History:
% --- 2012-06-01
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 



% Each time a bit changes (here rising only), the triggers save the
% sample index in Time vector, and the new bit value in decimal basis.
% To recover the full markers' type and time position, one should
% reconstruct unique triggers from decomposed bits times and values.
Fs = DATA.RawData.AcquisitionTaskDescription.SamplingFrequency;
[UniqueTrig, Unique_m] = unique(DATA.AsynchronData.Time,'first');
NbTrig = length(UniqueTrig);
trigIndex =  double(UniqueTrig');
trigTime = trigIndex / Fs;        
trigVal = zeros(NbTrig,1);
for i = 1:NbTrig-1
    trigVal(i) = sum(DATA.AsynchronData.Value(Unique_m(i):Unique_m(i+1)-1));
end
trigVal(end) = sum(DATA.AsynchronData.Value(Unique_m(end):end));





