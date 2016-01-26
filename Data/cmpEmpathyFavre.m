function [SCORES] = cmpEmpathyFavre(RAW_SCORES)
% function [SCORES] = cmpEmpathyFavre(RAW_SCORES)
% -------------------------------------------------------------------------
% Analysis of dual-EEG dual-empathy experiment data (Gispa-Lab 2012)
% -------------------------------------------------------------------------
%
% This function computes empathy scores from answers to the 36 items that 
% compose Favre's questionnaire. It returns scores for 3 dimensions: 
% emotion contagion, empathy, and emotion cut.
%
% *** Inputs:
% - RAW_SCORES  --> vector of length 36 with integer values in range [0:4] 
%
% *** Outputs:
% - SCORES      --> structure with scores for each dimension + global one.
%
% *** References:
% Favre, Joly, Reynaud, & Salvador, 2009
%
% *** History:
% --- 2013-07-16
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 

% transpose to range [0 1]
RAW_SCORES  = double(RAW_SCORES)/4;

% sum over 12 item for each dimension
SCORES.empathy     = sum(RAW_SCORES([1 6 9 14 16 17 19 21 27 28 34 36]));
SCORES.contagion   = sum(RAW_SCORES([2 4 5 7 12 13 18 20 23 24 32 35]));
SCORES.emoCut      = sum(RAW_SCORES([3 8 10 11 15 22 25 26 29 30 31 33]));
SCORES.global      = 2*SCORES.empathy + 1*SCORES.contagion - 1*SCORES.emoCut; % this wheighting is arbitrary

end

