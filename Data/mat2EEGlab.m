% -------------------------------------------------------------------------
% Analysis of dual-EEG dual-empathy experiment data (Gispa-Lab 2012)
% -------------------------------------------------------------------------
%
%
% This script computes EEG-lab .set structure from raw EEG data stored in
% .mat format.
% 
% History:
% --- 2013-04-23
% Created by Jonas Chatel-Goldman @ GIPSA-Lab 


clear all
close all
addpath(genpath([fileparts(pwd) '\FUNCTIONS_JCG'])); % add all subfolders below main functions folder
cfgANA.DataPath = 'C:\Documents and Settings\jonas\Mes documents\Dual_empathy\DATA\EEG-LAB_STUDY\';


% -------------------------- Parametrization
% Analysis parameters
cfgANA.coupleID         = 'verybest';   % can be 'all', 'best', 'verybest', 'eyesopen', 'eyesclosed'  or a vector with couples IDs
cfgANA.dataType         = 'EEG_only';   % 'all' or 'EEG_only' or 'EXT_only' : choice of data on which to perform analysis 

% Data importation parameters
cfgANA.manualArtifRej   = 0;            % '1' to consider manual artifact rejection from cartool marker file (i.e. remove these periods from .set files)
cfgANA.trialLim         = [20 120];   	% put 'all' or trial [start end] in seconds, useful for removing trial begin and/or end.
cfgANA.nbSmoothSamples  = 4;            % specifies nb of samples to smooth at begin and end of artifact-free periods

% EEG data parameters
cfgEEG.Fs               = 128 ;         % sampling rate (assumed identical for each subject)
cfgEEG.NbChannels       = 46;           % number of channels per subject (assumed identical for each subject)
cfgEEG.TrialLength      = 120;          % Trial length in seconds
cfgEEG.chan_str         = { 'FP1','FP2','AFZ','F7','F3','FZ','F4','F8','T7',...
                            'C3','CZ','C4','T8','P7','P3','PZ','P4','P8',...
                            'O1','O2','Resp','Pulse','GSR','FP1','FP2','AFZ', ...
                            'F7','F3','FZ','F4','F8','T7','C3','CZ','C4','T8',...
                            'P7','P3','PZ','P4','P8','O1','O2','Resp', 'Pulse','GSR'};
cfgEEG.EEG_channels     = [1:20,24:43]; % indexes of EEG channels
cfgEEG.EXT_channels     = [21:23,44:46];% indexes of external channels (Resp, Pulse, GSR)



%% --- Various initializations
% RawDataPath contains path with couples' EEG raw .mat files, (root folder
% from where STUDY files are stored)
cfgANA.RawDataPath = [ fileparts(cfgANA.DataPath(1:end-1)) '\']; 

% Here we initialize couples ID based on experimental priors
if strcmp(cfgANA.coupleID ,'all')    
    cfgANA.coupleID = [2:5,7:16];              % keep all couple data (couple n°6 has no data)
elseif strcmp(cfgANA.coupleID ,'best')
    cfgANA.coupleID = [7 9 10 11 13 14 15 16]; % choice of best couples is based on subjective experimental assessment
elseif strcmp(cfgANA.coupleID ,'verybest')
    cfgANA.coupleID = [9 10 14 15 16];         % choice of best couples is based on subjective experimental assessment
elseif strcmp(cfgANA.coupleID ,'eyesopen')
    cfgANA.coupleID = [2 3 4 5 7];             % subjects with eyes open during the experiment
elseif strcmp(cfgANA.coupleID ,'eyesclosed')
    cfgANA.coupleID = [9 10 11 13 14 15 16];   % subjects with eyes open during the experiment
end

% Data segmentation
cfgEEG.CondTrigVal      = [8:12,16:20]; % Trigger indexes used for all conditions in the experiment
cfgEEG.NbCondTot        = length(cfgEEG.CondTrigVal);
effect_name = { 'Touch-Congruent' ; 'Touch-noCongruent' ; 'Touch-Control' ; ...
                'noTouch-Congruent' ; 'noTouch-noCongruent' ; 'noTouch-Control'};
effect_cond_ix = {[17 20] ; [18 19] ; 16; [9 12] ; [10 11] ; [8]};  

% Channel selection
if strcmp(cfgANA.dataType,'EEG_only')
    cfgEEG.chan_str = cfgEEG.chan_str(cfgEEG.EEG_channels);
    cfgEEG.chanSel = cfgEEG.EEG_channels;
elseif strcmp(cfgANA.dataType,'EXT_only')
    cfgEEG.chan_str = cfgEEG.chan_str(cfgEEG.EXT_channels);
    cfgEEG.chanSel = cfgEEG.EXT_channels;
end
cfgEEG.nbChan  = length(cfgEEG.chanSel);           % total number of channels
cfgEEG.S       = [cfgEEG.nbChan/2  ; cfgEEG.nbChan/2];   % size of each subject channel sets
cfgEEG.chanSet = [1:cfgEEG.nbChan/2 ; cfgEEG.nbChan/2+1:cfgEEG.nbChan];  % subject1 = 1st row, subject2 = 2nd row

% read electrode coordinates
try
    cfgEEG.eloc = readlocs('EEGLAB_electrodes_Dual_Empathy.loc','filetype','loc'); 
catch
    error('You forgot to launch eeglab first! (again...)');  % (do not forget to launch EEGlab first)
end


fprintf('\n \n');


%% --- Create EEGLAB structures
        
% -- Creating one .mat file with raw EEG for each subject / condition
% here we use data that was imported first using 'Import_gTEC_data.m' 

for couple_ix = cfgANA.coupleID

    % load data
    dataFile = ['couple' int2str(couple_ix) '_raw.mat'];
    dataPath = [cfgANA.RawDataPath '\Couple_' int2str(couple_ix) '\' ];
    disp(['Loading current datafile:    ' dataFile '...']);
    try
        load([dataPath dataFile]); 
    catch
        error(['INEXISTANT DATA FILE: ' dataFile]); 
    end

    % keep only selected data type (EEG, EXT or both)
    EEG.data = EEG.data(cfgEEG.chanSel,:);

    % load markers for manual artifact rejection
    if cfgANA.manualArtifRej == 1
        mrkFile = [dataPath 'couple' int2str(couple_ix) '_filtered.ep.MANUAL-ARTIF-REJ.mrk'];
        if ~exist(mrkFile,'file')
            error('Marker file for manual artifact rejection does not exist!');
        end
        EEG.mask_mrk  = import_cartool_mrk(mrkFile, size(EEG.data,2));
        disp('Marker file loaded for manual artifact rejection.');
        EEG.lengthOK = length(EEG.mask_mrk(EEG.mask_mrk==1)); % in samples
    else
        EEG.mask_mrk  = true(1, size(EEG.data,2));    % select everything
        disp('No manual artifact rejection (everything selected).')
    end

    % format trigger to specified conditions and trial begin/end
    % - EEG.trigVal(:,1) -> trigger name (condition id)
    % - EEG.trigVal(:,2) -> start sample 
    % - EEG.trigVal(:,3) -> end sample
    % - EEG.trigIndex    -> cell array for indexing 'EEG.data' to condition
    %                       periods, {nb conditions} ( [1,trialLength] , ID )
    % - EEG.maskedTrigIndex -> idem with mask set using manual artifact rej
    % - EEG.mask_mrk_inCond -> mask with artifact free in condition periods
    startOffset         = cfgANA.trialLim(1) * cfgEEG.Fs ;
    endOffset           = cfgANA.trialLim(2) * cfgEEG.Fs ;
    EEG.trialLength     = cfgEEG.Fs*(cfgANA.trialLim(2)-cfgANA.trialLim(1)) ;  % data length per trial in samples
    tempTrigVal         = EEG.trigVal(:,2);
    EEG.trigVal(:,2)    = tempTrigVal + startOffset;
    EEG.trigVal(:,3)    = tempTrigVal + endOffset -1;
    condToKeep          = findVecInVec(EEG.trigVal(:,1), cfgEEG.CondTrigVal); 
    EEG.trigVal         = EEG.trigVal(condToKeep,:);  % this rejects inter trial rest periods
    EEG.trigIndex     	= zeros(cfgEEG.NbCondTot, EEG.trialLength); 
    EEG.trigIndex    	= num2cell(EEG.trigIndex,2);
    EEG.maskedTrigIndex = zeros(cfgEEG.NbCondTot, EEG.trialLength); 
    EEG.maskedTrigIndex = num2cell(EEG.trigIndex,2);
    EEG.mask_mrk_inCond = false(1,size(EEG.data,2));
    for cond_ix = 1:cfgEEG.NbCondTot
        EEG.trigIndex{cond_ix}      = EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3); % EEG.trigVal(cond_ix,2) + (EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3));
        EEG.trigIndex{cond_ix,2}    = EEG.trigVal(cond_ix,1) ;      % condition ID
        EEG.maskedTrigIndex{cond_ix}= EEG.trigVal(cond_ix,2) + ...
                                            find( EEG.mask_mrk(EEG.trigVal(cond_ix,2):EEG.trigVal(cond_ix,3)) );
        EEG.maskedTrigIndex{cond_ix,2} = EEG.trigVal(cond_ix,1) ;   % condition ID     
        EEG.mask_mrk_inCond(EEG.maskedTrigIndex{cond_ix}) = true;
    end
    EEG.lengthOK = length(find(EEG.mask_mrk_inCond)); % in samples
    if cfgANA.manualArtifRej == 1
        % weight EEG with smoothing windows to remove time discontinuities when selecting artifact periods
        EEG.mask_mrk_smooth = smoothMask(EEG.mask_mrk,cfgANA.nbSmoothSamples,2);	
        EEG.data = EEG.data .* (ones(cfgEEG.nbChan,1) * EEG.mask_mrk_smooth);          % smoothing EEG at edges of artifact periods
        disp([  'Total data length kept for analysis: ' int2str(EEG.lengthOK/cfgEEG.Fs) ...
            ' seconds (' int2str(100*EEG.lengthOK/(cfgEEG.NbCondTot*cfgEEG.TrialLength*cfgEEG.Fs)) '%)']);
    end

    % perform EEG data segmentation to condition periods
    for effect_ix = 1:size(effect_name,1)
        temp_ix     = ismember(cell2mat(EEG.maskedTrigIndex(:,2)),effect_cond_ix{effect_ix});   
        mask_effect = cat(2,EEG.maskedTrigIndex{temp_ix,1});
        for sub_ix = 1:length(cfgEEG.S)  
            % write segmented .mat file for each subject to disk 
            EEG_effect  = EEG.data(cfgEEG.chanSet(sub_ix,:),mask_effect);
            outFileName = ['couple' int2str(couple_ix) '_sub' int2str(sub_ix) '_' effect_name{effect_ix} '_raw.mat'];
            outFilePath = [cfgANA.DataPath 'couple' int2str(couple_ix) '_sub' int2str(sub_ix) '\'];
            disp(['Writing ' outFileName '...']);
            if ~isdir(outFilePath)
                mkdir(outFilePath);
            end
            save([outFilePath outFileName], 'EEG_effect');


            % import data to eeglab .set format
            subName = ['couple' int2str(couple_ix) '_sub' int2str(sub_ix) ];
            setName = [subName '_' effect_name{effect_ix}];
            EEGOUT	= pop_importdata(   'setname',setName,'data',EEG_effect,'dataformat','array',...
                                        'subject',subName,'condition',effect_name{effect_ix},...
                                        'group',['couple' int2str(couple_ix)],'chanlocs',cfgEEG.eloc,...
                                        'nbchan',cfgEEG.S(sub_ix),'srate',cfgEEG.Fs);

            % store eeglab .set files in disk    
            outFileName = ['couple' int2str(couple_ix) '_sub' int2str(sub_ix) '_' effect_name{effect_ix} '_raw.set'];
            pop_saveset( EEGOUT, 'filename',outFileName,'filepath',outFilePath); 
        end  
    end        
end

%     % -- Finally, create a big STUDY dataset
%     [ STUDY ALLEEG ] = pop_study([], [], 'gui', 'on'); 
%     % call it e.g. Couples_9-10-14-15-16_noArtifRej

clear dataFile couple_ix EEG startOffset endOffset condToKeep temp_ix setName subName
clear mrkFile tempTrigVal cond_ix couple_ix dataFile dataPath effect_ix mask_effect
clear sub_ix EEG_effect outFileName outFilePath effect_name effect_cond_ix






















