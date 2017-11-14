% function [OUTEEG] = pop_prop_landscape2(EEG, varargin)
function [OUTEEG] = pop_prop_landscape2(EEG,varargin)
% 

% TODO: 
% - connect selection rect with labels and button selection
% - manage mouse click events and connect with rectangles selection ginput
% - output updated IC labelling structure



% *** Description:
% Display a summary image ('landscape') of component properties for 20 ICs:
% - scalp topography
% - ERP image
% - component time series (selected window)
% - activity power spectrum
% - IC number, Percent Data Variance Accounted For, 
% - Residual Variance of dipole fitting
% 
% In addition, allows to manually label MULTIPLE components at once. 
% Possible labels are: 
% - 'brain'
% - 'eyeblink'
% - 'eyemov'
% - 'myogenic'
% - 'cardiac'
% - 'line'
% - 'channel'
% - 'other'

% Each time a label is set, it is stored in cell array 'EEG.icalabels'.
% If field 'EEG.icalabels' already exist, it is used to init checkboxes.
% Output OUTEEG structure is saved automatically upon figure closing.
%
% *** Usage (example):
%   >> EEG = pop_prop_landscape( EEG );
%   >> pop_prop_landscape( EEG,'icnum',1,'spec_opt',{'freqrange',[0 50]}); 
%
% *** Required command line input:
% - EEG     	--> EEGLAB dataset structure.
%
% *** Optionnal command line keywords:
% - 'spec_opt'	--> cell array of options for spectopo().
% - 'ID'        --> number or string for the dataset ID.
% - 'UI_control'--> 'interactive' to activate user interface controls to 
%                   scroll and annotate ICs (default). 
%                   'scroll' just display without user input for labelling
%                   'no_ui' only displays one component without any UI.
% - 'reset_labels' --> '1' to reset labels (default: '0')
%
% % *** Note:
% For the dipole plot, you need EEG.dipfit with one dipole and 
% EEG.dipfit2 with two dipoles. After setting dipfit, this can be obtained
% as follows:
%   EEG = pop_dipfit_settings( EEG, dipfit_settings_opts{:});
%   EEG = pop_multifit(EEG, [], 'dipoles', 1, 'threshold', 100);
% 	EEG_temp = pop_multifit(EEG, [], 'dipoles', 2, 'threshold', 100);
%  	EEG.dipfit2 = EEG_temp.dipfit;
%
% % *** Note 2:
% A bug was resolved in dipplot in line 1161 (see implementation at the end
% of this file).
%
% % *** Note 3:
% To speed display x100, function  erpimage() was modified as follows:
% line 3475, commented: axcopy(gcf);
% 
%  *** History
% Jonas Chatel-Goldman (2016)


%% initialize
% assess inputs
% if EEG structure is a handle, we are calling the function within itself (callback)
if ishandle(EEG)
    input   = varargin;
    EEG     = input{1};
    g       = input{2};
    currPanel = input{3};
    disp_landscape(currPanel);
   
% normal call    
else
    g.nIC_disp = 15; % number of component to disply simultaneously
    currPanel = 1; 
    try
       options = varargin;
       if ~isempty( varargin ), 
           for i = 1:2:numel(options)
               g.(options{i}) = options{i+1};
           end
       else g= []; end;
    catch
       disp('pop_prop_landscape() error: calling convention {''key'', value, ... } error'); return;
    end;
    try, g.spec_opt; 	catch, g.spec_opt   = {'freqrange',[0 100]};    end;
    try, g.ID;          catch, g.ID         = 'NO_ID';                  end;
    try, g.UI_control;  catch, g.UI_control = 'interactive';            end;
    try, g.reset_labels; catch, g.reset_labels = false;                 end;
    try, g.long_trial;  catch, g.long_trial = 'off';                    end;
    try, EEG.icaweights; catch, error('ICA must be computed first');    end;
    if isempty(EEG.chanlocs)
        error('channels location must be specified (field EEG.chanlocs is empty)'); 
    end
    g.nbIC = size(EEG.icaweights,1); % total number of IC in this dataset
    g.selected_ICs = false(1,g.nbIC);

    % reset labels if required
    if g.reset_labels || ~isfield(EEG,'icalabels') || ~iscell(EEG.icalabels)
        EEG.icalabels = cell(g.nbIC,1);
        % example: 
        % EEG.icalabels{1} = {'brain'};
        % EEG.icalabels{2} = {'other','not_sure'}
    end
    % recompute ICA activations if not present
    if isempty(EEG.icaact) 
        try
            EEG.icaact = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights,1)); 
        catch
            error('Could not obtain IC time courses. Was ICA already computed?');
        end
    end
    % recompute all power spectra if not present
    if ~isfield(EEG,'icaspec') || isempty(EEG.icaspec) || isempty(EEG.freqs) 
        [EEG.icaspec,EEG.freqs] = spectopo(EEG.icaact, 0, EEG.srate, 'overlap', 128, 'percent', 100, 'plot', 'off');
    end
    % recompute all ERP images if not present
    if EEG.trials>1 && (~isfield(EEG,'ERP') || isempty(EEG.ERP))
        EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
        ei_smooth = 1;
        EEG.ERP = cell(g.nbIC,1);
        for IC_ix = 1:g.nbIC
            icaacttmp = eeg_getdatact(EEG, 'component', IC_ix) * norm(EEG.icawinv(:,IC_ix));
            offset = nan_mean(icaacttmp(:));
            era    = nan_mean(squeeze(icaacttmp)')-offset;
            era_limits = get_era_limits(era);
            EEG.ERP{IC_ix} = erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
                '', ei_smooth, 1, 'caxis', 2/3, 'cbar','off', 'yerplabel', '','erp_vltg_ticks',era_limits,'noxlabel','NoShow','yes'); % strange automatic threshold for epochs > 10s with this option: 'limits',[g.epoch_lims(1) g.epoch_lims(2)]           
        end
    end

    % initialize figure
    basename = ['Component_Landscape'];
    g.BACKCOLOR = [.93 .96 1];
    g.hdl.mainfig = figure('name', ['SetID_' num2str(g.ID) '_' basename],...
        'color', g.BACKCOLOR,...
        'numbertitle', 'off',...
        'visible', 'on',...
        'PaperPositionMode','auto',...
        'ToolBar', 'none',...
        'MenuBar','none', ...
        'CloseRequestFcn', @closereq_cb);

    
    set(g.hdl.mainfig,'WindowStyle','docked');

    % store user information for subsequent callbacks
    gui_handle = guihandles(g.hdl.mainfig); 
    gui_handle.EEG = EEG;
    gui_handle.g = g;
    guidata(g.hdl.mainfig,gui_handle)     

    % initialize UI controls
    if ~strcmp(g.UI_control,'no_ui')
        nPanels = ceil(g.nbIC/g.nIC_disp); % number of panels necessary to display all components
        if nPanels>1
            g.hdl.sliderSelect  = uicontrol(g.hdl.mainfig,'style','slider','Min',1,'Max',nPanels,'SliderStep',[1/(nPanels-1),5/(nPanels-1)],...
                                'units','normalized','Position',[.25 .01 .5 .05], 'BackgroundColor','w', 'HandleVisibility', 'off', ...
                                'Value',currPanel,'HandleVisibility','off','Callback',@(hObject,callbackdata)slider_cb());  
        end
        % take key arrows as input to scroll components THIS DOESN'T WORK (TO FIX)
        % set(g.hdl.mainfig, 'WindowKeyPressFcn', @(hObject,evt)readkey_cb(evt.Key));
        % set(g.hdl.mainfig, 'WindowKeyPressFcn', @readkey_cb);
    end

    % Call display function
    disp_landscape(currPanel);

    % Pause until figure is closed (helas mandatory to output EEG struct!)
    waitfor(g.hdl.mainfig);   
    
end


%% Display component properties landscape
function disp_landscape(currPanel) % nested function
    
    % pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    EEG = gui_handle.EEG;
    
    % IC to scroll
    g.IC_scroll = 1+g.nIC_disp*(currPanel-1):g.nIC_disp*currPanel;
    g.IC_scroll(g.IC_scroll>g.nbIC) = []; 
    
    % erase previous plots
    clf(g.hdl.mainfig);
    
    % show ui control for manual labelling
    if strcmp(g.UI_control,'interactive') 
        g = disp_uicontrol_interactive(g);
    end
    
    disp('***********');
    disp(['DISPLAYING COMPONENTS PROPERTIES: IC' int2str(g.IC_scroll(1)) '-' int2str(g.IC_scroll(end))]);
   
    % plot properties for each IC
    x = 3; % nb lines
    y = 5; % nb columns   
    g.plotSel = false(1,length(g.IC_scroll));
    for plot_ix = 1:length(g.IC_scroll)
        current_IC = g.IC_scroll(plot_ix);
        g.hdl.h_sub(plot_ix) = subplot(x,y,plot_ix); % create subplot       
        cla(g.hdl.h_sub(plot_ix)); 
        set(g.hdl.h_sub(plot_ix),'Visible','off','Box','on');
   
        
        % plot erpimage
        if EEG.trials > 1 % epoched data
            h_erp = axes; % dummy axis necessary
            mindat = min(min(EEG.ERP{current_IC}));
            maxdat = max(max(EEG.ERP{current_IC}));
            maxdat =  max(abs([mindat maxdat])); % make symmetrical about 0
            mindat = -maxdat;
            imagesc(EEG.ERP{current_IC}',[mindat,maxdat]); hold on;
            colormap jet;
            set(gca,'XTick',[],'YTick',[],'Ylabel',[]);
            set(gca,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 .6 .6] + [-.1 -.02 0 0]); % reduced position in lower left quarter
            %annotation('line',[? ?],[? ?],'Color','k','Linewidth',2); % don't work?...
            %line([0 0],[0 size(EEG.ERP{current_IC},2)],'Color','k','Linewidth',2); % don't work?...
            
        end
        
        
        % plot spectrum
      	h_spec = axes; % dummy axis necessary
        plot(EEG.freqs,EEG.icaspec(current_IC,:),'color','r','Linewidth',2);
        axis on tight; box on; grid on;
        set(gca,'XTick',[],'YTick',[],'Ylabel',[],'Xlabel',[]);
        xlim(gca, g.spec_opt{2});
        set(h_spec,'XTick',[],'YTick',[],'Ylabel',[]);
        set(h_spec,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 .6 .6] + [-.02 -.02 0 0]);  % reduced position in lower right quarter
        
        
        % plot scalp map
      	h_scalp = axes; % dummy axis necessary
        scalpmap_norm = EEG.icawinv(:,current_IC)/sqrt(sum(EEG.icawinv(:,current_IC).^2));
        [htopo,Zi,plotrad] = topoplot( scalpmap_norm, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
            'shading', 'interp', 'numcontour', 0,'electrodes','on'); 
        set(h_scalp,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 .6 .6] + [-.1 .12 0 0]);   % reduced position in upper left quarter

        % plot time series
        [dat, time] = choose_ts_seg(EEG,current_IC);
        h_time = axes; % dummy axis necessary
        plot(time,dat'); 
        set(h_time,'Ylim',[-max(abs(dat)) max(abs(dat))],'XTick',[],'YTick',[],'Ylabel',[]); %doesn't work? 'Box','on','XGrid','on','YGrid','on',
        set(h_time,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 .6 .6] + [-.02 .12 0 0]);   % reduced position in upper right quarter

        
        % add text information: icnum and pvaf
        try
            pvaf = num2str(EEG.pvaf(current_IC), '%4.2f');
        catch
            pvaf = 'N/A';
        end
        str_id = ['IC#' int2str(current_IC) '/' int2str(g.nbIC)];
        str_infos = ['RV: ' num2str(EEG.dipfit.model(current_IC).rv*100, '%.1f') '%  |  pvaf: ' pvaf '%'];
        h_text = axes;
        text(0,1,str_id,'Units', 'Normalized','Fontsize', 8,'Interpreter','none','BackgroundColor','w','Margin',1);     
        set(h_text,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 0 0] + [-.04 0.12 0 0]);   % reduced position in upper right quarter
        h_text = axes; 
        text(0,1,str_infos,'Units', 'Normalized','Fontsize', 8,'Interpreter','none');
        set(h_text,'Position',get(g.hdl.h_sub(plot_ix),'Position') .* [1 1 0 0] + [-.075 -.03 0 0]);   % reduced position in upper right quarter
                
        % create selection rectangle
        currPos = g.hdl.h_sub(plot_ix).Position + [-.1 -.045 .03 .08];
        g.hdl.rec_Sel(plot_ix) = annotation('rectangle',currPos,'PickableParts','all','ButtonDownFcn',@(hObject,callbackdata)click_selectionRect(plot_ix),'EdgeColor','none');
       
        
        % draw that
        drawnow;
    end
        
    % push data to current figure
    gui_handle.EEG = EEG;
    gui_handle.g = g;
    guidata(gcf,gui_handle)  
    
    % display current IC selection (must be done AFTER data push)
  	disp_selectionRect();
    
    % TEMP 
%     OUTEEG = EEG;   
end

% --- Close request (nested) function, executes when closing the figure window
function closereq_cb(hObject,eventdata) 
    gui_handle = guidata(gcf);
    OUTEEG = gui_handle.EEG;
    delete(g.hdl.mainfig); % close GUI
end


% --- Slider selection (nested) function, executes when scrolling
function slider_cb()
    current_val = round(get(gcbo,'Value')); 
    gui_handle = guidata(gcf);
    EEG = gui_handle.EEG;
    g = gui_handle.g;
    g.selected_ICs = false(1,g.nbIC); 
    pop_prop_landscape2(gcf,EEG,g,current_val);
end

end


% --- executes when an IC is clicked
function click_selectionRect(plot_ix)
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    EEG = gui_handle.EEG;
    
     % toogle selection 
    if isscalar(plot_ix) 
        current_IC = g.IC_scroll(plot_ix);
        g.plotSel(plot_ix) = ~g.plotSel(plot_ix);
        g.selected_ICs(current_IC) = ~g.selected_ICs(current_IC);
    end
    
    % unselect all checkboxes
    for checkbox_ix = 1:length(g.hdl.checkBox)
        g.hdl.checkBox(checkbox_ix).Value = false; 
    end
   
    % display the intersection of all selected ICs labels
    if sum(g.plotSel)>0
        IC_sel = find(g.selected_ICs);
        labels_inter = EEG.icalabels{IC_sel(1)};
        if length(IC_sel)>1
            IC_sel(1)=[];
            for multiple_sel_ix = IC_sel
                try
                    labels_inter = intersect(labels_inter,EEG.icalabels{multiple_sel_ix});
                end
            end
        end
        if ~isempty(labels_inter)
            [~,ia] = intersect({g.hdl.checkBox(:).Tag},labels_inter);
            for checkbox_ix = ia'
                g.hdl.checkBox(checkbox_ix).Value = true;
            end
        end
    end
    
  	% push new 'g' structure to current figure
    gui_handle.g = g;
    gui_handle.EEG = EEG;
    guidata(gcf,gui_handle);
    
	% display current IC selection (must be done AFTER data push)
  	disp_selectionRect();
end


% --- Executes on button press in checkBoxes
% this function updates labels of currently selected ICs (TODO: VERIFY THAT!)
function checkBox_cb()
    % pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    EEG = gui_handle.EEG;
    currCheckbox = gcbo;
    % update state of checkboxes
    if sum([g.hdl.checkBox(:).Value]) == 0
        g.checkboxState = false;
    else
        g.checkboxState = true;
    end
    % checkbox was just unselected
	if fix(currCheckbox.Value) == 0
        % remove current label to all selected labels
        for IC_sel_ix = find(g.selected_ICs)
            [~,ia] = intersect(EEG.icalabels{IC_sel_ix},currCheckbox.Tag);
            EEG.icalabels{IC_sel_ix}(ia) = [];
        end
	% checkbox was just selected
    else
        g.selected_labels = [g.selected_labels, currCheckbox.Tag]; % add current label to all selected labels
        % add this label to all selected ICs
        for IC_sel_ix = find(g.selected_ICs)
            [C,ia] = intersect(EEG.icalabels{IC_sel_ix},currCheckbox.Tag);
            if isempty(ia)
                EEG.icalabels{IC_sel_ix}(end+1) = {currCheckbox.Tag};
            end
        end
    end
    guidata(gcf,gui_handle);
  	% push data to current figure
    gui_handle.g = g;
    gui_handle.EEG = EEG;
    guidata(gcf,gui_handle);
end


% --- display selection rectangles
function disp_selectionRect()
    % pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    EEG = gui_handle.EEG;
    
    for plot_ix = 1:length(g.IC_scroll)
        current_IC = g.IC_scroll(plot_ix);
        if g.plotSel(plot_ix) % if this IC is selected
          	set(g.hdl.rec_Sel(plot_ix),'LineStyle','-','LineWidth',2,'EdgeColor','k');
        else
            if isempty(EEG.icalabels{current_IC})  % surround in red if it has no label
                set(g.hdl.rec_Sel(plot_ix),'LineStyle','-','LineWidth',1,'EdgeColor','r');
            else    % otherwise do not display selection rectangle
            	set(g.hdl.rec_Sel(plot_ix),'LineStyle','none');
            end
        end
    end
    % push new 'g' structure to current figure
    gui_handle.g = g;
    guidata(gcf,gui_handle);
end

% --- Executes when pushing button select all
function selectAll()
    % first unselect all (so all checkbox are unselected without reseting labels)
    unselectAll();
    % pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    % select all ICs
    g.plotSel = true(1,length(g.plotSel));
    g.selected_ICs = false(1,length(g.selected_ICs));
    g.selected_ICs(g.IC_scroll) = true;
    % push data to current figure
    gui_handle.g = g;
    guidata(gcf,gui_handle);   
    % update display (must be done AFTER data push)
    click_selectionRect(g.selected_ICs);
    disp_selectionRect();
end

% --- Executes when pushing button unselect all
function unselectAll()
    % pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    % unselect all ICs
    g.plotSel = false(1,length(g.plotSel));
    g.selected_ICs = false(1,length(g.selected_ICs));
    % push data to current figure
    gui_handle.g = g;
    guidata(gcf,gui_handle);   
    % unselect all labels
    clearLabels();
    % update display (must be done AFTER data push)
    disp_selectionRect();
end

% --- Executes when pushing button unselect all
function clearLabels()
	% pull data from current figure
    gui_handle = guidata(gcf);
    g = gui_handle.g;
    EEG = gui_handle.EEG;
    % clear all selected labels and unselect all checkboxes
    g.selected_labels=[];  
    for checkbox_ix = 1:length(g.hdl.checkBox)
        g.hdl.checkBox(checkbox_ix).Value = 0;
    end
    % remove all labels for currently selected ICs
    for IC_ix = find(g.selected_ICs)
        EEG.icalabels{IC_ix} = {};
    end
    % push data to current figure
    gui_handle.g = g;
    gui_handle.EEG = EEG;
    guidata(gcf,gui_handle);
end

% --- show ui control for manual labelling
function g_out = disp_uicontrol_interactive(g_in)
    % get g structure
    g = g_in;
    
    % setup checkboxes for manual artifact annotation
    pos_checkBox = [.9 .9 .05 .05]; % position of first left checkbox
    space = .09; % vertical space between checkboxes
    
    % checkbox brain
    tag = 'brain';
    g.hdl.checkBox(1) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                         
    annotation('textbox',[pos_checkBox(2)-.01,pos_checkBox(2)-.01,.05,.02],'String','brain','LineStyle','none')
    % checkbox eye blink
    tag = 'eyeblink';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(2) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.01,pos_checkBox(2)-.01,.05,.02],'String','blink','LineStyle','none')
    % checkbox eye movement
    tag = 'eyemov';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(3) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.015,pos_checkBox(2)-.01,.05,.02],'String','eyemov','LineStyle','none')
    % checkbox cardiac
    tag = 'cardiac';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(4) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.012,pos_checkBox(2)-.01,.05,.02],'String','heart','LineStyle','none')
    % checkbox myogenic
    tag = 'myogenic';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(5) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.013,pos_checkBox(2)-.01,.1,.02],'String','muscle','LineStyle','none')
    % checkbox channel noise
    tag = 'channel';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(6) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.018,pos_checkBox(2)-.01,.1,.02],'String','channel','LineStyle','none')
    % checkbox line
    tag = 'line';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(7) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.008,pos_checkBox(2)-.01,.1,.02],'String','line','LineStyle','none')
    % checkbox other
    tag = 'other';
    pos_checkBox(2) = pos_checkBox(2)-space;
    g.hdl.checkBox(8) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
                                'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
    annotation('textbox',[pos_checkBox(1)-.012,pos_checkBox(2)-.01,.05,.02],'String','other','LineStyle','none')
%     % checkbox not sure
%     tag = 'not_sure';
%     pos_checkBox(2) = pos_checkBox(2)-space;
%     g.hdl.checkBox(9) = uicontrol(g.hdl.mainfig,'style','checkbox','Tag',tag,'Min',0,'Max',1,...
%                                 'units','normalized','Position',pos_checkBox, 'BackgroundColor',g.BACKCOLOR, ...
%                                 'Value',0,'Callback',@(hObject,callbackdata)checkBox_cb());                                 
% %     annotation('textbox',[pos_checkBox(1)-.025,pos_checkBox(2)-.01,.1,.02],'String','(not sure)','LineStyle','none')
    
  	% button unselect all
    pos_button = [.85 .14 .125 .05];
    g.hdl.button(1) = uicontrol(g.hdl.mainfig,'style','pushbutton','units','normalized','String','UNSELECT ALL IC',...
                                'Position',pos_button, 'BackgroundColor',g.BACKCOLOR, ...
                                'Callback',@(hObject,callbackdata)unselectAll()); 

  	% button select all
    pos_button = pos_button + [0 -.065 0 0];
    g.hdl.button(2) = uicontrol(g.hdl.mainfig,'style','pushbutton','units','normalized','String','SELECT ALL IC',...
                                'Position',pos_button, 'BackgroundColor',g.BACKCOLOR, ...
                                'Callback',@(hObject,callbackdata)selectAll()); 
    
   % button clear labels
    pos_button = pos_button + [0 -.065 0 0];
    g.hdl.button(3) = uicontrol(g.hdl.mainfig,'style','pushbutton','units','normalized','String','CLEAR LABELS',...
                                'Position',pos_button, 'BackgroundColor',g.BACKCOLOR, ...
                                'Callback',@(hObject,callbackdata)clearLabels());                          
   
    % create a binary scalar that contains the state of all checkboxes (true if at least one checkboxe is checked)
    % checkboxes = findobj(gcf,'Type','uicontrol','-and','style','checkbox');    % browse all checkboxe type objects
    g.checkboxState = false;  
    g.selected_labels = {};
    
    % save new g structure
    g_out = g;
end


function era_limits=get_era_limits(era)
%function era_limits=get_era_limits(era)
%
% Returns the minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)
%
% Inputs:
% era - [vector] Event related activation or potential
%
% Output:
% era_limits - [min max] minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)

    mn=min(era);
    mx=max(era);
    mn=orderofmag(mn)*round(mn/orderofmag(mn));
    mx=orderofmag(mx)*round(mx/orderofmag(mx));
    era_limits=[mn mx];
end

function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable
    val=abs(val);
    if val>=1
        ord=1;
        val=floor(val/10);
        while val>=1,
            ord=ord*10;
            val=floor(val/10);
        end
        return;
    else
        ord=1/10;
        val=val*10;
        while val<1,
            ord=ord/10;
            val=val*10;
        end
        return;
    end
end

function [dat, time] = choose_ts_seg(EEG,ic)
    srate = EEG.srate;
    pnts = EEG.pnts;
    epoched = EEG.trials>1;
    time2show = 5; % seconds
    pnts2show = round(time2show*srate);
    trials2show = 1;

    % if the data is epoched
    if epoched

        % if the epochs are shorter than the desired time to show
        if EEG.xmax-EEG.xmin < time2show

            % sort epochs by variance
            v = squeeze(var(EEG.icaact(ic,:,:),[],2));
            [~,vsind] = sort(v);

            % find the median variance epoch and select it and the epochs
            % about it
            med = round(length(vsind)/2);
            chunks2show = vsind(med-floor((trials2show-1)/2):med+ceil((trials2show-1)/2));
            dat = squeeze(EEG.icaact(ic,:,chunks2show));
            time = EEG.times/1000;

        % if the epochs are longer than the desired time to show
        else
            % separate each epoch into non-overlapping windows and stack them
            tinds = 1:pnts2show:pnts-pnts2show;
            num_ts = length(tinds)*EEG.trials;
            ts = zeros(num_ts,pnts2show);
            index = reshape(1:length(tinds)*pnts2show,pnts2show,length(tinds))';
            for it = 1:EEG.trials
                blockind = (it-1)*length(tinds) + (1:length(tinds));
                temp_ica = EEG.icaact(ic,:,it);
                ts(blockind,:) = temp_ica(index);
            end
            times = EEG.times(index)/1000;

            % sort epochs by variance
            v = var(ts,[],2);
            [~,vsind] = sort(v);

            % find the median variance epoch and select it and the epochs
            % about it
            med = round(length(vsind)/2);
            chunks2show = vsind(med-floor((trials2show-1)/2):med+ceil((trials2show-1)/2));
            dat = ts(chunks2show,:)';
            time = linspace(0,time2show,pnts2show);

            tind = randi(pnts-pnts2show);
            time = linspace(0,time2show,pnts2show);
            dat = squeeze(EEG.icaact(ic,tind:tind+pnts2show-1,randperm(EEG.trials,trials2show)));
        end

    % if the data is epoched
    else

        % if the data is shorter than the desired time to show (shouldn't be)
        if EEG.xmax-EEG.xmin < time2show
            time = EEG.times/1000;
            dat = EEG.icaact(ic,:)';

        % if the data is longer than the desired time to show
        else
            % separate data into non-overlapping windows
            tinds = 1:pnts2show:pnts-pnts2show;

            % sort by variance
            v = size(tinds);
            for it = 1:length(tinds)
                v(it) = var(EEG.icaact(ic,tinds(it) + (1:pnts2show)));
            end
            [~,vsind] = sort(v);

            % find the median variance window and select it and the windows
            % about it
            med = round(length(vsind)/2);
            chunks2show = vsind(med-floor((trials2show-1)/2):med+ceil((trials2show-1)/2));
            dat = [];
            for it = 1:trials2show
                dat = [dat EEG.icaact(ic,tinds(chunks2show(it)) + (1:pnts2show))'];
            end
            time = linspace(0,time2show,pnts2show);
        end
    end
end

