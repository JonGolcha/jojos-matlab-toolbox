% function [OUTEEG] = pop_prop_landscape(EEG, varargin)
function [OUTEEG] = pop_prop_landscape(EEG,varargin)
% 
% *** Description:
% Display a summary image ('landscape') of component properties: 
% - scalp topography
% - dipole plot with residual variance and dipole moment ratio
% - ERP image
% - component time series
% - activity power spectrum
% - IC number and Percent Data Variance Accounted For 
% 
% In addition, allows to manually label components using checkboxes. 
% Possible labels are: 
% - 'brain'
% - 'eye'
% - 'muscle'
% - 'heart'
% - 'line'
% - 'channel'
% - 'other'
% - 'not_sure'
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
% - 'icnum'    	--> IC index to plot first
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
%  *** History
% updated by Ramon Martinez-Cancino and Luca Pion-Tonachini (2015)
% modified, cleaned and UI added (IC scrolling and annotations) by 
% Jonas Chatel-Goldman (2016)


%% initialize
% assess inputs
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
try, g.icnum;       catch, g.icnum = 1;            end;
try, g.spec_opt; 	catch, g.spec_opt   = {'freqrange',[0 100]};    end;
try, g.ID;          catch, g.ID         = 'NO_ID';                  end;
try, g.UI_control;  catch, g.UI_control = 'interactive';            end;
try, g.reset_labels; catch, g.reset_labels = false;                 end;
try, g.long_trial;  catch, g.long_trial = 'off';                    end;
try, EEG.icaweights; catch, error('ICA must be computed first');    end;
if isempty(EEG.chanlocs)
    error('channels location must be specified (field EEG.chanlocs is empty)'); 
end
nbIC = size(EEG.icaweights,1);

if g.reset_labels || ~isfield(EEG,'icalabels') || ~iscell(EEG.icalabels)
    EEG.icalabels = cell(nbIC,1);
    % example: 
    % EEG.icalabels{1} = {'brain'};
    % EEG.icalabels{2} = {'other','not_sure'}
end

% initialize figure
basename = ['Component_' int2str(g.icnum) ];
BACKCOLOR = [.93 .96 1];
hdl.mainfig = figure('name', ['SetID_' num2str(g.ID) '_' basename],...
    'color', BACKCOLOR,...
    'numbertitle', 'off',...
    'visible', 'on',...
    'PaperPositionMode','auto',...
    'ToolBar', 'none',...
    'MenuBar','none', ...
    'CloseRequestFcn', @closereq_cb);


set(hdl.mainfig,'WindowStyle','docked');
% pos = get(hdl.mainfig,'position');
% set(hdl.mainfig,'Position', [pos(1) pos(2)-700+pos(4) 1200 700]); % old small window
delete(gca);

% store user information for subsequent callbacks
myhandles = guihandles(hdl.mainfig); 
myhandles.EEG = EEG;
myhandles.g = g;
guidata(hdl.mainfig,myhandles)     
    
% initialize UI controls
if ~strcmp(g.UI_control,'no_ui')
    if nbIC~=1
        hdl.sliderSelect  = uicontrol(hdl.mainfig,'style','slider','Min',1,'Max',nbIC,'SliderStep',[1/(nbIC-1),5/(nbIC-1)],...
                                'units','normalized','Position',[.045 .03 .25 .05], 'BackgroundColor','w', ...
                                'Value',g.icnum,'HandleVisibility','off','Callback',@(hObject,callbackdata)slider_cb());  
    end
    % take key arrows as input to scroll components THIS DOESN'T WORK (TO FIX)
	% set(hdl.mainfig, 'WindowKeyPressFcn', @(hObject,evt)readkey_cb(evt.Key));
    % set(hdl.mainfig, 'WindowKeyPressFcn', @readkey_cb);
end

% Call display function
disp_landscape(g.icnum);

% Pause until figure is closed (mandatory to output EEG struct!)
if ~strcmp(g.UI_control,'no_ui')
    waitfor(hdl.mainfig);   
end


%% Display component properties landscape
function disp_landscape(current_IC) % nested function
    disp('***********');
    disp(['DISPLAYING COMPONENTS PROPERTIES: IC' int2str(g.icnum)]);

    % erase previous plots
	clf(hdl.mainfig);
    
    % show ui control for manual labelling
    if strcmp(g.UI_control,'interactive') 
        disp_uicontrol_interactive(hdl,current_IC,BACKCOLOR);
    end
    
    % recall EEG structure
    gui_handle = guidata(hdl.mainfig);
  	EEG = gui_handle.EEG;
    
    % plot time series
    if isempty(EEG.icaact)
        try
            EEG.icaact = eeg_getdatact(EEG, 'component', 1:size(EEG.icaweights,1)); 
        catch
            error('Could not obtain IC time courses. Was ICA already computed?');
        end
    end
    [dat, time] = choose_ts_seg(EEG,current_IC);
    datax = axes('position',[0.3712 0.75 0.5641 0.2002],'units','normalized');
    plot(time,dat'); 
    axis;
    axis tight
    ylim(datax, [-max(abs(dat)) max(abs(dat))])
    axis on;
    box on;
    grid on;
    title('Component Time Series','FontSize',10);
    xlabel(datax,'Time (s)','fontsize', 10);
    ylabel(datax,'\muV','FontSize',10);

    % plot scalp map
    axes('position',[0.0143 0.65 0.3121 0.3267],'units','normalized');
    scalpmap_norm = EEG.icawinv(:,current_IC)/sqrt(sum(EEG.icawinv(:,current_IC).^2));
    [htopo,Zi,plotrad] = topoplot( scalpmap_norm, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
        'shading', 'interp', 'numcontour', 0,'electrodes','on'); axis square;

    % add text information: icnum and pvaf
    try
        pvaf = num2str(EEG.pvaf(current_IC), '%4.2f');
    catch
        pvaf = 'N/A';
    end
    h = title({sprintf('IC %d of %d',current_IC,size(EEG.icawinv,2)); '{\bfData Var. Accounted For}:'; [pvaf '%'] }, ...
              'fontsize', 10,'fontweight','normal','Units','Normalized');
    set(h,'position',get(h,'position')+[-0.55 -1.13 0]);

    % add text information: ID
    text(0,1,[num2str(g.ID,'%06d') '_' num2str(current_IC,'%03d')],'Units', 'Normalized','Fontsize', 12,'Interpreter','none');


    % plot erpimage
    axes('position',[0.0643 0.2 0.2421 0.3850],'units','normalized');
    % eeglab_options;
    if EEG.trials > 1 % epoched data
        axis off
        EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
        if EEG.trials < 6
            ei_smooth = 1;
        else
            ei_smooth = 1;
        end
        icaacttmp = eeg_getdatact(EEG, 'component', current_IC) * norm(EEG.icawinv(:,current_IC));
        offset = nan_mean(icaacttmp(:));
        era    = nan_mean(squeeze(icaacttmp)')-offset;
        era_limits=get_era_limits(era);
        [~,~,~,~,axhndls] = erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
            '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','erp_vltg_ticks',era_limits); % strange automatic threshold for epochs > 10s with this option: 'limits',[g.epoch_lims(1) g.epoch_lims(2)]
        title('Epoched Data', 'fontsize', 10 );
           
    else % continuous data
        EI_TITLE = 'Continuous Data';
        ERPIMAGELINES = 200; % show 200-line erpimage
        while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
            ERPIMAGELINES = 0.9 * ERPIMAGELINES;
        end
        ERPIMAGELINES = round(ERPIMAGELINES);
        if ERPIMAGELINES > 2   % give up if data too small
            if ERPIMAGELINES < 10
                ei_smooth = 1;
            else
                ei_smooth = 1;
            end
            erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
            erpimageframestot = erpimageframes*ERPIMAGELINES;
            eegtimes = linspace(0, EEG.srate/1000, erpimageframes);
            icaacttmp = eeg_getdatact(EEG, 'component', current_IC) * norm(EEG.icawinv(:,current_IC));
            offset = nan_mean(icaacttmp(:));
            [~,~,~,~,axhndls] = erpimage(reshape(icaacttmp(:,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset,ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel','');
            ylabel(axhndls{1}, 'Data');
            xlabel(axhndls{1}, 'Time (s)');
            set(axhndls{1},'XTick',0:.05:25);
            if verLessThan('matlab', '8.3')
                set(axhndls{1},'XTickLabel', sprintf('%.1f|',[0:.2*erpimageframes/EEG.srate:erpimageframes/EEG.srate]));
            else
                set(axhndls{1},'XTickLabel', sprintf('%.1f\n',[0:.2*erpimageframes/EEG.srate:erpimageframes/EEG.srate]));
            end
        else
            axis off;
            text(0.1, 0.3, [ 'No erpimage plotted' 10 'for too small continuous data']);
        end;
    end;
    set(axhndls{2},'position', get(axhndls{2},'position') - [0.00 0 0.02 0]);
    lab = text(1.3,.85,'RMS \muVolts per channel');
    set(lab, 'rotation', -90)
    drawnow;

    % plot spectrum
    try
        hfreq1 = axes('position', [0.6567 0.2 0.2785 0.2], 'units', 'normalized');
        spectopo( EEG.icaact(current_IC,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,current_IC), 'percent', 50, g.spec_opt{:} );
        axis on tight;
        box on;
        grid on;
        xlim(hfreq1, [3 80]);
        xlabel('Frequency (Hz)', 'fontsize', 10)
        hfreq1.YLabel.Position = hfreq1.YLabel.Position + [5 25 0];
        hfreq1.YLabel.FontSize = 10;
        set(hfreq1, 'Xtick', [3 10:10:80],'fontsize', 8)

        hfreq2 = copyobj(hfreq1,gcf);
        hfreq2.Position = [0.6567 0.45 0.2785 0.2];
        hfreq2.XLabel = [];
        hfreq2.YLabel = [];
        xind = hfreq2.Children.XData>=3 & hfreq2.Children.XData<=40;
        xlim(hfreq2, [3 30]);
        title(hfreq2,'Activity Power Spectrum','units','normalized', 'fontsize', 10);

    catch
        axis off;
        lasterror
        text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
    end;
    drawnow;


    % dipplot
    % 1 dipole
    if isfield(EEG,'dipfit') && ~isempty(EEG.dipfit)
        try
            rv1 = num2str(EEG.dipfit.model(current_IC).rv*100, '%.1f');
        catch
            rv1 = 'N/A';
        end
        temp = axes('position', [0.4 0.2 0.2 0.15*3], 'units', 'normalized');
        imagesc(0)
        colormap(temp,[0 0 0])
        axis(temp, 'off')
        pos1 = [0.4 0.2 0.1 0.1557]; % [0.39 0.1109 0.1 0.1557]
        % axial
        % dipfitdefs
        ax(1) = axes('position', pos1, 'units', 'normalized');
        axis equal off
        dipplot(EEG.dipfit.model(current_IC), 'meshdata', EEG.dipfit.hdmfile, 'mri', EEG.dipfit.mrifile, ...
            'normlen', 'on', 'coordformat', 'MNI', 'axistight', 'on', 'gui', 'off', 'view', [0 0 1], 'pointout', 'on'); % bug in dipplot: modif JCG line 1161
        temp = get(ax(1),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        % coronal
        ax(2) = axes('position', [pos1(1) pos1(2)+.15 pos1(3) pos1(4)], 'units', 'normalized');
        axis equal off
        copyobj(allchild(ax(1)),ax(2));
        view([0 -1 0])
        axis equal off
        temp = get(ax(2),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        % sagital
        ax(3) = axes('position', [pos1(1) pos1(2)+.3 pos1(3) pos1(4)] , 'units', 'normalized');
        title(ax(3),'1 dipole')
        axis equal off
        copyobj(allchild(ax(1)),ax(3));
        view([1 0 0])
        axis equal off
        temp = get(ax(3),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        set(gcf, 'color', BACKCOLOR)
        drawnow;
    else
        % message no dipfit
        annotation('textbox',[.45,.45,.2,.05],'String','No dipole to plot!','LineStyle','none','FontSize',10)
%         text(0.2, 0.2, ['No dipole to plot!']);
    end
    % 2 dipole
    if isfield(EEG,'dipfit2') && ~isempty(EEG.dipfit2)
        try
            rv2 = num2str(EEG.dipfit2.model(current_IC).rv*100, '%.1f');
        catch
            rv2 = 'N/A';
        end
        % axial
        pos2 = [0.5 0.2 0.1 0.1557]; 
        ax(4) = axes('position', pos2, 'units', 'normalized');
        axis equal off
        dipplot(EEG.dipfit2.model(current_IC),'meshdata', EEG.dipfit.hdmfile, 'mri', EEG.dipfit.mrifile, ...
            'normlen', 'on', 'coordformat', 'MNI', 'axistight', 'on', 'gui', 'off', 'view', [0 0 1], 'pointout', 'on'); % bug in dipplot: modif JCG line 1161
        temp = get(ax(4),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        set(temp(6),'markersize',15,'color', 'm')
        set(temp(7),'linewidth',2,'color', 'm')
        % coronal
        ax(5) = axes('position', [pos2(1) pos2(2)+.15 pos2(3) pos2(4)], 'units', 'normalized');
        axis equal off
        copyobj(allchild(ax(4)),ax(5));
        view([0 -1 0])
        axis equal off
        temp = get(ax(5),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        set(temp(6),'markersize',15,'color', 'm')
        set(temp(7),'linewidth',2,'color', 'm')
        % sagital
        ax(6) = axes('position', [pos2(1) pos2(2)+.3 pos2(3) pos2(4)], 'units', 'normalized');
        axis equal off
        copyobj(allchild(ax(4)),ax(6));
        view([1 0 0])
        axis equal off
        temp = get(ax(6),'children');
        set(temp(4),'markersize',15)
        set(temp(5),'linewidth',2)
        set(temp(6),'markersize',15,'color', 'm')
        set(temp(7),'linewidth',2,'color', 'm')
        title(ax(6),'2 dipoles')
        hh = gcf;
        try
            hh.CurrentAxes = ax(1);
        end
        set(gcf, 'color', BACKCOLOR)
        % dipole text
        dmr = norm(EEG.dipfit2.model(current_IC).momxyz(1,:))/norm(EEG.dipfit2.model(current_IC).momxyz(2,:));
        if dmr<1
            dmr = 1/dmr; 
        end
        text(-70,-150,['RV: ' rv1 '%'])
        text(200,-163,{['RV: ' rv2 '%'];['DMR:' num2str(dmr,'%.1f')]})
        drawnow;
    else
        % message no dipfit
    end
end

% --- Close request (nested) function, executes when closing the figure window
function closereq_cb(hObject,eventdata) 
    gui_handle = guidata(gcf);
    OUTEEG = gui_handle.EEG;
    delete(hdl.mainfig); % close GUI
end

% --- Slider selection (nested) function, executes when scrolling
function slider_cb()
    current_val = round(get(gcbo,'Value'));
    gui_handle = guidata(gcf);
    if current_val ~= gui_handle.g.icnum
        gui_handle.g.icnum = current_val;
        guidata(gcf,gui_handle);
        set(gcbo,'Value',gui_handle.g.icnum);
        disp_landscape(gui_handle.g.icnum);
    end
end

% THIS DOESN'T WORK (TO FIX)
% --- Executes on arrow button press
% function readkey_cb(hObject,eventdata) 
%     if strcmp(eventdata.Key, 'rightarrow')==1
%         set(hdl.sliderSelect,'Value',round(get(gcbo,'Value'))+1);
%     elseif strcmp(eventdata.Key, 'leftarrow')==1
%         set(hdl.sliderSelect,'Value',round(get(gcbo,'Value'))-1);
%     end
%     slider_cb();
% end

end


% --- show ui control for manual labelling
function disp_uicontrol_interactive(hdl,current_IC,BACKCOLOR)
    % setup checkboxes for manual artifact annotation
    pos_checkBox = [.44 .035 .05 .05]; % position of first left checkbox
    space = .07; % horizontal space between checkboxes
	gui_handle = guidata(gcf);
  	EEG = gui_handle.EEG;
    % annotation('textbox',[pos_checkBox(1)-.1, pos_checkBox(2)+.03,.05,.02],'String','Labelling:','LineStyle','none','FontSize',12)
    
    % checkbox brain
    hdl.checkBox_brain = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'brain')),'Callback',@(hObject,callbackdata)checkBox_cb('brain'));                         
    annotation('textbox',[pos_checkBox(1)-.01,pos_checkBox(2)-.01,.05,.02],'String','brain','LineStyle','none')
    % checkbox eye
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_eye = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'eye')),'Callback',@(hObject,callbackdata)checkBox_cb('eye'));                         
    annotation('textbox',[pos_checkBox(1)-.01,pos_checkBox(2)-.01,.05,.02],'String','eye','LineStyle','none')
    
    % checkbox muscle
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_muscle = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'muscle')),'Callback',@(hObject,callbackdata)checkBox_cb('muscle'));                         
    annotation('textbox',[pos_checkBox(1)-.015,pos_checkBox(2)-.01,.05,.02],'String','muscle','LineStyle','none')
    
    % checkbox heart
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_heart = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'heart')),'Callback',@(hObject,callbackdata)checkBox_cb('heart'));                         
    annotation('textbox',[pos_checkBox(1)-.01,pos_checkBox(2)-.01,.05,.02],'String','heart','LineStyle','none')
    
    % checkbox line noise
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_line = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'line')),'Callback',@(hObject,callbackdata)checkBox_cb('line'));                         
    annotation('textbox',[pos_checkBox(1)-.01,pos_checkBox(2)-.01,.1,.02],'String','line','LineStyle','none')
    
    % checkbox channel noise
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_channel = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'channel')),'Callback',@(hObject,callbackdata)checkBox_cb('channel'));                         
    annotation('textbox',[pos_checkBox(1)-.02,pos_checkBox(2)-.01,.1,.02],'String','channel','LineStyle','none')
    
    % checkbox other
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_other = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'other')),'Callback',@(hObject,callbackdata)checkBox_cb('other'));                         
    annotation('textbox',[pos_checkBox(1)-.015,pos_checkBox(2)-.01,.05,.02],'String','other','LineStyle','none')
    
    % checkbox not sure
    pos_checkBox(1) = pos_checkBox(1)+space;
    hdl.checkBox_not_sure = uicontrol(hdl.mainfig,'style','checkbox','Min',0,'Max',1,...
                                'units','normalized','Position',pos_checkBox, 'BackgroundColor',BACKCOLOR, ...
                                'Value',sum(strcmp(EEG.icalabels{current_IC},'not_sure')),'Callback',@(hObject,callbackdata)checkBox_cb('not_sure'));                         
    annotation('textbox',[pos_checkBox(1)-.025,pos_checkBox(2)-.01,.1,.02],'String','(not sure)','LineStyle','none')
end


% --- Executes on button press in checkBoxes
function checkBox_cb(label)
    gui_handle = guidata(gcf);
	if fix(get(gcbo,'Value')) == 0
      	if sum((strcmp(gui_handle.EEG.icalabels{gui_handle.g.icnum},label)))>0
            gui_handle.EEG.icalabels{gui_handle.g.icnum}((strcmp(gui_handle.EEG.icalabels{gui_handle.g.icnum},label)))=[];
        end
    else
       if sum((strcmp(gui_handle.EEG.icalabels{gui_handle.g.icnum},label)))==0
           gui_handle.EEG.icalabels{gui_handle.g.icnum}(end+1) = {label};
       end
    end
    guidata(gcf,gui_handle);
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









% Code for solving bug in dipplot() (add this at line 1161)
%    if indx>size(dat.imgs,1)
%        indx = size(dat.imgs,1);
%    end
%    if indx<1
%        indx = 1;
%    end
%    if indy>size(dat.imgs,2)
%        indy = size(dat.imgs,2);
%    end
%    if indy<1
%        indy = 1;
%    end
%    if indz>size(dat.imgs,3)
%        indz = size(dat.imgs,3);
%    end
%    if indz<1
%        indz = 1;
%    end


