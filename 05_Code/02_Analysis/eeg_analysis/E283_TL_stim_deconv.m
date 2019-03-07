% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'TL_dc_Stim';
%model           = 'IM_STcsi';
% % model           = 'IM_STcsi_Sdxyendampang';
model           = 'IM_STcsi_Sdxy';
if ismac 
    run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')        
else
    run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')   
end   
%%
% subject configuration and data
%{ 
for tk = p.subj
    tk
%  for tk = p.subj;
% tk = str2num(getenv('SGE_TASK_ID'));
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jpo/trabajo/E283/');
    end
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname',p.analysisname);    % single experiment/session parameters 
   
    % get relevant epochevents
     load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'&origstart','<1000';'&latposStim','>-1000'},...
                                [800 100],{-1,2,'origstart','>0'}); 
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.type        = cell(1,length(events.start));
    epochevents.type(events.type==1) = repmat({'fix'},1,sum(events.type==1));
    epochevents.type(events.type==2) = repmat({'sac'},1,sum(events.type==2));

    epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    events.angle(events.angle<0) = 360+events.angle(events.angle<0);           % transform angles to 0-360
    epochevents.angle       = events.angle;
    epochevents.amp         = events.amp;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini; 
    epochevents.pxdiff(2:2:end) = epochevents.pxdiff(1:2:end);         % fixation vector is the same sa previous saccade
    epochevents.pydiff(2:2:end) = epochevents.pydiff(1:2:end);
    epochevents.side        = nan(1,length(events.start));    
    epochevents.cross       = nan(1,length(events.start));    
    epochevents.inst        = nan(1,length(events.start)); 
    
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 900]);
    events                      = struct_elim(events,find(~ismember(events.value,[1:6,9,10,13,14,96])),2,0);
    epochevents.latency         = [epochevents.latency,events.time];
    ETttype                     = cell(1,length(events.value));
    ETttype(events.value==96)   = repmat({'image'},1,sum(events.value==96));
    ETttype(events.value<96)    = repmat({'stim'},1,sum(events.value<96));
    epochevents.type            = [epochevents.type,ETttype];
    ETside                      = nan(1,length(events.value));
    ETside(ismember(events.value,[1 3 5 9 13])) = -1;    % left
    ETside(ismember(events.value,[2 4 6 10 14])) = 1;
    epochevents.side            = [epochevents.side,ETside];
    ETcross                     = nan(1,length(events.value));
    ETcross(ismember(events.value,[1 2 9 10])) = -1;    % uncross
    ETcross(ismember(events.value,[3:6 13 14])) = 1;
    epochevents.cross           = [epochevents.cross,ETcross];
    ETinst                      = nan(1,length(events.value));
    ETinst(ismember(events.value,[1 2 5 6])) = 1;  % instructive
    ETinst(ismember(events.value,[9 10 13 14])) = -1;
    epochevents.inst            = [epochevents.inst,ETinst];
    
    epochevents.pxini           = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini           = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend           = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend           = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff          = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff          = [epochevents.pydiff,nan(1,length(events.value))];  
    epochevents.orderposStim    = [epochevents.orderposStim,nan(1,length(events.value))];  
    epochevents.angle           = [epochevents.angle,nan(1,length(events.value))];
    epochevents.amp             = [epochevents.amp,nan(1,length(events.value))];
     
    % getting the data in EEGlab format
    [EEG,winrej]            = getDataDeconv(cfg_eeg,epochevents,200,1); 
    mirindx                 = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
             
    if any(strfind(p.analysisname,'CI'))
        LstimTimes = epochevents.latency(find(strcmp(epochevents.type,'stim') & epochevents.side == -1)-1);
        for tt  = 1:length(LstimTimes)
            EEGst = find(EEG.times>LstimTimes(tt),1);
            mirsamples = EEGst-floor(EEG.srate*.4):EEGst+ceil(EEG.srate*1.6);
            EEG.data(:,mirsamples) = EEG.data(mirindx,mirsamples);
            ETst  = find(epochevents.latency>LstimTimes(tt)-400 & epochevents.latency<LstimTimes(tt)+1000 & ...
                (strcmp(epochevents.type,'sac') | strcmp(epochevents.type,'fix')));
            if ~isempty(ETst)
                [epochevents.pxini(ETst) epochevents.pxend(ETst) epochevents.pxdiff(ETst)] = deal(-epochevents.pxini(ETst),-epochevents.pxend(ETst),-epochevents.pxdiff(ETst));
            end
        end
    end
    if any(strfind(p.analysisname,'mirr'))
        EEG.data = EEG.data-EEG.data(mirindx,:);   
    end
    % deconvolution design
    cfgDesign           = [];
%      cfgDesign.eventtypes= {'sac','image','stim'};
%      cfgDesign.formula   = {'y~pxdiff+pydiff+pxend+pyend','y~1','y~cross*inst*side'};
%      model               = 'Sxydxyend_IM_STsci';
    if strcmp(model,'IM_STcsi')
        cfgDesign.eventtypes= {'image','stim'};
        cfgDesign.formula   = {'y~1','y~cross*inst*side'};
      elseif strcmp(model,'IM_STcsi_Sdxydendxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)+spl(pxend,10)+spl(pyend,10)',...
                                'y~1','y~cross*inst*side'};
    elseif strcmp(model,'IM_STcsi_Sdxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)',...
                                'y~1','y~cross*inst*side'};
    elseif strcmp(model,'IM_STci_Sdxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)',...
                                'y~1','y~cross*inst'};
 
    end
    EEG                 = uf_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-.5,.8];tic
    EEG                 = uf_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
   % cfgfit.method       = 'glmnet';
   % cfgfit.glmnetalpha  = 0;
    %EEG                 = uf_glmfit(EEG,cfgfit);
   
    EEG                 = uf_glmfit(EEG);
    unfold              = uf_condense(EEG);

    mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
    save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    clear unfold
end 
%}
% ufpredict = uf_predictContinuous(unfold,'predictAt',{{'amp',[.1 .5 1 2.5 5 7.5 10]}});
% ufmarginal = uf_addmarginal(ufpredict)
% uf_plotParam(ufpredict,'channel',35,'plotParam',{'(Intercept)','amp'},'add_intercept',1,'baseline',[-.2 0])
%%
% %2nd level analysis
% clear
% E283_params                                 % basic experimental parameters               % 
 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', p.analysisname);
 end
    
stimB       = [];
stimBspl    = [];
ufpredict   = [];
bslcor    = [-.5 -.25];
subjparvalues = {};
for tk = p.subj
     cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
     load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    if strcmp(model,'IM_STcsi_Sdxydendxy')
%         ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'amp',[.1 .5 1 2.5 5 7.5 10]},{'angle',[0:45:315]}});
        pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
        pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
        pxeM        = 12;   pxe_pred    = [-fliplr((2.^[-1:4]/2^4)*pxeM) 0 (2.^[-1:4]/2^4)*pxeM];
        pyeM        = 9;    pye_pred    = [-fliplr((2.^[-1:4]/2^4)*pyeM) 0 (2.^[-1:4]/2^4)*pyeM];
        
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred},{'pxend',pxe_pred},{'pyend',pye_pred}});
    elseif  strcmp(model,'IM_STcsi_Sdxy')  ||  strcmp(model,'IM_STci_Sdxy')      
%         pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
%         pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         pxdM        = 12;   pxd_pred    = [-pxdM:1.5:pxdM];
         pydM        = 12;   pyd_pred    = [-pydM:1.5:pydM];
        
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
        subjparvalues{tk,1} = [unfold.unfold.splines{1}.paramValues(~isnan(unfold.unfold.splines{1}.paramValues))];
        subjparvalues{tk,2} = [unfold.unfold.splines{2}.paramValues(~isnan(unfold.unfold.splines{2}.paramValues))];
    end
        % remove splines from coefficients estimates
     if any(strfind(model,'spl')) ||  any(strfind(model,'Sdxy'))   
        rmspl                   = strmatch('spline',{unfold.param.type});
        unfold.beta(:,:,rmspl)  = [];
        unfold.param(rmspl)     = [];
        
        rmnospl                   = setdiff(1:length(ufpredict.param),strmatch('spline',{ufpredict.param.type}));
        ufpredict.beta(:,:,rmnospl) = [];
        ufpredict.param(rmnospl)    = [];
    end
    if ~isempty(bslcor) 
        stimB       = cat(4,stimB,permute(unfold.beta-...
            repmat(mean(unfold.beta(:,find(unfold.times>bslcor(1) & unfold.times<bslcor(2)),:),2),1,size(unfold.beta,2),1),[1,3,2]));
        if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta-...
                repmat(mean(ufpredict.beta(:,find(ufpredict.times>bslcor(1) & ufpredict.times<bslcor(2)),:),2),1,size(ufpredict.beta,2),1),[1,3,2]));
        end
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
        if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta(:,:,:),[1,3,2]));
        end
    end

end
%
% run 2nd level stat and save
stattype = 'signpermT';
mc = 'cluster';
E283_run2nd_save

%%
%%
% plot betas
% bslcor    = [];
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
if any(strfind(p.analysisname,'mirr'))
    plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.12  .6 .015],[],1)
else
    plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.2  .7 .05],[],0)
end

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
% cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
  
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop

chnstoPlot      = {{'FCz','Cz','CPz'}}
linecolors      = [0 0 0];
Bnames         = {'stim_cross'}; Blabels         = {''};
axLim           = [-.2 .8 -1 1];
plotBetasChannels(cfg_eeg,'boot_se',result,chnstoPlot,Bnames,Blabels,linecolors,pathfig,axLim,'stim_cross')

chnstoPlot      = {{'C3','C5'},{'C4','C6'},{'FCz','Cz','CPz'}};
% chnstoPlot      = {{'FCz','Cz','CPz'}};
%Bnames         = {'image_2_Intercept','stim_3_Intercept','sac_Intercept'} 
%Blabels        = {'image','touch','saccade'} 
chnstoPlot      = {{'FCz','Cz','CPz'}}
Bnames         = {'stim_cross','stim_crossXXinst','stim_inst'} 
Blabels        = {'cross','crossXinst','inst'} 
axLim           = [-.2 .7 -1.5 3];
% plotBetasChannels(cfg_eeg,result,chnstoPlot,Bnames,pathfig,axLim,'preTarget')
lineColors      = cbrewer('qual','Set1',9)
plotBetasChannels(cfg_eeg,'mean',result,chnstoPlot,Bnames,Blabels,[1 1 1],lineColors,pathfig,axLim,'stimLock_intercepts')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bslcor    = [-.5 -.15];
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],'splines');
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],'splines');
end

if any(strfind(p.analysisname,'mirr'))
    plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3 .7 .025],[],1,1)
else
    plotSplinesTopos(cfg_eeg,resultSplines,'mean',pathfig,[-.3  .7 .05],[],0,1)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chnstoPlot      = {{'P5','P7','PO7'},{'P6','P8','PO8'}};
% Bnames        = {'diff:sac_pxdiff_    -10:sac_pxdiff_     10',...
%                  'diff:sac_pxdiff_     -5:sac_pxdiff_      5',...
%                  'diff:sac_pxdiff_   -2.5:sac_pxdiff_    2.5',...
%                  'diff:sac_pxdiff_  -1.25:sac_pxdiff_   1.25',...
%                  'diff:sac_pxdiff_ -0.625:sac_pxdiff_  0.625',...
%                  'diff:sac_pxdiff_-0.3125:sac_pydiff_ 0.3125'};

Bnames = {resultSplines.coeffs{strmatch('sac_pxdiff',resultSplines.coeffs)}};
Blabels = {[-12:1.5:12]};
axLim           = [-.3 .7 -4 8];
lineColors      = cbrewer('qual','Set1',9)
lineColors      = cbrewer('div','RdYlGn',length(Bnames));
% plotBetasChannels(cfg_eeg,resultSplines,chnstoPlot,Bnames,pathfig,axLim,'preTarget')
plotBetasChannels(cfg_eeg,'mean',resultSplines,chnstoPlot,Bnames,Blabels,lineColors,pathfig,axLim,'spline')

%%
% checking correlaiton between cross coefficient and performance

crossB = squeeze(result.B(:,strcmp(result.coeffs,'stim_cross'),:,:));