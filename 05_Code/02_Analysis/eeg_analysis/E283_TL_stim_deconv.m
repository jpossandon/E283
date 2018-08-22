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
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'&origstart','<1500';'&latposStim','>-1000'},...
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
    [EEG,winrej]            = getDataDeconv(cfg_eeg,epochevents,200); 
    mirindx                 = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
             
    if any(strfind(p.analysisname,'CI'))
        LstimTimes = epochevents.latency(strcmp(epochevents.type,'stim') & epochevents.side == -1);
        for tt  = 1:length(LstimTimes)
            EEGst = find(EEG.times>LstimTimes(tt),1);
            mirsamples = EEGst-floor(EEG.srate*.4):EEGst+ceil(EEG.srate*1.6);
            EEG.data(:,mirsamples) = EEG.data(mirindx,mirsamples);
            ETst  = find(epochevents.latency>LstimTimes(tt)-400 & epochevents.latency<LstimTimes(tt)+900 & ...
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
    elseif strcmp(model,'IM_STcsi_Sdxyendampang')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {sprintf('y~pyend+pxend+spl(amp,10)+circspl(angle,10,%3.1f,%3.1f)',nanmin(epochevents.angle),nanmax(epochevents.angle)),...
                                 'y~1','y~cross*inst*side'};
    elseif strcmp(model,'IM_STcsi_Sdxydendxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)+spl(pxend,10)+spl(pyend,10)',...
                                'y~1','y~cross*inst*side'};
    elseif strcmp(model,'IM_STcsi_Sdxy')
        cfgDesign.eventtypes= {'sac','image','stim'};
        cfgDesign.formula   = {'y~1+spl(pxdiff,10)+spl(pydiff,10)',...
                                'y~1','y~cross*inst*side'};
    end
    EEG                 = uf_designmat(EEG,cfgDesign);
    cfgTexp             = [];
    cfgTexp.timelimits  = [-1,1];tic
    EEG                 = uf_timeexpandDesignmat(EEG,cfgTexp);toc
    EEG                 = uf_continuousArtifactExclude(EEG,struct('winrej',winrej));
    EEG                 = uf_glmfit(EEG);
 
    unfold              = uf_condense(EEG);

    mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
    save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
    clear unfold
end 
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
bslcor    = [];
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
    elseif  strcmp(model,'IM_STcsi_Sdxy') 
       % pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
       % pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         pxdM        = 10;   pxd_pred    = [-pxdM:.5:pxdM];
         pydM        = 10;   pyd_pred    = [-pydM:.5:pydM];
        
        ufpredict               = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
        subjparvalues{tk,1} = [unfold.unfold.splines{1}.paramValues(~isnan(unfold.unfold.splines{1}.paramValues))];
        subjparvalues{tk,2} = [unfold.unfold.splines{2}.paramValues(~isnan(unfold.unfold.splines{2}.paramValues))];
    end
        % remove splines from coefficients estimates
     if strcmp(model,'IM_STcsi_Sdxydendxy') ||  strcmp(model,'IM_STcsi_Sdxy')   
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

% run 2nd level stat and save
E283_run2nd_save

%%
%%
% plot betas
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr']);
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')]);
end
if any(strfind(p.analysisname,'mirr'))
    plotBetasTopos(cfg_eeg,result,pathfig,[-.12  .6 .015],[],1)
else
    plotBetasTopos(cfg_eeg,result,pathfig,[-.12  .6 .03],[],0)
end

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SELECTED CHANNELS WITH SUBJECT VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.analysisname  = 'deconvTF';
% cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
  
cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop

%chnstoPlot      = {{'Cz','CPz'},{'Pz','POz'},{'O1','Oz','O2'},...
 %   {'P5','P7','PO7'},{'P6','P8','PO8'},{'C3','CP3'}};
chnstoPlot      = {{'FCz','Cz','CPz'},{'O1','Oz','O2'},{'P1','Pz','P2','PO3','POz','PO4'}};
chnstoPlot      = {{'C3','C5'},{'C4','C6'},{'P5','P7','PO7'},{'P6','P8','PO8'}};
Bnames        = {'stim_side'};% Bnames        = {'fixpreT0_Intercept','fixpre_3_Intercept','fixpre_orderPreT_2','fixpre_orderPreT_3','fixpre_orderPreT_4','fixpre_orderPreT_5'}
axLim           = [-1 .1 -10 10];
plotBetasChannels(cfg_eeg,result,chnstoPlot,Bnames,pathfig,axLim,'preTarget')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(bslcor)
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],'splines');
else
    pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],'splines');
end
mkdir(pathfig)

if any(strfind(p.analysisname,'mirr'))
   plotinterval = [-.12  .6 .015];
   half =1;
else
    plotinterval = [-.12  .6 .03];
    half = 0;
end
setAbsoluteFigureSize
splnames = cellfun(@(x) x{1},regexp(resultSplines.coeffs,'^([.*_]?.*_)','match'), 'UniformOutput', false)

for splType=unique(splnames)'
    splb        = strmatch(splType,splnames);
    betas.dof   = 1;
    betas.n     = size(resultSplines.beta,4);
    betas.time  = resultSplines.times;
    stat.time   = resultSplines.times;
    %     collim      =[-6 6]; 
    for splVal = splb'
        betas.avg   = squeeze(mean(resultSplines.beta(:,splVal,:,:),4));
        collim      =[-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
    
        for pint = 1:size(plotinterval,1)
            fh       = topomitlines(cfg_eeg,stat,betas,plotinterval(pint,:),collim,half);
            figsize  = [17.6 17.6*fh.Position(4)/fh.Position(3)];
               doimage(gcf,pathfig,'pdf',[resultSplines.coeffs{splVal} '_' strjoin('_',{num2str(plotinterval(pint,1)),num2str(plotinterval(pint,2))})],figsize,1)
        end
    end
end

% difference at same vector
for splType=unique(splnames)'
    splb        = strmatch(splType,splnames);
    betas.dof   = 1;
    betas.n     = size(resultSplines.beta,4);
    betas.time  = resultSplines.times;
    stat.time   = resultSplines.times;
    %     collim      =[-6 6]; 
    for splVal = 1:floor(length(splb)/2)
        betas.avg   = squeeze(mean(resultSplines.beta(:,splb(end-splVal+1),:,:)-resultSplines.beta(:,splb(splVal),:,:),4));
        collim      =[-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
    
        for pint = 1:size(plotinterval,1)
            fh       = topomitlines(cfg_eeg,stat,betas,plotinterval(pint,:),collim,half);
            figsize  = [17.6 17.6*fh.Position(4)/fh.Position(3)];
              doimage(gcf,pathfig,'pdf',[resultSplines.coeffs{splb(end-splVal+1)} '_minus_' resultSplines.coeffs{splb(splVal)} '_' strjoin('_',{num2str(plotinterval(pint,1)),num2str(plotinterval(pint,2))})],figsize,1)
        end
    end
end
%%
% test 2d splines
%     setAbsoluteFigureSize
% %
% for t =.11
% timetoplot = [t  t+.02 .02];
% 
%     betas.dof   = 1;
%     betas.n     = 1;
%     betas.avg   = ufpredict.beta(:,:,2:122)+ repmat(ufpredict.beta(:,:,1),1,1,121);
%     betas.time  = ufpredict.times;
%     indxT       = find(betas.time>=timetoplot(1),1);
% 
%     betas.param = ufpredict.param(2:122);
%     stat.time = betas.time;
%     collim = [-8 8]
% %     collim      =[-3*std(reshape(squeeze(betas.avg(:,indxT,:)),1,size(betas.avg,1)*size(betas.avg,3))) 3*std(reshape(squeeze(betas.avg(:,indxT,:)),1,size(betas.avg,1)*size(betas.avg,3)))]; 
%  fh = topogrid(cfg,stat,betas,timetoplot,collim)  
%  fh.Name = num2str(t)
% end

