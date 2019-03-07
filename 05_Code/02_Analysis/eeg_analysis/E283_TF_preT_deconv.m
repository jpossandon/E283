% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'deconvTFpreT';
%model           = 'IM_STcsi';
 model = 'S_orderPrecat_xdifspl_IM_STsc_eff';
%
% subject configuration and data
 
if ismac 
run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')        
else
run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')   
end    
%%
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
    % check eyedad is the same as full data

    nmovs = 7;
    eyedata.events.nextToTarg = any([eyedata.events.nextToTargH;eyedata.events.nextToTargV;eyedata.events.nextToTargD]);
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT',['<' num2str(nmovs)]},...
        [800 100],{-1,2,'origstart','>0';-2,1,'origstart','>0'});
    events.DToTarg = events.DToTarg./45; 
    events.onTarg(2:3:end)  = events.onTarg(1:3:end);
    events.DToTarg(2:3:end) = events.DToTarg(1:3:end);         
    events = struct_elim(events,[1:3:length(events.type)],2);
    
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
     epochevents.dur         = events.dur;
    epochevents.type        = cell(1,length(events.start));
     epochevents.trial      = events.trial;
    epochevents.orderPreT   =  events.orderPreT;
    epochevents.DToTarg     = events.DToTarg;
    if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff') || strcmp(model,'S_orderPrecat_evidence_xdifspl_IM_STsc_eff') || strcmp(model,'S_orderPrecat_evidencespl_xdifspl_IM_STsc_eff')  || strcmp(model,'S_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
        epochevents.type(events.type==2 & events.orderPreT==0 & events.onTarg) = repmat({'sacpreT0onTarg'},1,sum(events.type==2 & events.orderPreT==0 & events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT==0 & ~events.onTarg) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT==0 & ~events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT>0) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT>0));
    elseif strcmp(model,'S_orderPrecat_evidencepre1spl_xdif_IM_STsc_eff')
        epochevents.type(events.type==2 & events.orderPreT==0 & events.onTarg) = repmat({'sacpreT0onTarg'},1,sum(events.type==2 & events.orderPreT==0 & events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT==0 & ~events.onTarg) = repmat({'sacpreT0'},1,sum(events.type==2 & events.orderPreT==0 & ~events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT>0) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT>0));
   
    elseif strcmp(model,'S_orderPrecat_xdifspl_prefixdurspl_IM_STsc_eff')
        epochevents.type(events.type==2 & events.orderPreT==0 & events.onTarg) = repmat({'sacpreT0onTarg'},1,sum(events.type==2 & events.orderPreT==0 & events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT==0 & ~events.onTarg) = repmat({'sacpreT0'},1,sum(events.type==2 & events.orderPreT==0 & ~events.onTarg));
        epochevents.type(events.type==2 & events.orderPreT>0) = repmat({'sacpre'},1,sum(events.type==2 & events.orderPreT>0));
    end
    
    if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff') || strcmp(model,'F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
        epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
        epochevents.type(events.type==1 & events.orderPreT==1 & events.onTarg)    = repmat({'fixpreT1onTarg'},1,sum(events.type==1  & events.orderPreT==1 & events.onTarg));
        epochevents.type(events.type==1 & ((events.orderPreT==1 & ~events.onTarg) | events.orderPreT>1))   = repmat({'fixpre'},1,sum(events.type==1 & ((events.orderPreT==1 & ~events.onTarg) | events.orderPreT>1)));
        
    elseif strcmp(model,'F0_xdifspl_F1_xdifspl_DToTargspl_F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff') 
        epochevents.type(events.type==1 & events.orderPreT==0) = repmat({'fixpreT0'},1,sum(events.type==1 & events.orderPreT==0));
        epochevents.type(events.type==1 & events.orderPreT==1 & events.onTarg)    = repmat({'fixpreT1onTarg'},1,sum(events.type==1  & events.orderPreT==1 & events.onTarg));
        epochevents.type(events.type==1 & events.orderPreT>1)   = repmat({'fixpre'},1,sum(events.type==1 & events.orderPreT>1));
       
        epochevents.type(events.type==1 & events.orderPreT==1 & ~events.onTarg)    = repmat({'fixpreT1'},1,sum(events.type==1 & events.orderPreT==1 & ~events.onTarg));
        
    end
    
    
%     epochevents.orderposStim= events.orderposStim;
    epochevents.pxini       = (events.posinix-960)/45;
    epochevents.pyini       = (events.posiniy-540)/45;
    epochevents.pxend       = (events.posendx-960)/45;
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini;
    epochevents.pxdiff(2:2:end)  = epochevents.pxdiff(1:2:end);
    epochevents.pydiff(2:2:end)  = epochevents.pydiff(1:2:end);
    epochevents.amp       = events.amp;
    epochevents.side        = nan(1,length(events.start));
    epochevents.cross       = nan(1,length(events.start));
    epochevents.inst        = nan(1,length(events.start));
    
    if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff') || ...
            strcmp(model,'S_orderPrecat_xdifspl_prefixdurspl_IM_STsc_eff') ||  strcmp(model,'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc') || ...
            strcmp(model,'S_orderPrecat_evidence_xdifspl_IM_STsc_eff') || strcmp(model,'S_orderPrecat_evidencepre1spl_xdif_IM_STsc_eff') || ...
            strcmp(model,'S_orderPrecat_evidencespl_xdifspl_IM_STsc_eff') || strcmp(model,'S_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
        epochevents = struct_elim(epochevents,[2:2:length(epochevents.type)],2);
    end
    if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff') || strcmp(model,'F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff') || ...
        strcmp(model,'F0_xdifspl_F1_xdifspl_DToTargspl_F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff');
        epochevents = struct_elim(epochevents,[1:2:length(epochevents.type)],2);
    end
    
    if any(epochevents.orderPreT==nmovs)
        epochevents = struct_elim(epochevents,find(epochevents.orderPreT==nmovs),2);
    end
        
    
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 900]);
    events                  = struct_elim(events,find(~ismember(events.value,[1:6,9,10,13,14,96])),2,0);
    epochevents.latency     = [epochevents.latency,events.time];
    ETttype                 = cell(1,length(events.value));
    ETttype(events.value==96)   = repmat({'image'},1,sum(events.value==96));
    ETttype(events.value<96)    = repmat({'stim'},1,sum(events.value<96));
    epochevents.type            = [epochevents.type,ETttype];
    ETside                  = nan(1,length(events.value));
    ETside(ismember(events.value,[1 3 5 9 13])) = -1;    % left
    ETside(ismember(events.value,[2 4 6 10 14])) = 1;
    epochevents.side        = [epochevents.side,ETside];
    ETcross                 = nan(1,length(events.value));
    ETcross(ismember(events.value,[1 2 9 10])) = -1;    % uncross
    ETcross(ismember(events.value,[3:6 13 14])) = 1;
    epochevents.cross       = [epochevents.cross,ETcross];
    ETinst                 = nan(1,length(events.value));
    ETinst(ismember(events.value,[1 2 5 6])) = 1;  % instructive
    ETinst(ismember(events.value,[9 10 13 14])) = -1;
    epochevents.inst       = [epochevents.inst,ETinst];
    epochevents.trial      = [epochevents.trial,events.trial];
        
     epochevents.dur       = [epochevents.dur,nan(1,length(events.value))];
    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];
    epochevents.orderPreT       = [epochevents.orderPreT nan(1,length(events.value))];
    epochevents.amp             = [epochevents.amp,nan(1,length(events.value))];
    epochevents.DToTarg        = [epochevents.DToTarg,nan(1,length(events.value))];
    % getting the data in EEGlab format
    % without trials that have a prefixation on target
    if any(strfind(model,'S_'))
        epochevents = struct_elim(epochevents,find(ismember(epochevents.trial,epochevents.trial(strcmp(epochevents.type,'sacpreT0onTarg')))),2);
    elseif any(strfind(model,'F_'))
        epochevents = struct_elim(epochevents,find(ismember(epochevents.trial,epochevents.trial(strcmp(epochevents.type,'fixpreT1onTarg')))),2);
    end
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,250,1);
    mirindx         = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']);
    
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
    if strcmp(model,'F_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'F_orderPrecat_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1onTarg','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cross*inst*side'};
    end
    
    if strcmp(model,'F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1onTarg','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    
    if strcmp(model,'F_orderPrecat_dirtgt_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0T1onTarg','fixpreT1','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(dirtotarg)+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cross*inst*side'};
    end
    if strcmp(model,'F0_xdifspl_F1_xdifspl_DToTargspl_F_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'fixpreT0','fixpreT1','fixpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    
    if strcmp(model,'S_orderPrecat_xdifspl_IM_STsc') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'sacpre','image','stim'};
        cfgDesign.formula   = {'y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    if strcmp(model,'S_orderPrecat_evidence_xdifspl_IM_STsc_eff')
       cfgDesign.eventtypes = {'sacpreT0onTarg','sacpre','image','stim'};
       cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)+evidence+evidence2+evidence3+evidence4','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    if strcmp(model,'S_orderPrecat_evidencespl_xdifspl_IM_STsc_eff')
        cfgDesign.eventtypes = {'sacpreT0onTarg','sacpre','image','stim'};
        cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)+spl(evidence,5)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    if strcmp(model,'S_orderPrecat_evidencepre1spl_xdif_IM_STsc_eff')
       cfgDesign.eventtypes = {'sacpreT0onTarg','sacpreT0','sacpre','image','stim'};
       cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+spl(pxdiff,10)+spl(pydiff,10)+spl(evidence,6)+spl(evidence2,6)+spl(evidence3,6)','y ~1+spl(pxdiff,10)+spl(pydiff,10)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    if strcmp(model,'S_orderPrecat_DToTargspl_xdifspl_IM_STsc_eff')
       cfgDesign.eventtypes = {'sacpreT0onTarg','sacpre','image','stim'};
       cfgDesign.formula   = {'y ~1+spl(pxdiff,10)+spl(pydiff,10)','y ~1+cat(orderPreT)+spl(pxdiff,10)+spl(pydiff,10)+spl(DToTarg,10)','y~1','y~cat(cross)*cat(inst)*cat(side)'};
    end
    
    cfgDesign.codingschema = 'effects';
    %freqbands       = {'alfa','beta'};
    %bplim           = [9 15;16 26];
    
    freqbands       = {'alfa','beta'};
    bplim           = [9 15;16 25];
    filtPnts        = [368,254]; % check this
    trSt            = strmatch('S 96',{EEG.event.type});
    evLat           = [EEG.event.latency];
    lat_trSt        = round(evLat(trSt));
    bsl             = [-.450 0];
    winsize         = .3/(1/EEG.srate);
%     for tt = winsize/2+1:1:EEG.pnts-winsize/2-1
%         auxd = EEG.data(:,tt-winsize/2:tt+winsize/2-1);

    for fb = 1:length(freqbands)
      
      
        % this two are quite equivalent but holber path is faster
        EEGaux              = pop_eegfiltnew(EEG, bplim(fb,1), bplim(fb,2), [], 0, [], 0);
        EEGaux.data         = abs(hilbert(EEGaux.data')').^2; %'power'
          %foi = bplim(fb,1):bplim(fb,2)
         % EEGaux= EEG;
%         for ch = 1:EEG.nbchan
%             [spectrum_mtmconvol,ntaper,foi,toi] = ft_specest_mtmconvol(EEG.data(ch,:),EEG.times/1000,...
%                 'timeoi',EEG.times/1000,'timwin',3./foi ,'taper','hanning',...
%                 'pad',[],'freqoi',foi,'polyorder',0,'dimord','chan_time_freqtap');
%             EEGaux.data(ch,:) =  sum(abs(spectrum_mtmconvol).^2,3);
%         end
        bslsample           = floor(bsl.*EEGaux.srate);
        sampltobsl          = [-.4 1.6].*EEGaux.srate;
        for ll = 1:length(lat_trSt)
            trlbslpow(:,ll)       = mean(EEGaux.data(:,bslsample(1)+lat_trSt(ll):bslsample(2)+lat_trSt(ll)),2);
        end
        EEGaux.data = 10*log10(EEGaux.data./repmat(mean(trlbslpow,2),1,EEGaux.pnts));
    %        EEGaux.data(:,sampltobsl(1)+lat_trSt(ll):sampltobsl(2)+lat_trSt(ll)) = 10*log10(EEGaux.data(:,sampltobsl(1)+lat_trSt(ll):sampltobsl(2)+lat_trSt(ll))./...
     %           repmat(trlbslpow,1,diff(sampltobsl)+1)); 
    
        %EEGaux.data         = 10*log10(EEGaux.data./repmat(nanmean(EEGaux.data,2),1,size(EEGaux.data,2))); 
        if any(strfind(p.analysisname,'mirr'))
            EEGaux.data = EEGaux.data-EEGaux.data(mirindx,:); 
        end
        EEGaux              = uf_designmat(EEGaux,cfgDesign);
        cfgTexp             = [];
        cfgTexp.timelimits  = [-.5,.8];tic
        EEGaux              = uf_timeexpandDesignmat(EEGaux,cfgTexp);toc
        EEGaux              = uf_continuousArtifactExclude(EEGaux,struct('winrej',winrej));
        EEGaux              = uf_glmfit(EEGaux);
 
        unfold.(freqbands{fb})          = uf_condense(EEGaux);
        
    end
  for fb = 1:length(freqbands)
      unfold.(freqbands{fb}).bplim = bplim(fb,:);
  end
      mkdir(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm'))
      save(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
      clear unfold
end 
%%
% %2nd level analysis
 if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname',p.analysisname);
 end
    
freqbands       = {'alfa','beta'};
%freqbands       = {'alfa'};
for fb = 1:2%:length(freqbands)
    stimB       = [];
    bslcor      = [];
    stimBspl    = [];
    ufpredict   = [];
    
    for tk = p.subj
         cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/' model '/glm/' cfg_eeg.sujid '_' model],'unfold')
        unfold = unfold.(freqbands{fb});
         if  strcmp(model,'IM_STcsi_Sdxy') || strcmp(model,'IM_STcsi_SdxyCI') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
            pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
            pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         
            ufpredict       = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
         end
    
        if strcmp(model,'IM_STcsi_Sdxy') || strcmp(model,'IM_STcsi_SdxyCI') || strcmp(model,'S_orderPrecat_xdifspl_IM_STsc_eff')
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
        %  stimB       = cat(4,stimB,permute(unfold.(freqbands{fb}).beta./...
         %   repmat(mean(unfold.(freqbands{fb}).beta(:,find(unfold.(freqbands{fb}).times>bslcor(1) & ...
          %      unfold.(freqbands{fb}).times<bslcor(2)),:),2),1,size(unfold.(freqbands{fb}).beta,2),1),[1,3,2]));
        % FALTA aQUI      
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
    stattype = 'signpermT';
    mc       = 'cluster';
    E283_run2nd_save
end


%%
% resultALL = result;
for fb = 1:2
    load(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',['glmALL_' freqbands{fb} stattype '_' mc]),'result')
    if ~isempty(bslcor)
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],(freqbands{fb}));
    else
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],(freqbands{fb}));
    end

%     auxresult = resultALL.(freqbands{fb});
    if any(strfind(p.analysisname,'mirr'))
        plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.3  .7 .025],[],1);
    else
        plotBetasTopos(cfg_eeg,result,'mean',pathfig,[-.3  .7 .05],[],0);
    end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultSplinesALL = resultSplines;
freqbands       = {'alfa','beta'};
for fb = 1:2
    resultSplines = resultSplinesALL.(freqbands{fb});
    if ~isempty(bslcor)
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],['splines',freqbands{fb}]);
    else
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],['splines',freqbands{fb}]);
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
end
%%

%chnstoPlot      = {{'Cz','CPz'},{'Pz','POz'},{'O1','Oz','O2'},...
 %   {'P5','P7','PO7'},{'P6','P8','PO8'},{'C3','CP3'}};
chnstoPlot      = {{'FCz','Cz','CPz'},{'O1','Oz','O2'}};
Bnames        = {'image_Intercept','stim_2_Intercept','stim_side'};
axLim           = [-1 .9 -10 10];
for fb = 1:2
    result = resultALL.(freqbands{fb});
    plotBetasChannels(cfg_eeg,result,chnstoPlot,Bnames,axLim)
end