% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'phaseCI';
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
    [trl,events]           = define_event(cfg_eeg,eyedata,2,{'&origstart','>0';'&origstart','<1000'},...
                                [800 100],{-1,1,'origstart','>0'}); 
    epochevents             = [];
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.type        = cell(1,length(events.start));
    epochevents.type(events.type==1) = repmat({'fix'},1,sum(events.type==1));
    epochevents.type(events.type==2) = repmat({'sac'},1,sum(events.type==2));

    epochevents.pxini       = (events.posinix-960)/45;            
    epochevents.pyini       = (events.posiniy-540)/45;     
    epochevents.pxend       = (events.posendx-960)/45;            
    epochevents.pyend       = (events.posendy-540)/45;
    epochevents.pxdiff      = epochevents.pxend-epochevents.pxini;  
    epochevents.pydiff      = epochevents.pyend-epochevents.pyini; 
    epochevents.side        = nan(1,length(events.start));    
    epochevents.cross       = nan(1,length(events.start));    
    epochevents.inst        = nan(1,length(events.start)); 
    
    [trl,events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','>0'},[1500 1000]);
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
    fprintf('\nMinimum intertrial latency subject %d:%d',tk,min(diff(events.time(events.value==96))))

    epochevents.pxini       = [epochevents.pxini,nan(1,length(events.value))];    
    epochevents.pyini       = [epochevents.pyini,nan(1,length(events.value))];       
    epochevents.pxend       = [epochevents.pxend,nan(1,length(events.value))];          
    epochevents.pyend       = [epochevents.pyend,nan(1,length(events.value))];  
    epochevents.pxdiff      = [epochevents.pxdiff,nan(1,length(events.value))];  
    epochevents.pydiff      = [epochevents.pydiff,nan(1,length(events.value))];  
    % getting the data in EEGlab format
    [EEG,winrej] = getDataDeconv(cfg_eeg,epochevents,250,0);  
     mirindx         = mirrindex({EEG.chanlocs.labels},[cfg_eeg.analysisfolder '/01_Channels/mirror_chans']); 
       
    if any(strfind(p.analysisname,'CI'))
        LstimTimes = epochevents.latency(find(strcmp(epochevents.type,'stim') & epochevents.side == -1)-1);
        for tt  = 1:length(LstimTimes)
            EEGst = find(EEG.times>LstimTimes(tt),1);
            mirsamples = EEGst-floor(EEG.srate*.45):EEGst+ceil(EEG.srate*1.6);
            EEG.data(:,mirsamples) = EEG.data(mirindx,mirsamples);
            ETst  = find(epochevents.latency>LstimTimes(tt) & epochevents.latency<LstimTimes(tt)+1000 & ...
                (strcmp(epochevents.type,'sac') | strcmp(epochevents.type,'fix')));
            if ~isempty(ETst)
                [epochevents.pxini(ETst) epochevents.pxend(ETst) epochevents.pxdiff(ETst)] = deal(-epochevents.pxini(ETst),-epochevents.pxend(ETst),-epochevents.pxdiff(ETst));
            end
        end
    end
   
        
    freqbands       = {'alfa','beta'};
    bplim           = [9 15;16 25];
    filtPnts        = [368,254]; % check this
     LstimTimes = epochevents.latency(find(strcmp(epochevents.type,'stim')));
    side = epochevents.side(find(strcmp(epochevents.type,'stim')));
    cross = epochevents.cross(find(strcmp(epochevents.type,'stim')));
    info = epochevents.inst(find(strcmp(epochevents.type,'stim')));
    trSt            = strmatch('S 96',{EEG.event.type});
    evLat           = [EEG.event.latency];
    lat_trSt        = round(evLat(trSt));
    bsl             = [-.450 0];
    winsize         = .3/(1/EEG.srate);

    for fb = 1:length(freqbands)
      
      
        % this two are quite equivalent but holber path is faster
        EEGaux              = pop_eegfiltnew(EEG, bplim(fb,1), bplim(fb,2), [], 0, [], 0);
        hilb                = hilbert(EEGaux.data')';
        phasehilb           = angle(hilb);
        EEGaux.data         = abs(hilb).^2; %'power'
        for tt = 1:length(LstimTimes)
            ETst        = find(epochevents.latency>LstimTimes(tt) & epochevents.latency<LstimTimes(tt)+850 & ...
                (strcmp(epochevents.type,'sac')));
            stimSaux = find(EEG.times>=LstimTimes(tt),1,'First')
            phasedata.phase = [phasedata.phase,phasehilb(:,stimSaux)];
            phasedata.pow = [phasedata.phase,EEGaux.data(:,stimSaux)];
            sign(epochevents.pxend(ETst))
            epochevents.latency(ETst)-LstimTimes(tt)
        end
          %foi = bplim(fb,1):bplim(fb,2)
         % EEGaux= EEG;
%         for ch = 1:EEG.nbchan
%             [spectrum_mtmconvol,ntaper,foi,toi] = ft_specest_mtmconvol(EEG.data(ch,:),EEG.times/1000,...
%                 'timeoi',EEG.times/1000,'timwin',3./foi ,'taper','hanning',...
%                 'pad',[],'freqoi',foi,'polyorder',0,'dimord','chan_time_freqtap');
%             EEGaux.data(ch,:) =  sum(abs(spectrum_mtmconvol).^2,3);
%         end
        bslsample           = bsl.*EEGaux.srate;
       % sampltobsl          = [-.4 1.6].*EEGaux.srate;
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
freqbands       = {'alfa'};
for fb = 1:length(freqbands)
    stimB       = [];
    bslcor      = [];
    stimBspl    = [];
    ufpredict   = [];
    
    for tk = p.subj
         cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/' model '/glm/' cfg_eeg.sujid '_' model],'unfold')
        
         if  strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
            pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
            pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
         
            ufpredict.(freqbands{fb})       = uf_predictContinuous(unfold.(freqbands{fb}),'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
         end
    
        if strcmp(model,'IM_STcsi_Sdxy')| strcmp(model,'IM_STcsi_SdxyCI')
            rmspl                   = strmatch('spline',{unfold.(freqbands{fb}).param.type});
            unfold.(freqbands{fb}).beta(:,:,rmspl)  = [];
            unfold.(freqbands{fb}).param(rmspl)     = [];

            rmnospl                   = setdiff(1:length(ufpredict.(freqbands{fb}).param),strmatch('spline',{ufpredict.(freqbands{fb}).param.type}));
            ufpredict.(freqbands{fb}).beta(:,:,rmnospl) = [];
            ufpredict.(freqbands{fb}).param(rmnospl)    = [];
        end
        if ~isempty(bslcor) 
            stimB       = cat(4,stimB,permute(unfold.(freqbands{fb}).beta-...
            repmat(mean(unfold.(freqbands{fb}).beta(:,find(unfold.(freqbands{fb}).times>bslcor(1) & unfold.(freqbands{fb}).times<bslcor(2)),:),2),1,size(unfold.(freqbands{fb}).beta,2),1),[1,3,2]));
            if ~isempty(ufpredict)
                stimBspl    = cat(4,stimBspl,permute(ufpredict.(freqbands{fb}).beta-...
                    repmat(mean(ufpredict.(freqbands{fb}).beta(:,find(ufpredict.(freqbands{fb}).times>bslcor(1) & ufpredict.(freqbands{fb}).times<bslcor(2)),:),2),1,size(ufpredict.(freqbands{fb}).beta,2),1),[1,3,2]));
            end
        else
            stimB = cat(4,stimB,permute(unfold.(freqbands{fb}).beta(:,:,:),[1,3,2]));
            if ~isempty(ufpredict)
                stimBspl    = cat(4,stimBspl,permute(ufpredict.(freqbands{fb}).beta(:,:,:),[1,3,2]));
            end
        end
    end
    if strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
        coeffs_spl      = strcat([ufpredict.(freqbands{fb}).param.event]','_',strrep(strrep(strrep({ufpredict.(freqbands{fb}).param.name},':','XX'),'(',''),')','')','_',num2str([ufpredict.(freqbands{fb}).param.value]'));
        resultSplines.(freqbands{fb}).beta    = stimBspl;
        resultSplines.(freqbands{fb}).coeffs  = coeffs_spl;
        resultSplines.(freqbands{fb}).splines = ufpredict.(freqbands{fb}).unfold.splines;
        resultSplines.(freqbands{fb}).times   = ufpredict.(freqbands{fb}).times;
    end
    load(cfg_eeg.chanfile)
    result.(freqbands{fb}) = regmodel2ndstat(stimB,unfold.(freqbands{fb}).times,elec,1000,'signpermT','cluster');
    coeffs_nospl  = strcat([unfold.(freqbands{fb}).param.event]','_',strrep(strrep(strrep({unfold.(freqbands{fb}).param.name},':','XX'),'(',''),')','')');
    result.(freqbands{fb}).coeffs = coeffs_nospl;
end


mkdir(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm'))
if ~isempty(bslcor) 
    if strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALLbslcorr'),'result','resultSplines')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALLbslcorr'),'result')
    end
else
    if strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALL'),'result','resultSplines')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALL'),'result')
    end
end%

%%
resultALL = result;
for fb = 1%:2
    
    if ~isempty(bslcor)
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy') '_bslcorr'],(freqbands{fb}));
    else
        pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'figures',[datestr(now,'ddmmyy')],(freqbands{fb}));
    end

    auxresult = resultALL.(freqbands{fb});
    if any(strfind(p.analysisname,'mirr'))
        plotBetasTopos(cfg_eeg,auxresult,pathfig,[-.12  .6 .015],[],1);
    else
        plotBetasTopos(cfg_eeg,auxresult,pathfig,[-.3  .7 .05],[],0);
    end
end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultSplinesALL = resultSplines;
freqbands       = {'alfa','beta'};
for fb = 1:2
    resultSplines = resultSplinesALL.(freqbands{fb})
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