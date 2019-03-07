clear
E283_params                                 % basic experimental parameters               %
p.analysisname  = 'TL_dc_preT';
model           = 'S_orderPrecat_xdifspl_IM_STsc_eff';
cfg_eeg         = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop

glmpath = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname,model,'glm','glmALL_boottrimet_cluster')
cmap3 = cbrewer('qual','Set1',6);
for tk = p.subj(9:end)
    cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
        'filename',filename,...
        'EDFname',filename,...
        'event',[filename '.vmrk'],...
        'clean_name','final',...
        'analysisname',p.analysisname);    % single experiment/session parameters
    
    % get relevant epochevents
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata
    nmovs = 7;
    eyedata.events.nextToTarg = any([eyedata.events.nextToTargH;eyedata.events.nextToTargV;eyedata.events.nextToTargD]);
      
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT','==0'},...
        [800 100],{-1,2,'origstart','>0';-2,1,'onTarg','==0'});
    
    events.prefixdur        = nan(1,length(events.start));
    events.prefixdur(2:3:end)  = events.dur(1:3:end);
    events = struct_elim(events,[1:3:length(events.type)],2);
    events = struct_elim(events,[2:2:length(events.type)],2);
    prefixdur = events.prefixdur;
    [ERPall,toelim]     = getERPsfromtrl({cfg_eeg},{trl},[-.8 -.6],'no','ICAem',42,'yes');
    prefixdur(toelim{1}) = [];
    load(glmpath,'result')
    indxB = strmatch('sacpre_orderPreT_0',result.coeffs);
    kernel = squeeze(mean(mean(result.B(:,indxB,find(result.clusters(1).time>-.1 & result.clusters(1).time<0),:),3),4))';
    for tt=1:size(ERPall.ICAem.trial,1)
        preTcomp(tt,:) = kernel*squeeze(ERPall.ICAem.trial(tt,:,:));
    end
    figure
    sujquant = [0 prctile(prefixdur,[20 40 60 80]) max(prefixdur)];
    plot(ERPall.ICAem.time,mean(preTcomp),'k','LineWidth',2), hold on
    for sp=1:5
        auxindx = find(prefixdur>sujquant(sp)&prefixdur<sujquant(sp+1));
    h(sp) = plot(ERPall.ICAem.time,mean(preTcomp(auxindx,:)),'Color',cmap3(sp,:),'LineWidth',2);
    M{sp}= num2str(round(mean(prefixdur(auxindx))));
    end
    legend(h,M,'Location','SouthWest')
    hline(0)
    vline(0)
end
  