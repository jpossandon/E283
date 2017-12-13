clear
E283_params
s=1;
%%
allFix = struct('fix0',struct('ERP',[],'durs',[]),'fix1',struct('ERP',[],'durs',[]),'fix2',struct('ERP',[],'durs',[]),'fix3',struct('ERP',[],'durs',[]),'fix4',struct('ERP',[],'durs',[]),'fix5',struct('ERP',[],'durs',[]));
for tk = p.subj
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end

    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname','fixLock');       % single experiment/session parameters 
    cfg_eeg.eegfolder       = ['/Volumes/nibaldo/trabajo/E283/data/' filename '/'];
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERPs/'])
   
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])
   
    % trial definitions
    p.bsl           = [-.25 0];
    p.rref          = 'yes';
    p.analysis_type = {'ICAem'};
    p.keep          = 'no';
    p.lpfreq        = 45;
    p.interval      = [-.2 .64 .02];
    p.colorlim      = [-6 6];
    p.ch_tplot      = [4,5,12,13,44,46];
    p.win           = [250 750];
   
    % revist and prerevisit lag 2
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==2';'refix','==0';'onTarg','==0'},p.win,{-2 1 'revisit' '==0'});
    [ERP.revisit2,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.revisit2.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' revisit lag2 / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_revisit_Lag2'],1)    
        
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==3';'refix','==0';'onTarg','==0'},p.win,{-2 1 'revisit' '==0'});
    [ERP.revisit3,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.revisit3.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' revisit lag3 / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_revisit_Lag3'],1)    
        
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','>4';'refix','==0';'onTarg','==0'},p.win,{-2 1 'revisit' '==0'});  
     [ERP.revisit4,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.revisit4.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' revisit lag4 / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_revisit_Lag4'],1)       
        
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==0';'refix','==0';'onTarg','==0'},p.win);
     [ERP.nore,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.nore.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' nore / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_nore'],1)       
        
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==0';'refix','==0';'onTarg','==0'},p.win,{+2 1 'revisit' '==2'}) ;
    [ERP.prerevisit2,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.prerevisit2.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' prerevisit2 / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_prerevisit2'],1)  
        
    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==0';'refix','==0';'onTarg','==0'},p.win,{+2 1 'refix' '==1'}) ;
     [ERP.prerefix,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
    fh = plot_topos(cfg_eeg,ERP.prerefix.(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' prerefix / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_prerefix'],1)     
    
    for ft = 0:5
        if ft == 0, oT =1;else,oT=0;end
        [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT',sprintf('==%d',ft);'onTarg',sprintf('==%d',oT)},p.win);
        p.keep          = 'yes';
        [ERPall,toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
        
        fixdurs             = events.dur;
        fixdurs(toelim{1})  = [];
        [Y,sortFix]         = sort(fixdurs);
        p.(['fix' num2str(ft)]) = fixdurs
         fh = newtplot(squeeze(mean(ERPall.(p.analysis_type{1}).trial(sortFix,p.ch_tplot,:),2)),...
                 ERPall.(p.analysis_type{1}).time,15,Y'/1000,p.colorlim,'ERPvsFix')
         doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_fixpreT_' num2str(ft) '_centralCH'],1)
        allFix.(['fix' num2str(ft)]).ERP        = cat(1,allFix.(['fix' num2str(ft)]).ERP,squeeze(mean(ERPall.(p.analysis_type{1}).trial(sortFix,p.ch_tplot,:),2)));
        allFix.(['fix' num2str(ft)]).durs        = [allFix.(['fix' num2str(ft)]).durs,Y];
        
         p.keep          = 'no';
        [ERP.(['fix' num2str(ft)]),toelim] = getERPsfromtrl({cfg_eeg},{trl},p.bsl,p.rref,p.analysis_type{1},p.lpfreq,p.keep);
         fh = plot_topos(cfg_eeg,ERP.(['fix' num2str(ft)]).(p.analysis_type{1}),p.interval,p.bsl,p.colorlim,[cfg_eeg.sujid ' fixpreT # ' num2str(ft) '/ bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
     
         doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],...
             'tiff',[ cfg_eeg.sujid '_fixpreT_' num2str(ft)],1)
    end
  save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERPs/' cfg_eeg.sujid '_fixTL'],'ERP','p')
end

%%
% allFix plots
% figure
% for ft = 0:5
%     [Y,sortFix]         = sort(allFix.(['fix',num2str(ft)]).durs);
%     fh = newtplot(allFix.(['fix',num2str(ft)]).ERP,...
%                  ERPall.(p.analysis_type{1}).time,50,Y'/1000,p.colorlim,'ERPvsFix')
% end

%grand average
clear
E283_params
at                  = 1;
p.analysis_type     = {'ICAem'}; 
cfgr                = [];
s=1
for tk = [6 7 11 12 15 19 22 20 21 23 26 27 28 29 31]; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','fixLock','expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','fixLock');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERPs/' cfg_eeg.sujid '_fixTL'],'ERP','p')
    % fix
    ERPall(s) = ERP;
    s=s+1;
end

fERP    = fields(ERPall);
for ff=1:length(fERP)
    str_GA = 'GA.(fERP{ff}) = ft_timelockgrandaverage([]';
    for ss = 1:length(ERPall)
        str_GA = [str_GA, ',ERPall(' num2str(ss) ').' fERP{ff} '.' p.analysis_type{1} ''];
    end
    str_GA = [str_GA,');'];
    eval(str_GA);
end
%%
% GA FIGURE
fERP    = fields(GA);
 p.interval      = [-.2 .64 .02];
p.colorlim  = [-6 6];
  p.bsl           = [-.25 -0];
for ff=1:length(fERP)

    fh            = plot_topos(cfg_eeg,GA.(fERP{ff}),p.interval,p.bsl,p.colorlim,['all ' fERP{ff}  ' / bsl: ' sprintf('%2.2f to %2.2f /',p.bsl(1),p.bsl(2))],1);
    doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/'],'tiff',['all_' fERP{ff}],1)
% saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/GA/all_imlock_' p.analysis_type{1} '_' fERP{ff}],'fig')
% close(fh)
end


%%
load(cfg_eeg.chanfile)
cfgp                = [];
cfgp.showlabels     = 'no'; 
cfgp.fontsize       = 12; 
cfgp.elec           = elec;
cfgp.interactive    = 'yes';
% cfgp.clim      = [-.6 .6];
figure,ft_multiplotER(cfgp,GA.nore,GA.fix0,GA.fix1,GA.fix2,GA.prerevisit2)