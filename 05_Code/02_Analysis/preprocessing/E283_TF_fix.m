% E275
% Basic TF analysis aligned to saccade
% - Simple TF charts
% - Selection fo peak frequencies at theta, alpha and beta bands
% - GlM analysis per frequency

%%
% TFR
% sge               = str2num(getenv('SGE_TASK_ID'));
clear
% Analysis parameters
E283_params
p.times_tflock              = [500 800];
p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
p.bsl                       = [-.375 -.125]; 
p.reref                     = 'yes';
p.keep                      = 'no';
p.collim                    = [0 2];
p.cfgTFR.channel            = 'all';	
p.cfgTFR.keeptrials         = 'no';	                
p.cfgTFR.foi                = 8:1:30;%2.^(3:.125:5);%(2.^([3:.25:5.25]));% %6:1:40	
p.cfgTFR.method             = 'mtmconvol';%'wavelet';%%
p.cfgTFR.taper              = 'hanning';%'dpss';
   
p.cfgTFR.pad                = 4;
p.cfgTFR.t_ftimwin          = .250*ones(1,length(p.cfgTFR.foi));      %    single taper

p.cfgTFR.output             = 'pow';	
p.cfgTFR.toi                = (-p.times_tflock(1):20:p.times_tflock(2))/1000;	
%%
for tk =    p.subj; % subject number

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
                                            'analysisname','fixTFR');       % single experiment/session parameters 

    
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])
    load('/Users/jossando/trabajo/E283/analysis/eyedata/alleyedataFULL')                         
   
   % load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_fix'],'TFRav')
    % trial definitions
%    [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
%             'revisit','==0';'refix','==0';'onTarg','==0'},p.times_tflock,...
%             {+1,2,'origstart','>0'})
    
    at=1
%     [TFRav.fixL] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)<0),:)},...
%         p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
%     [TFRav.fixR] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)>0),:)},...
%         p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
% 
%     [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
%             'revisit','==0';'refix','==0';'onTarg','==0'},p.times_tflock,{+2, 1, 'revisit', '==2';+2, 1, 'onTarg', '==0'}) ;
%    [TFRav.prerevisit] = getTFRsfromtrl({cfg_eeg},{trl},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%     
eyedata.events.nextToTarg = nansum(eyedata.events.nextToTarg)
    for tfix = 0:1
        
         if tfix == 0, 
             oT =1;
          [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT',sprintf('==%d',tfix);'onTarg',sprintf('==%d',oT)},p.times_tflock);
         else
             oT=0;
             [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT',sprintf('==%d',tfix);'onTarg','==0';'nextToTarg','==1'},p.times_tflock,{+1,2,'origstart','>0'});
         end
         
         [TFRav.(['fixLn' num2str(tfix)])] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)<0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
         [TFRav.(['fixRn' num2str(tfix)])] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)>0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
    
    end
     [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT','>1';'onTarg','==0';'nextToTarg','==1'},p.times_tflock,{+1,2,'origstart','>0'});
       [TFRav.nextLtarget] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)<0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
         [TFRav.nextRtarget] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)>0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);   
       save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_fix'],'TFRav','cfg_eeg','p')
end

%%
% grand averages
at                  = 1;
E283_params
p.bsl                =  [-.375 -.125]; 
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
cfgr                = [];
% p.bsl               = [-1 -.5]; 
cfgr.baseline       = p.bsl;
cfgr.baselinetype   = 'db';
ss=1;
for tk = p.subj; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','fixTFR','expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','fixTFR');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_fix'],'TFRav')
    fTFR    = fields(TFRav);
    for ff=1:length(fTFR)
        auxstr = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
         if isfield(auxstr,'elec')
             auxstr = rmfield(auxstr,'elec');
         end
        faux(ss,ff) = auxstr;
    end
    ss = ss +1;
end

for ff=1:length(fTFR)
    str_GA = 'GA.(fTFR{ff}) = ft_freqgrandaverage([]'
    for ss = 1:size(faux(:,ff),1)
        str_GA = [str_GA, ',faux(', num2str(ss), ',ff)'];
    end
    str_GA = [str_GA,');'];
    eval(str_GA)
end

%%
% cfgr                = [];
% cfgr.baseline       = p.bsl;
% cfgr.baselinetype   = 'db';
% [freq1]             = ft_freqbaseline(cfgr, TFRallt_Lsac.(p.analysis_type{at}));
% [freq2]             = ft_freqbaseline(cfgr, TFRallt_Rsac.(p.analysis_type{at}));
% 
% %%

cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';
GA.diff_fixRn1_nextLtarget        = ft_math(cfgs,GA.fixRn1,GA.nextLtarget);
GA.diff_fixRn1_nextRtarget        = ft_math(cfgs,GA.fixRn1,GA.nextRtarget);
GA.diff_fixLn1_nextLtarget        = ft_math(cfgs,GA.fixLn1,GA.nextLtarget);
GA.diff_fixLn1_nextRtarget        = ft_math(cfgs,GA.fixLn1,GA.nextRtarget);

%%
% Fieldtrip fast plotting
 
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
cfgp.baseline       = 'no'
%  cfgp.baseline       = [-.375 -.125]; ;
%      cfgp.baselinetype   = 'db';
% cfgp.trials     = 51:70
%          cfgp.ylim           = [10 20];
%           cfgp.xlim           = [-.4 .4];
%   cfgp.zlim           = [4 15]
%      cfgp.zlim           = [-1 1];
%   data = GA.diff_fixLn1_nextRtarget;
data = GA.fixRn0
  figure
        ft_multiplotTFR(cfgp,data)
%  cfgp.comment = 'yes'
%      figure,ft_topoplotTFR(cfgp,data)