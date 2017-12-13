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
p.times_tflock              = [800 800];
p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
p.bsl                       = [-.8 -.5]; 
p.reref                     = 'yes';
p.keep                      = 'no';
p.collim                    = [0 2];
p.cfgTFR.channel            = 'all';	
p.cfgTFR.keeptrials         = 'no';	                
p.cfgTFR.method             = 'mtmconvol';
% p.cfgTFR.taper              = 'hanning'
% p.cfgTFR.width              = 5; 
p.cfgTFR.pad                = 3;
p.cfgTFR.output             = 'pow';	
p.cfgTFR.foi                = 4:1:40;	

p.cfgTFR.t_ftimwin          = 3./p.cfgTFR.foi;
p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
plottp(p.cfgTFR)

% p.cfgTFR.t_ftimwin          = .5*ones(1,length(p.cfgTFR.foi));
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
cfg_eeg.eegfolder       = ['/Volumes/nibaldo/trabajo/E283/data/' filename '/'];
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
   load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])
    load('/Users/jossando/trabajo/E283/analysis/eyedata/alleyedataFULL')                         
   eyedata.events.orderPreT = data.orderPreT(data.subject==tk);
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_fix'],'TFRav')
    % trial definitions
   [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'revisit','==0';'refix','==0';'onTarg','==0'},p.times_tflock,{+1,2,'origstart','>0'})
    
    at=1
%     [TFRav.fix] = getTFRsfromtrl({cfg_eeg},{trl},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%     
    [TFRav.fixL] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)<0),:)},...
        p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
    [TFRav.fixR] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)>0),:)},...
        p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
% 
%     [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
%             'revisit','==0';'refix','==0';'onTarg','==0'},p.times_tflock,{+2, 1, 'revisit', '==2';+2, 1, 'onTarg', '==0'}) ;
%    [TFRav.prerevisit] = getTFRsfromtrl({cfg_eeg},{trl},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%     
    for tfix = 0:3
        
         if tfix == 0, 
             oT =1;
          [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT',sprintf('==%d',tfix);'onTarg',sprintf('==%d',oT)},p.times_tflock);
         else
             oT=0;
             [trl,events]  = define_event(cfg_eeg,eyedata,1,{'&origstart',sprintf('>%d',0);...
            'orderPreT',sprintf('==%d',tfix);'onTarg',sprintf('==%d',oT)},p.times_tflock,{+1,2,'origstart','>0'});
         end
       
%       
%         [TFRav.(['fixn' num2str(tfix)])] = getTFRsfromtrl({cfg_eeg},{trl},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%       
        [TFRav.(['fixLn' num2str(tfix)])] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)<0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
       [TFRav.(['fixRn' num2str(tfix)])] = getTFRsfromtrl({cfg_eeg},{trl(find(events.posendx(2:2:end)-events.posinix(2:2:end)>0),:)},...
            p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR); 
    
    end
       save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_fix'],'TFRav','cfg_eeg','p')
end

%%
% grand averages
at                  = 1;
E283_params
p.bsl                       = [-.8 -.5]; 
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
freq1           = GA.prerevisit;
freq2          = GA.fix;
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';
difffreq        = ft_math(cfgs,freq1,freq2);


%%
% Fieldtrip fast plotting
 
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
% cfgp.baseline       = [-.75 .25] ;
%     cfgp.baselinetype   = 'db';
% cfgp.trials     = 51:70
%          cfgp.ylim           = [10 20];
          cfgp.xlim           = [-.4 .4];
%   cfgp.zlim           = [4 15]
%      cfgp.zlim           = [-1 1];
% %   data =TFRallt_RCsac.(p.analysis_type{at});
  data = difffreq;
% data.powspLFRallt_LU.ICAemUvsCa.powspctrm)4
       figure
        ft_multiplotTFR(cfgp,data)
%  cfgp.comment = 'yes'
%      figure,ft_topoplotTFR(cfgp,data)