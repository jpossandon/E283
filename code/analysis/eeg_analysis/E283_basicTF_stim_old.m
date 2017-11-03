% E283
% Basic TF analysis aligned to stimulus
% - Simple TF charts
% - Selection fo peak frequencies at theta, alpha and beta bands
% - GlM analysis per frequency

%%
% TFR
% sge               = str2num(getenv('SGE_TASK_ID'));
clear
for tk = [2,4,5,6,7]; % subject number

    % Analysis parameters
    p.times_tflock              = [500 1500];
    p.analysis_type             = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
    p.bsl                       = [-.5 -.25]; 
    p.reref                     = 'yes';
    p.keep                      = 'yes';
    p.collim                    = [0 2];
    p.cfgTFR.channel            = 'all';	
    p.cfgTFR.keeptrials         = 'yes';	                
    p.cfgTFR.method             = 'mtmconvol';
    p.cfgTFR.taper              = 'hanning';
    % p.cfgTFR.width              = 5; 
    p.cfgTFR.output             = 'pow';	
    % p.cfgTFR.foi                = 4:2:35;	
    % p.cfgTFR.t_ftimwin          = .512*ones(1,length(p.cfgTFR.foi));
    p.cfgTFR.toi                = (-p.times_tflock(1):20:p.times_tflock(2))/1000;	

    p.cfgTFR.foi               = 4:1:40;	
    p.cfgTFR.pad                = 4;
    p.cfgTFR.t_ftimwin         = .250.*ones(1,length(p.cfgTFR.foi));
%     p.cfgTFR.t_ftimwin          = 3./p.cfgTFR.foi;
%     p.cfgTFR.tapsmofrq          = 0.5*p.cfgTFR.foi;
%     plottp(p.cfgTFR)


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
                                            'analysisname','stimlockTFR');       % single experiment/session parameters 

    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])                         
    
    fieldstoav      = {'LU_I','RU_I','LC_I','RC_I','LU_unI','RU_unI','LC_unI','RC_unI','IM','U_I','C_I','U_unI','C_unI'};
    trigs           = [1,2,5,6,9,10,13,14,96];
    for f = 1:9
        [trls.(fieldstoav{f}),events]  = define_event(cfg_eeg,eyedata,'ETtrigger',{'value',['==' num2str(trigs(f))]},p.times_tflock);            
    end
%     [trls.trlLU_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==1'},p.times_tflock);            
%     [trls.trlRU_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==2'},p.times_tflock);            
%     [trls.trlLC_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==5'},p.times_tflock);            
%     [trls.trlRC_I,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==6'},p.times_tflock);            
%     [trls.image,events]                = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==96'},p.times_tflock);  
    
    % Locked to stimulus
    at  = 1;
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])

    for f = 1:9
        eval(['[TFR_' fieldstoav{f} '] = getTFRsfromtrl({cfg_eeg},{trls.' fieldstoav{f} '},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);'])
     
    end
%     [TFR_LU_I] = getTFRsfromtrl({cfg_eeg},{trls.trlLU_I},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%     [TFR_RU_I] = getTFRsfromtrl({cfg_eeg},{trls.trlRU_I},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
% 
%     [TFR_LC_I] = getTFRsfromtrl({cfg_eeg},{trls.trlLC_I},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
%     [TFR_RC_I] = getTFRsfromtrl({cfg_eeg},{trls.trlRC_I},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);
% 
%     [TFR_IM] = getTFRsfromtrl({cfg_eeg},{trls.image},p.bsl,p.reref,p.analysis_type{at},p.keep,p.cfgTFR);

%      mirroring left stimulation data to make contra/ipsi plots
    mirindx         = mirrindex(TFR_LU_I.(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
    TFR_LU_I_mirr      = TFR_LU_I;
    TFR_LC_I_mirr      = TFR_LC_I;
    TFR_LU_I_mirr.(p.analysis_type{at}).powspctrm = TFR_LU_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR_LC_I_mirr.(p.analysis_type{at}).powspctrm = TFR_LC_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR_U_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RU_I.(p.analysis_type{1}),TFR_LU_I_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR_C_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RC_I.(p.analysis_type{1}),TFR_LC_I_mirr.(p.analysis_type{1})); %SAME

    TFR_LU_unI_mirr      = TFR_LU_unI;
    TFR_LC_unI_mirr      = TFR_LC_unI;
    TFR_LU_unI_mirr.(p.analysis_type{at}).powspctrm = TFR_LU_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR_LC_unI_mirr.(p.analysis_type{at}).powspctrm = TFR_LC_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR_U_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RU_unI.(p.analysis_type{1}),TFR_LU_unI_mirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR_C_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RC_unI.(p.analysis_type{1}),TFR_LC_unI_mirr.(p.analysis_type{1})); %SAME

    
    for f = 1:length(fieldstoav)
        TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], eval(['TFR_' fieldstoav{f} '.(p.analysis_type{1})']));
    end
    save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')

end

%%
% grand averages
at                  = 1;
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
cfgr                = [];
p.bsl               = [-.5 -.25]; 
cfgr.baseline       = p.bsl;
cfgr.baselinetype   = 'db';
s=1
for tk = [2,4,5,6,7]; % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockTFR','expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname','stimlockTFR');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_' p.analysis_type{at}],'TFRav','cfg_eeg','p')
    fTFR    = fields(TFRav);
    for ff=1:length(fTFR)
        faux(s,ff) = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
    end
    s=s+1;
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
[freq1]             = ft_freqbaseline(cfgr,GA.LU_I);
[freq2]             = ft_freqbaseline(cfgr,GA.LC_I);


% load(cfg_eeg.chanfile)
% statUC = freqpermBT(freq1,freq2,elec);

freq1av   = ft_freqdescriptives([], GA.U_unI);
freq2av   = ft_freqdescriptives([], GA.C_unI);
freq3av   = ft_freqdescriptives([], GA.U_I);
freq4av   = ft_freqdescriptives([], GA.C_I);


cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';
difffreq1       = ft_math(cfgs,freq1av,freq2av);
difffreq2       = ft_math(cfgs,freq3av,freq4av);
difffreq       = ft_math(cfgs,difffreq1,difffreq2);
difffreq       = ft_math(cfgs,GA.RU_unI,GA.RC_unI);

% difffreq.mask   = statUC.mask;
%%
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
% cfgp.trials     =4
%                  cfgp.baseline       = p.bsl;
%                   cfgp.baselinetype   = 'db';
% cfgp.ylim           = [0 40];
%  cfgp.xlim           = [-.75 1.25];
       cfgp.zlim           = [-2 2];
%     cfgp.maskparameter = 'mask';
%       cfgp.maskalpha = .3
%     data = TFRav.C_I.ICAem;
%       data = GA.LC_I;
        data =difffreq;

figure,ft_multiplotTFR(cfgp,data)

%%
% Comparisons stat
cfgr                = [];
cfgr.baseline       = p.bsl;
cfgr.baselinetype   = 'db';
[freq1]             = ft_freqbaseline(cfgr, TFR.RU.(p.analysis_type{at}));
[freq2]             = ft_freqbaseline(cfgr, TFR.RC.(p.analysis_type{at}));
    
% load(cfg_eeg.chanfile)
% statLR = freqpermBT(freq1,freq2,elec)
