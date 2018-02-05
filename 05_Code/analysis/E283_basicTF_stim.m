% E283
% Basic TF analysis aligned to stimulus
% - Simple TF charts
% - Selection fo peak frequencies at theta, alpha and beta bands
% - GlM analysis per frequency
clear
E283_params                                 % basic experimental parameters               % 
fmodel                      = 11;            % wich glm model
E283_models 
%%
% TFR
% sge               = str2num(getenv('SGE_TASK_ID'));
stimB_alpha = []; stimB_alpha_mirr = [];
stimB_beta = [];  stimB_beta_mirr = [];
for tk = p.subj; % subject number

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

    E283_base_trl_event_def_stim                                        % trial configuration  

    
    % Locked to stimulus
    at  = 1;
    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/TFR_' p.analysis_type{at} '/'])

    for e = 1:8
        [aux,toelim]= getTFRsfromtrl({cfg_eeg},{trls.(p.trls_stim{e})},p.bsl,p.rref,p.analysis_type{at},p.keep,p.cfgTFR);
         aux.ICAem.toelim = toelim;
         eval(['TFR_' p.trls_stim{e} '=aux;'])
    end
    %      mirroring left stimulation data to make contra/ipsi plots
    mirindx         = mirrindex(TFR_LU_I.(p.analysis_type{1}).label,[cfg_eeg.expfolder '/channels/mirror_chans']); 
    TFR_LU_Imirr      = TFR_LU_I;
    TFR_LC_Imirr      = TFR_LC_I;
    TFR_LU_Imirr.(p.analysis_type{at}).powspctrm = TFR_LU_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR_LC_Imirr.(p.analysis_type{at}).powspctrm = TFR_LC_I.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR_U_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RU_I.(p.analysis_type{1}),TFR_LU_Imirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR_C_I.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RC_I.(p.analysis_type{1}),TFR_LC_Imirr.(p.analysis_type{1})); %SAME

    TFR_LU_unImirr      = TFR_LU_unI;
    TFR_LC_unImirr      = TFR_LC_unI;
    TFR_LU_unImirr.(p.analysis_type{at}).powspctrm = TFR_LU_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);
    TFR_LC_unImirr.(p.analysis_type{at}).powspctrm = TFR_LC_unI.(p.analysis_type{at}).powspctrm(:,mirindx,:,:);

    cfgs            = [];
    cfgs.parameter  = 'powspctrm';
    TFR_U_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RU_unI.(p.analysis_type{1}),TFR_LU_unImirr.(p.analysis_type{1})); % ERASEME: there was an error here apeenof TFR.RU with TFR.LU instead of TFR.LUmirr
    TFR_C_unI.(p.analysis_type{1})           = ft_appendfreq(cfgs, TFR_RC_unI.(p.analysis_type{1}),TFR_LC_unImirr.(p.analysis_type{1})); %SAME

    fieldstoav      = {'LU_I','RU_I','LC_I','RC_I','LU_unI','RU_unI','LC_unI','RC_unI','LU_Imirr','LC_Imirr','U_I','C_I',...
        'LU_unImirr','LC_unImirr','U_unI','C_unI'};
    for f = 1:length(fieldstoav)
        TFRav.(fieldstoav{f}).(p.analysis_type{1}) = ft_freqdescriptives([], eval(['TFR_' fieldstoav{f} '.(p.analysis_type{1})']));
    end
    save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid p.analysisname],'TFRav','cfg_eeg','p')            

    % betas for glmb by freq bands
    freqsbands = [8 13;18 25];
    bandslabel = {'alpha','beta'};
    
    cfgr.baseline     = p.bsl;
    cfgr.baselinetype = 'relative';
    
    for fb = 1:size(freqsbands,1)
        for f = 1:8
            aux = eval(['TFR_' fieldstoav{f}]);
            indxfreqs                           = find(aux.(p.analysis_type{1}).freq>freqsbands(fb,1) & aux.(p.analysis_type{1}).freq<freqsbands(fb,2));
            aux.(p.analysis_type{1}).trial      = squeeze(nanmean(aux.(p.analysis_type{1}).powspctrm(:,:,indxfreqs,:),3));
            aux.(p.analysis_type{1})            = rmfield(aux.(p.analysis_type{1}),'powspctrm');
            aux.(p.analysis_type{1}).freq       = mean(freqsbands(fb,:));
            aux.(p.analysis_type{1}).dimord     = 'rpt_chan_time';
            meanVals                            = repmat(nanmean(aux.(p.analysis_type{1}).trial(:,:,[aux.(p.analysis_type{1}).time>p.bsl(1) & aux.(p.analysis_type{1}).time<p.bsl(2)],:), 3), [1 1 size(aux.(p.analysis_type{1}).trial,3) 1]);
            aux.(p.analysis_type{1}).trial      = aux.(p.analysis_type{1}).trial-meanVals;
            TFRaux.(bandslabel{fb}).(fieldstoav{f}) = aux;
            aux.(p.analysis_type{1}).avg        = squeeze(mean(aux.(p.analysis_type{1}).trial));
            fh = plot_topos(cfg_eeg, aux.(p.analysis_type{1}),p.interval,[],[-4 4],[filename ' ' bandslabel{fb} ' ' fieldstoav{f}]);
            doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'],'png',[bandslabel{fb} '_' fieldstoav{f}],1)
%             saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' filename '_' bandslabel{fb} '_' fieldstoav{f}],'fig')
%             close(fh)
        end
        
        TFR.(bandslabel{fb})(1) = aux;
        TFR.(bandslabel{fb})(1).(p.analysis_type{1}).trial = cat(1,TFRaux.(bandslabel{fb}).LU_I.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).LC_I.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).LU_unI.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).LC_unI.(p.analysis_type{1}).trial);
        TFR.(bandslabel{fb})(1).(p.analysis_type{1}).toelim = {[TFRaux.(bandslabel{fb}).LU_I.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.LU_I,1)+TFRaux.(bandslabel{fb}).LC_I.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.LU_I,1)+size(trls.LC_I,1)+TFRaux.(bandslabel{fb}).LU_unI.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.LU_I,1)+size(trls.LC_I,1)+size(trls.LU_unI,1)+TFRaux.(bandslabel{fb}).LC_unI.(p.analysis_type{1}).toelim{:}]};
        
        TFR.(bandslabel{fb})(2) = aux;
        TFR.(bandslabel{fb})(2).(p.analysis_type{1}).trial = cat(1,TFRaux.(bandslabel{fb}).RU_I.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).RC_I.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).RU_unI.(p.analysis_type{1}).trial,...
                                                                    TFRaux.(bandslabel{fb}).RC_unI.(p.analysis_type{1}).trial);
        TFR.(bandslabel{fb})(2).(p.analysis_type{1}).toelim = {[TFRaux.(bandslabel{fb}).RU_I.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.RU_I,1)+TFRaux.(bandslabel{fb}).RC_I.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.RU_I,1)+size(trls.RC_I,1)+TFRaux.(bandslabel{fb}).RU_unI.(p.analysis_type{1}).toelim{:};...
                                                                size(trls.RU_I,1)+size(trls.RC_I,1)+size(trls.RU_unI,1)+TFRaux.(bandslabel{fb}).RC_unI.(p.analysis_type{1}).toelim{:}]};
    end
    p.mirror         = [];
    [modelos_stim]   = regmodelpermutef({cfg_eeg},TFR.alpha,p);
    stimB_alpha      = cat(4,stimB_alpha,modelos_stim.B);
    [modelos_stim]   = regmodelpermutef({cfg_eeg},TFR.beta,p);
    stimB_beta       = cat(4,stimB_beta,modelos_stim.B);
    tiempo = modelos_stim.time';
  
   p.mirror         = [1 0];
   [modelos_stim_mirr]   = regmodelpermutef({cfg_eeg},TFR.alpha,p);
   stimB_alpha_mirr      = cat(4,stimB_alpha_mirr,modelos_stim_mirr.B);
   [modelos_stim_mirr]   = regmodelpermutef({cfg_eeg},TFR.beta,p);
   stimB_beta_mirr       = cat(4,stimB_beta_mirr,modelos_stim_mirr.B);
   save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' datestr(now,'ddmmyy') 'betas_' p.analysisname],'stimB_alpha','stimB_beta','stimB_alpha_mirr','stimB_beta_mirr','cfg_eeg','tiempo','p')            

end
 
%%
%2nd level analysis glm


% meanVals = repmat(nanmean(stimB_alpha(:,:,[modelos_stim.time>p.bsl(1) & modelos_stim.time<p.bsl(2)],:), 3), [1 1 size(stimB_alpha,3) 1]);
% stimB_alpha = stimB_alpha-meanVals;
% meanVals = repmat(nanmean(stimB_beta(:,:,[modelos_stim.time>p.bsl(1) & modelos_stim.time<p.bsl(2)],:), 3), [1 1 size(stimB_beta,3) 1]);
% stimB_beta = stimB_beta-meanVals;

load(cfg_eeg.chanfile)
result_alpha = regmodel2ndstat(stimB_alpha,tiempo,elec,1000,'signpermT','cluster');
result_beta = regmodel2ndstat(stimB_beta,tiempo,elec,1000,'signpermT','cluster');
result_alpha_mirr = regmodel2ndstat(stimB_alpha_mirr,modelos_stim_mirr.time,elec,1000,'signpermT','cluster');
result_beta_mirr = regmodel2ndstat(stimB_beta_mirr,modelos_stim_mirr.time,elec,1000,'signpermT','cluster');
save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/' datestr(now,'ddmmyy') '_TFRglm_' p.analysisname],'result_alpha','stimB_alpha','result_beta','stimB_beta','result_alpha_mirr','stimB_alpha_mirr','result_beta_mirr','stimB_beta_mirr','p','cfg_eeg','tiempo');

%%
result = result_alpha;
stimB  = stimB_alpha;
for b=1:size(result.B,2)
    betas.dof   = 1;
    betas.n     = size(stimB,4);
    betas.avg   = squeeze(median(stimB(:,b,:,:),4));
    betas.time  = tiempo;
    collim      =[-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
    fh = plot_stat(cfg_eeg,result.clusters(b),betas,[],p.interval,collim,.05,sprintf('Beta:%s',p.coeff{b}),1);
% doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/'],'png',[datestr(now,'ddmmyy') 'glm' p.coeff{b} '_' p.analysisname],1)

    %     saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/figures/' datestr(now,'ddmmyy') 'glm_' p.coeff{b} '_' p.analysisname],'fig')
%     close(fh)
end

%%
% plot means of the different 8 condition
% LR cross inf LRxCross LRxInf CrossxInf LRxCrossxInf
% -1  -1   -1     1       1        1        -1
% -1  -1    1     1      -1       -1         1
% -1   1   -1    -1       1       -1         1
% -1   1    1    -1      -1        1        -1
%  1  -1   -1    -1      -1        1         1
%  1  -1    1    -1       1       -1        -1
%  1   1   -1     1      -1       -1        -1
%  1   1    1     1       1        1         1
condmx = [1 -1 -1 -1  1  1  1 -1;...
          1 -1 -1  1  1 -1 -1  1;...
          1 -1  1 -1 -1  1 -1  1;...
          1 -1  1  1 -1 -1  1 -1;...
          1  1 -1 -1 -1 -1  1  1;...
          1  1 -1  1 -1  1 -1 -1;...
          1  1  1 -1  1 -1 -1 -1;...
          1  1  1  1  1  1  1  1];
      
[ch,~,tims,ss] = size(stimB);
for c = 1:size(condmx,1)
    data.dof   = 1;
    data.n     = size(stimB,4);
    data.time  = tiempo;
    data.avg = mean(sum(stimB.*repmat(condmx(c,:),[ch,1,tims,ss]),2),4);
     plot_stat(cfg_eeg,result.clusters(end),data,[],p.interval,[-4 4],.05,sprintf('LR(%d) Cross(%d) Inf(%d)',condmx(c,2),condmx(c,3),condmx(c,4)),1);
end          

%%
% grand averages
at                  = 1;
p.analysis_type     = {'ICAem'}; %'plain' / 'ICAe' / 'ICAm' / 'ICAem' 
cfgr                = [];
p.bsl               = [-.75 -.25]; 
cfgr.baseline       = p.bsl;
cfgr.baselinetype   = 'db';
subj=  [2 4 5 6 7 11 12 15 19 20 21 22 23 24]
for tk = 1:length(subj); % subject number
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',subj(tk)),'analysisname','stimlockTFR','expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',subj(tk)),'analysisname','stimlockTFR');
    end
    load([cfg_eeg.analysisfolder cfg_eeg.analysisname '/tfr/'  cfg_eeg.sujid p.analysisname],'TFRav','cfg_eeg','p')
    fTFR    = fields(TFRav);
    for ff=1:length(fTFR)
        faux(tk,ff) = ft_freqbaseline(cfgr,TFRav.(fTFR{ff}).(p.analysis_type{1}));
    end
end

for ff=1:length(fTFR)
    str_GA = 'GA.(fTFR{ff}) = ft_freqgrandaverage([]';
    for ss = 1:size(faux(:,ff),1)
        str_GA = [str_GA, ',faux(', num2str(ss), ',ff)'];
    end
    str_GA = [str_GA,');'];
    eval(str_GA)
end
%%
[freq1]             = ft_freqbaseline(cfgr,TFR_U.ICAem);
[freq2]             = ft_freqbaseline(cfgr, TFR_C.ICAem);

% %%
% load(cfg_eeg.chanfile)
% statUC = freqpermBT(freq1,freq2,elec);

freq1av   = ft_freqdescriptives([], freq1);
freq2av   = ft_freqdescriptives([], freq2);

%%
cfgs            = [];
cfgs.parameter  = 'powspctrm';
cfgs.operation  = 'subtract';
difffreq       = ft_math(cfgs,GA.RU_I,GA.RC_I);
% difffreq.mask   = statUC.mask;
%%
load(cfg_eeg.chanfile)
cfgp            = [];
cfgp.showlabels = 'no'; 
cfgp.fontsize   = 12; 
cfgp.elec       = elec;
cfgp.interactive    = 'yes';
% % cfgp.trials     =4
%            cfgp.baseline       = p.bsl;
%            cfgp.baselinetype   = 'db';
% cfgp.ylim           = [0 40];
 cfgp.xlim           = [-.75 1.25];
       cfgp.zlim           = [-4 4];
%     cfgp.maskparameter = 'mask';
%       cfgp.maskalpha = .3
    data =GA.LC_I;
% 
 figure,ft_multiplotTFR(cfgp,data)
% 
% %%
% % Comparisons stat
% cfgr                = [];
% cfgr.baseline       = p.bsl;
% cfgr.baselinetype   = 'db';
% [freq1]             = ft_freqbaseline(cfgr, TFR.RU.(p.analysis_type{at}));
% [freq2]             = ft_freqbaseline(cfgr, TFR.RC.(p.analysis_type{at}));
%     
% % load(cfg_eeg.chanfile)
% % statLR = freqpermBT(freq1,freq2,elec)
