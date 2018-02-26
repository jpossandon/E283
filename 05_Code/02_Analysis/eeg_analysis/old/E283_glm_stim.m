% Time-series analysis locked to tactile stimulation
%%
% eeglab
clear
E283_params                                 % basic experimental parameters               % 
fmodel                      = 1;            % wich glm model
E283_models                                 % the models

%%
% subject configuration and data
stimB = []; stimB_mirr = [];
for tk = p.subj
%  for tk = p.subj;
% tk = str2num(getenv('SGE_TASK_ID'));
    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
                                            'filename',filename,...
                                            'EDFname',filename,...
                                            'event',[filename '.vmrk'],...
                                            'clean_name','final',...
                                            'analysisname','stimlock');    % single experiment/session parameters 
   
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata               

    mkdir([cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/'])
      E283_base_trl_event_def_stim                                        % trial configuration  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLM ANALYSIS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.npermute = 1;
p.analysis_typ = p.analysis_type{1}; 
% [modelos_stim]    = regmodelpermutef({cfg_eeg},p.analysis_type{1},{trls.left;trls.right},eval(p.model_cov),p.model_inter,p.bsl,p.rref,p.npermute,'effect',p.mirror);
[modelos_stim]    = regmodelpermutef({cfg_eeg},[],p);
stimB = cat(4,stimB,modelos_stim.B);

% save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/' cfg_eeg.sujid '/ERPglm_' cfg_eeg.sujid '_'  p.analysisname],'modelos_stim','p','trls','cfg_eeg');


end
%%
% % plotting betas each subject  
% % stimlock
% modelplot = modelos_stim;
% for b=1:size(modelplot.B,2)
%     betas.avg       = squeeze(modelplot.B(:,b,:));
%     collim          =[-6*std(betas.avg(:)) 6*std(betas.avg(:))];
%     betas.time      = modelplot.time;
%     betas.dof       = 1;
%     betas.n         = sum(modelplot.n);
%     
%     fh = plot_stat(cfg_eeg,modelplot.TCFEstat(b),betas,[],p.interval,collim,.05,sprintf('Beta:%s',p.coeff{b}),1);
%     saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/figures/' cfg_eeg.sujid '/' cfg_eeg.sujid '_glm_stimlock_' p.coeff{b} '_' p.analysis_type{1}],'fig')
%     close(fh)
% end

%%
%2nd level analysis

load(cfg_eeg.chanfile)
result = regmodel2ndstat(stimB,modelos_stim.time,elec,1000,'signpermT','cluster');
save([cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/' datestr(now,'ddmmyy') '_ERPglm_' p.analysisname],'result','stimB','p','cfg_eeg');

%%
Bs      = stimB;
res     = result;
for b=1:size(res.B,2)
    betas.dof   = 1;
    betas.n     = size(Bs,4);
    betas.avg   = squeeze(mean(Bs(:,b,:,:),4));
    betas.time  = res.clusters(1).time;
    collim      =[-6*std(betas.avg(:)) 6*std(betas.avg(:))]; 
    fh = plot_stat(cfg_eeg,result.clusters(b),betas,[],p.interval,collim,.05,sprintf('Beta:%s',p.coeff{b}),1);
 doimage(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/figures/'],'png',[ datestr(now,'ddmmyy') 'glm_' p.coeff{b} '_' p.analysisname],1)

    %     saveas(fh,[cfg_eeg.analysisfolder cfg_eeg.analysisname '/ERP/figures/' datestr(now,'ddmmyy') 'glm_' p.coeff{b} '_' p.analysisname],'fig')
%     close(fh)
end