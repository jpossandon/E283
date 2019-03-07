clear
E283_params
p.analysisname = 'TF_Stim_hann';
bands               = {'alpha','beta'};
MACpath = '/Users/jossando/trabajo/E283/';
at=1
for bb =1%:length(bands)
    stimB = [];
    for tk = p.subj; % subject number
        if ismac    
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname,'expfolder',MACpath); % this is just to being able to do analysis at work and with my laptop
        else
            cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'analysisname',p.analysisname);
        end
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/tfr/' cfg_eeg.sujid '_tfr_stim_ICAem'],'sujreg')
        stimB = cat(4,stimB,sujreg.(bands{bb}).B);
    end
    load(cfg_eeg.chanfile)
    result.(bands{bb})        = regmodel2ndstat(stimB,sujreg.times,elec,2000,'signpermT','cluster');
    result.(bands{bb}).coeffs = {'Intercept','HandCross','Info','Inter'};
model = 'STci_glm'
pathfig = fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,bands{bb},'figures',[datestr(now,'ddmmyy')]);
plotBetasTopos(cfg_eeg,result.(bands{bb}),pathfig,[-.12  .6 .03],[],0)
end

%%
