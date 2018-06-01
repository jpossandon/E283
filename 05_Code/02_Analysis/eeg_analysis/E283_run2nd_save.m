if ~isempty(stimBspl)
    coeffs_spl      = strcat([ufpredict.param.event]','_',strrep(strrep(strrep({ufpredict.param.name},':','XX'),'(',''),')','')','_',num2str([ufpredict.param.value]'));
    resultSplines.beta    = stimBspl;
    resultSplines.coeffs  = coeffs_spl;
    resultSplines.splines = ufpredict.unfold.splines;
    resultSplines.times   = ufpredict.times;
end
load(cfg_eeg.chanfile)
result        = regmodel2ndstat(stimB,unfold.times,elec,1000,'signpermT','cluster');
coeffs_nospl  = strcat([unfold.param.event]','_',strrep(strrep(strrep({unfold.param.name},':','XX'),'(',''),')','')');
result.coeffs = coeffs_nospl;
mkdir(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm'))

if ~isempty(bslcor) 
    if ~isempty(stimBspl)
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALLbslcorr'),'result','resultSplines')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALLbslcorr'),'result')
    end
else
    if ~isempty(stimBspl)
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALL'),'result','resultSplines')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm','glmALL'),'result')
    end
end%