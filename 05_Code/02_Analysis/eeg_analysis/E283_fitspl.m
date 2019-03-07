clear
E283_params                                  % basic experimental parameters               %
p.analysisname  = 'TL_dc_preT';
if ismac
    run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')
else
    run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')
end

if ismac
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
else
    cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname', p.analysisname);
end

stimB = [];
model = 'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff';
%  model               = 'S_preT_onT_T0_T1_T2_T3_pre_xdifspl_IM_STsc';
%bslcor      = [-.8 -.4];
  bslcor      = [];
stimSpl     = [];
stimBspl    = [];
ufpredict   = [];
binsize     = 1.5;
 bind        = [0:binsize:15]+binsize/2;
for tk = p.subj
    cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
    load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
   
    
    if  strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
            pDToTatg  = bind ;
       if strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{'3_DToTarg',pDToTatg}});
        
       end
      
    end
    
  
    % remove splines from coefficients estimates
    if any(strfind(model,'spl'))
         rmnospl                   = setdiff(1:length(ufpredict.param),strmatch('spline',{ufpredict.param.type}));
       
          ufpredict.beta(:,:,rmnospl) = [];
        ufpredict.param(rmnospl)    = [];
    end
    if ~isempty(bslcor)
         if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta-...
                repmat(mean(ufpredict.beta(:,find(ufpredict.times>bslcor(1) & ufpredict.times<bslcor(2)),:),2),1,size(ufpredict.beta,2),1),[1,3,2]));
        end
    else
        stimB = cat(4,stimB,permute(unfold.beta(:,:,:),[1,3,2]));
        if ~isempty(ufpredict)
            stimBspl    = cat(4,stimBspl,permute(ufpredict.beta(:,:,:),[1,3,2]));
        end
    end
    
end
coeffs_spl      = strcat([ufpredict.param.event]','_',strrep(strrep(strrep({ufpredict.param.name},':','XX'),'(',''),')','')','_',num2str([ufpredict.param.value]'));
resultSplines.B       = stimBspl;
resultSplines.coeffs  = coeffs_spl;
resultSplines.splines = ufpredict.unfold.splines;
resultSplines.times   = ufpredict.times;

%%
spltofit  = 'fixpre_3_DToTarg_' 
chtofit = {{'P1','Pz','P2','PO3','POz','PO4'}}
timetofit = [.4 .6];
load(cfg_eeg.chanlocs)
for ch = 1:length(chtofit)
    for ss = 1:size(resultSplines.B,4)

        auxChns         = find(ismember({chanlocs.labels},chtofit{ch}));
        spls        = strmatch(spltofit,resultSplines.coeffs);
        times       = find(resultSplines.times>timetofit(1) & resultSplines.times<timetofit(2));
        auxdata     = nansum(mean(resultSplines.B(auxChns,spls,times,ss)),3);
        fitdata     = [bind' [(auxdata-min(auxdata))/range(auxdata)*100]' 100*ones(length(spls),1)];
        options.sigmoidName = 'norm';   % choose a cumulative Gauss as the sigmoid  
        options.expType     = 'YesNo';
        resultfitsplsuj(ss) = psignifit(fitdata,options);
    end
end
    load('pbyorderfit.mat')
allfit = [resultfitsuj.Fit]
allfitspl = [resultfitsplsuj.Fit];
figure,plot(allfitspl(1,:),allfit(1,:),'.')