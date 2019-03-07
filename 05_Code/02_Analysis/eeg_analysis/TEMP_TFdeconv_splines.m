clear
E283_params                                 % basic experimental parameters               % 
p.analysisname  = 'deconvTFCI';
%model           = 'IM_STcsi';
model           = 'IM_STcsi_SdxyCI';
%

if ismac    
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jossando/trabajo/E283/','analysisname', p.analysisname); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('expfolder','/Users/jpo/trabajo/E283/','analysisname',p.analysisname);
 end
    
freqbands       = {'alfa','beta'};
%freqbands       = {'alfa'};
xx = [];
yy = [];
for fb = 1%:2%:length(freqbands)
   stimB       = [];
    bslcor      = [-.4 -.1];
    stimBspl    = [];
    ufpredict   = [];
    for tk = p.subj
         cfg_eeg             = eeg_etParams_E283(cfg_eeg,'sujid',sprintf('s%02dvs',tk));
        load([cfg_eeg.eeganalysisfolder cfg_eeg.analysisname '/' model '/glm/' cfg_eeg.sujid '_' model],'unfold')
        unfold = unfold.(freqbands{fb});
        xx = [xx,unfold.unfold.splines{1}.paramValues(~isnan(unfold.unfold.splines{1}.paramValues))];
        yy = [yy,unfold.unfold.splines{2}.paramValues(~isnan(unfold.unfold.splines{2}.paramValues))];
         if  strcmp(model,'IM_STcsi_Sdxy') | strcmp(model,'IM_STcsi_SdxyCI')
            pxdM        = 10;   pxd_pred    = [-fliplr((2.^[-1:4]/2^4)*pxdM) 0 (2.^[-1:4]/2^4)*pxdM];
            pydM        = 10;   pyd_pred    = [-fliplr((2.^[-1:4]/2^4)*pydM) 0 (2.^[-1:4]/2^4)*pydM];
            pxd_pred    = [[-8.7:2:.7],0,[.7:2:8.7]];
            pyd_pred    = [[-4.38:1:.38],0,[.38:1:4.38]];
            ufpredict       = uf_predictContinuous(unfold,'predictAt',{{'pxdiff',pxd_pred},{'pydiff',pyd_pred}});
         end
%     
        if strcmp(model,'IM_STcsi_Sdxy')| strcmp(model,'IM_STcsi_SdxyCI')
            rmspl                   = strmatch('spline',{unfold.param.type});
            unfold.beta(:,:,rmspl)  = [];
            unfold.param(rmspl)     = [];

            rmnospl                   = setdiff(1:length(ufpredict.param),strmatch('spline',{ufpredict.param.type}));
            ufpredict.beta(:,:,rmnospl) = [];
            ufpredict.param(rmnospl)    = [];
        end
        if ~isempty(bslcor) 
            stimB       = cat(4,stimB,permute(unfold.beta-...
            repmat(mean(unfold.beta(:,find(unfold.times>bslcor(1) & unfold.times<bslcor(2)),:),2),1,size(unfold.beta,2),1),[1,3,2]));
        %  stimB       = cat(4,stimB,permute(unfold.(freqbands{fb}).beta./...
         %   repmat(mean(unfold.(freqbands{fb}).beta(:,find(unfold.(freqbands{fb}).times>bslcor(1) & ...
          %      unfold.(freqbands{fb}).times<bslcor(2)),:),2),1,size(unfold.(freqbands{fb}).beta,2),1),[1,3,2]));
        % FALTA aQUI      
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
%     end
    resultdiff        = regmodel2ndstat(stimB,unfold.times,elec,2000,stattype,mc);

    end
         
 stattype = 'signpermT';
    mc       = 'cluster';
    if ~isempty(stimBspl)
        coeffs_spl      = strcat([ufpredict.param.event]','_',strrep(strrep(strrep({ufpredict.param.name},':','XX'),'(',''),')','')','_',num2str([ufpredict.param.value]'));
        resultSplines.B       = stimBspl;
        resultSplines.coeffs  = coeffs_spl;
        resultSplines.splines = ufpredict.unfold.splines;
        resultSplines.times   = ufpredict.times;
    end
end