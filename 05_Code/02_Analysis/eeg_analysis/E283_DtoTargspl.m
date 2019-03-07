clear
chnstoPlot      = {{'P1','Pz','P2','PO3','POz','PO4'}}
E283_params                                  % basic experimental parameters               %
p.analysisname  = 'TL_dc_preT';
if ismac
    run('/Users/jossando/trabajo/matlab/unfold/init_unfold.m')
else
    run('/Users/jpo/trabajo/matlab/unfold/init_unfold.m')
end

model = 'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff';

if ismac
    load('/Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/alleyedataFULLevidence')
else
    load('/Users/jpo/trabajo/E283/07_Analysis/03_Eye/eyedata/alleyedataFULLevidence')
end
% figure, hold on
ss = 1;
alldata = [];
tint0   = nan(76,800);
   
bslcor      = [-.8 -.4];
% bslcor = [];
for tk = p.subj
    tk
     if ismac
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),...
            'expfolder','/Users/jpo/trabajo/E283/');
    end
    
    filename                = sprintf('s%02dvs',tk);
    cfg_eeg                 = eeg_etParams_E283(cfg_eeg,...
        'filename',filename,...
        'EDFname',filename,...
        'event',[filename '.vmrk'],...
        'clean_name','final',...
        'analysisname',p.analysisname);    % single experiment/session parameters
    load(cfg_eeg.chanlocs)
    auxChns         = find(ismember({chanlocs.labels},chnstoPlot{1}));
    % get relevant epochevents
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])            % eyedata
    % check eyedad is the same as full data
       nmovs = 4;
    eyedata.events.nextToTarg = any([eyedata.events.nextToTargH;eyedata.events.nextToTargV;eyedata.events.nextToTargD]);
    [trl,events]           = define_event(cfg_eeg,eyedata,1,{'&origstart','>0';'orderPreT',['<' num2str(nmovs)]},...
        [800 100],{-1,2,'origstart','>0';-2,1,'origstart','>0'});
    events.DToTarg = events.DToTarg./45; 
    events.onTarg(2:3:end)  = events.onTarg(1:3:end);
    events.DToTarg(2:3:end) = events.DToTarg(1:3:end);  
    events = struct_elim(events,[1:3:length(events.type)],2);
    epochevents.latency     = events.start;                       % fixation start, here the important thing is the ini pos
    epochevents.dur         = events.dur;
    epochevents.type        = cell(1,length(events.start));
    epochevents.trial      = events.trial;
    epochevents.orderPreT   =  events.orderPreT;
    epochevents.DToTarg     = events.DToTarg;
   
   
   epochevents = struct_elim(epochevents,[1:2:length(epochevents.type)],2);
   
   load(fullfile(cfg_eeg.eeganalysisfolder,cfg_eeg.analysisname,model,'glm',[cfg_eeg.sujid,'_',model]),'unfold')
   

   if strcmp(model,'F_T0_T1_Tx_DToTargspl_spldiffxy_IM_STsc_eff')
        pDToTatg = [.1:.1:30];
        splnames = {'DToTarg','3_DToTarg'};
        for sn = 1:length(splnames)
            ufpredict               = uf_predictContinuous(unfold,'predictAt',...
                {{splnames{sn},pDToTatg}});
            if ~isempty(bslcor)
            ufpredict.beta = ufpredict.beta-...
                repmat(mean(ufpredict.beta(:,find(ufpredict.times>bslcor(1) & ufpredict.times<bslcor(2)),:),2),1,size(ufpredict.beta,2),1);
            end
            splines(sn).B = ufpredict.beta(:,:,strmatch(splnames{sn},{ufpredict.param.name}));
            splines(sn).DToTatg = pDToTatg;
            splines(sn).name = splnames{sn};
        end
        intnames = {'(Intercept)','2_(Intercept)','3_(Intercept)'};
        for sn = 1:length(intnames)
            intercept(sn).B = unfold.beta(:,:,strmatch(intnames{sn},{unfold.param.name}));
            intercept(sn).name = intnames{sn};
        end
   end
%    checksplines(ss,:) = mean(splines(2).B(auxChns,:,5));
   ntrials = unique(epochevents.trial);
   tdata   = nan(76,800,3,length(ntrials));
   tint    = nan(76,800,3,length(ntrials));
   tint0(:,601-200:601+199,ss)   = intercept(1).B;
        
    for tt = 1:length(ntrials)
        indxtr  = find(epochevents.trial==ntrials(tt));
        indxadd = fliplr(round(cumsum(diff(epochevents.latency(indxtr)),'reverse')/4)); % TODO check so it does not need rounding
        for nfix = 1:length(indxadd)
            lastsaclag = indxadd(1)-round(epochevents.dur(indxtr(end-1))/4);
        
            if indxadd(nfix)<400
            if nfix ==1
                thispl = find(splines(1).DToTatg>epochevents.DToTarg(indxtr(end-nfix)),1,'first');
                  tdata(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt)   = nansum(cat(3,tdata(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt),splines(1).B(:,:,thispl)),3);
                 tint(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt)   = nansum(cat(3,tint(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt),intercept(2).B),3);

                  %                  tdata(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt)   = nansum(cat(3,tdata(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt),splines(1).B(:,:,thispl)),3);
%                  tint(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt)   = nansum(cat(3,tint(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt),intercept(2).B),3);
%                 plot(601-indxadd(nfix)-200:601-indxadd(nfix)+199,splines(1).B(4,:,thispl))
            else
                thispl = find(splines(2).DToTatg>epochevents.DToTarg(indxtr(end-nfix)),1,'first');
                  tdata(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt)   = nansum(cat(3,tdata(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt),splines(2).B(:,:,thispl)),3);
                    tint(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt)   = nansum(cat(3,tint(:,601-indxadd(nfix)-200:601-indxadd(nfix)+199,nfix,tt),intercept(3).B),3);
           
                  %                  tdata(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt)   = nansum(cat(3,tdata(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt),splines(2).B(:,:,thispl)),3);
                tint(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt)   = nansum(cat(3,tint(:,601-indxadd(nfix)+lastsaclag-200:601-indxadd(nfix)+lastsaclag+199,nfix,tt),intercept(3).B),3);
            end
            else
                indxadd
            end
        end
    end
     alldata(:,ss) = nanmean(nansum(nanmean(tdata(auxChns,:,:,:)),3),4);
     alldata2(:,:,ss) = squeeze(nanmean(nanmean(tdata(auxChns,:,:,:)),4));
     allint(:,ss) = nanmean(nansum(nanmean(tint(auxChns,:,:,:)),3),4);
     allint2(:,:,ss) = squeeze(nanmean(nanmean(tint(auxChns,:,:,:)),4));
%     plot(-2400:4:796,nanmean(tdata(4,:,:),3))
blup = squeeze(nanmean(nanmean(tdata(auxChns,:,:,:)),4));
%  figure,plot(-2400:4:796,blup)
    ss = ss+1;
end
 figure,plot(-2400:4:796,squeeze(nanmean(alldata2,3)))
  hold on,plot(-2400:4:796,squeeze(nanmean(allint,2)))
%  hold on,plot(-2400:4:796,squeeze(nanmean(nansum(alldata2,3),4)))
  hold on,plot(-2400:4:796,squeeze(nanmean(alldata,2)))
  hold on,plot(-2400:4:796,squeeze(nanmean(allint,2))+squeeze(nanmean(alldata,2)))
  hold on, plot(-2400:4:796,squeeze(nanmean(nanmean(tint0(auxChns,:,:),3),1)))
  hold on, plot(-2400:4:796,squeeze(nanmean(nanmean(tint0(auxChns,:,:),3),1))'+squeeze(nanmean(allint,2))+squeeze(nanmean(alldata,2)))
 vline(0),hline(0)