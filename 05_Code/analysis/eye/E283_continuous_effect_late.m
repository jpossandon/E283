% analysis of effect of late stimulation in viewing behavior
% TODO (solved): covariate does not seem to work? Yes it does but it is a calue
% between 0 and 1
%%
clear 
if ismac
    path = '/Users/jossando/trabajo/E283/';
else
    path = '/net/store/users/jossando/E283/';
end
load([path 'analysis/eyedata/alleyedata.mat'])
load([path 'analysis/eyedata/allsampledata.mat'])
eyedata.events          = data;
eyedata.events.angle    = data.angle*pi/180;
subjects                = unique(eyedata.events.subject);

%%
E283_params
subjects            = p.subj;
horres              = 1920;
dt                  = 1500;                        % length of analysis in seconds;
et_sr               = 500;                         % sampling rate eye-tracker

cross               = [0 0 0 1 1 1 0 0 0 1 1 1];  % (no cross-no stim / no cross/left / no cross/right / no cross/both 
info                = [1 1 1 1 1 1 0 0 0 0 0 0];
% ??stimside            = [0 1 2 0 3 4];  %    cross-no stim  /  cross/left   /  cross/right   / cross/both)
 saclat              = cell(length(subjects),12);
 fixlat              = cell(length(subjects),12);
 fixelap              = cell(length(subjects),12);
conds               = {'u_no_inf','u_l_inf','u_r_inf','c_no_inf','c_l_inf','c_r_inf',...
    'u_no_uninf','u_l_uninf','u_r_uninf','c_no_uninf','c_l_uninf','c_r_uninf'};
% clear auxdata
earlylatlimit       = [105 200]; 
for s = 1:length(subjects)
     display(sprintf('Processing subject %d',subjects(s)))
     stimside            = [0 1 2 0 5 6 0 9 10 0 13 14]; 
  
    subjsample          = struct_select(sample,{'subject'},{['==' num2str(subjects(s))]},2);    % data from the current subject s
    
    indxstim            = find(stim.subject==subjects(s) & stim.trial>16);                      % stimulation for subjects s and after first 10 test trials
    val                 = stim.value(indxstim);                             % respective values for trigger, time and trial
    rnval               = length(val);
    t                   = stim.time(indxstim);
    tr                  = stim.trial(indxstim);  
    
    % control values, we take them for trials without stimulation after
    % image apearance (~ half trials) at the same time that stimulation
    % ocurred in the other trials
%     nostimtrl           = setdiff(17:max(stim.trial(stim.subject==subjects(s))),tr);
%     val                 = [val,zeros(1,rnval)];
%     tr                  = [tr,randsample(nostimtrl,rnval,'true')];
%     t                   = [t,randsample(t,rnval)];
        
    nreprog =1;                 % counter to separate trials in which the first saccade after stimulation ocurred after earlylatlimit
    for n = 1:length(t)
        
        crossaux                    = eyedata.events.block(find(eyedata.events.subject==subjects(s)  & eyedata.events.trial==tr(n),1,'first'));  % what is the crossing condition for this trial
        infoaux                    = eyedata.events.cue(find(eyedata.events.subject==subjects(s)  & eyedata.events.trial==tr(n),1,'first'));  % what is the cueing condition for this trial
        if isempty(crossaux)
            if ismember(val(n),[1,2,9,10])
                crossaux            = 0;
            else
                crossaux            = 1;
            end
            sprintf('No data for trial %d subject %d',tr(n),subjects(s))
        end
         if isempty(infoaux)
            if ismember(val(n),[1,2,5,6])
                infoaux            = 1;
            else
                infoaux            = 0;
            end
            sprintf('No data for trial %d subject %d',tr(n),subjects(s))
        end
            
        indxsac                     = find(eyedata.events.subject==subjects(s) & eyedata.events.trial==tr(n) & eyedata.events.type==2 & eyedata.events.start>t(n));  % the saccades for this specific subject, trial that ocurred after stimulation
        indxfix                     = find(eyedata.events.subject==subjects(s) & eyedata.events.trial==tr(n) & eyedata.events.type==1 & eyedata.events.start<t(n) & eyedata.events.end>t(n));  % the fixation for this specific subject, trial that started before and end after stimulation
        if ~isempty(indxsac)        % if there is a saccade
            auxlat                                          = eyedata.events.start(indxsac(1))-t(n); % latency from stimulation
            saclat{s,cross==crossaux & stimside == val(n) & info==infoaux}  = [saclat{s,cross==crossaux & stimside == val(n) & info==infoaux}; auxlat]; % get it for later
        end
        if ~isempty(indxfix)        % if there is a fixation
            auxlat                                          = eyedata.events.end(indxfix(1))-t(n); % stimulation to end of fixation
            fixlat{s,cross==crossaux & stimside == val(n) & info==infoaux}  = [fixlat{s,cross==crossaux & stimside == val(n) & info==infoaux}; auxlat]; % get it for later
            auxlat                                          = t(n)-eyedata.events.start(indxfix(1)); % start fixation to stimulation
            fixelap{s,cross==crossaux & stimside == val(n) & info==infoaux}  = [fixelap{s,cross==crossaux & stimside == val(n) & info==infoaux}; auxlat]; % get it for later
        end
        timsample                   = subjsample.time>t(n) & subjsample.time<=t(n)+dt;              % we look into 1 second after stimulation TODO: remove data that is after 5 sec ??
        timsamplebsl                = subjsample.time>t(n)-50 & subjsample.time<=t(n);                % and compare to baseline
        
        aux                         = subjsample.pos(1,subjsample.trial==tr(n) & timsample);           % get respective x pos
        fix_pre_posx                = -horres/2+round(nanmean(subjsample.pos(1,subjsample.trial==tr(n) & timsamplebsl))); % mean position in the baseline period
       
        if fix_pre_posx>5000 , fix_pre_posx = NaN;end,    
        
        auxdata(n,1,1:dt*et_sr/1000)            = NaN;                                                           % next three lines is to fill with NANs for trials with less data than 1.5 second
        auxdata(n,1,1:length(aux))              = aux;
        auxdata(n,1,auxdata(n,1,:)>horres | auxdata(n,1,:)<0) = NaN;
        
%         auxbsl                          = subjsample.pos(1,subjsample.trial==tr(n) & timsamplebsl);   % this is before stimulation so it is not necesarry to fill with NaNs
%         auxbsl(auxbsl>horres | auxbsl<0)  = NaN;
%         auxdata(n,1,:)                  = auxdata(n,1,:)-nanmean(auxbsl);                             % baseline correction
%         auxdata(n,1,:)                  = auxdata(n,1,:);%-nanmean(auxbsl);                             % baseline correction, but we do not baseline correct when we are usin prestim pos as a covariate
        
        
        if ~isempty(indxsac)            % data only for trials with saccades after stimulation and that occurred after earlylatlimit
            if auxlat>earlylatlimit(1) & auxlat<earlylatlimit(2)
                auxdata_reprog(nreprog,1,:)              = auxdata(n,1,:);
            end
        end
        
        datac(s).(conds{cross==crossaux & stimside==val(n) & info==infoaux})(n,:)   = auxdata(n,1,1:dt*et_sr/1000);
            
        % the next lines generete the dummy coding for
        % left,right,bilateral, against cros-nostim
        % this does not work whene tere is stimulation always
%         if any(val==0)
%             XY(n,1) = crossaux;  % crossing
%             XY(n,2) = infoaux;  % cue informativens
% 
% 
%             if ismember(val(n),[1,5,9,13])
%                 XY(n,3) = 1;
%             elseif ismember(val(n),[2,6,10,14])
%                 XY(n,4) = 1;
%             end
%             XY(n,10) = fix_pre_posx;
%         end

%         if ~isempty(indxsac)          % data only for trials with saccades after stimulation and that occurred after earlylatlimit
%             if auxlat>earlylatlimit(1) & auxlat<earlylatlimit(2) 
%                 datac_reprog(s).(conds{cross==crossaux & stimside==val(n)})(nreprog,:)   = auxdata(n,1,1:499);   
%             
%                 XY_reprog(nreprog,1) = crossaux;
%                 if val(n)>0 && val(n)<3
%                     XY_reprog(nreprog,val(n)+1) = 1;
%                 elseif val(n)>2 && val(n)<5   
%                     XY_reprog(nreprog,val(n)-1) = 1;
%                 end
%                  XY_reprog(nreprog,6) = fix_pre_posx;
%                 nreprog = nreprog+1;
%             end
%         end
        % effect coding without comparison to no-stim
        if n<rnval+1
            if crossaux
                XY_eff(n,1)     = 1;
            else
                XY_eff(n,1)     = -1;
            end
            if infoaux
                XY_eff(n,2)     = 1;
            else
                XY_eff(n,2)     = -1;
            end
          
            if ismember(val(n),[1,5,9,13]) 
                XY_eff(n,3) = -1;
            elseif ismember(val(n),[2,6,10,14])
                XY_eff(n,3) = 1;
            end
            XY_eff(n,8)     = fix_pre_posx;
        end
    end
       
    condsF = fields(datac);
    for st=1:length(condsF)
        datac(s).(condsF{st})            = nanmean(datac(s).(condsF{st}));
%         datac_reprog(s).(conds{st})     = nanmean(datac_reprog(s).(conds{st}));
    end
    %         if any(val==0)
%     XY(:,5) = XY(:,1).*XY(:,2);
%     XY(:,6) = XY(:,1).*XY(:,3);
%     XY(:,7) = XY(:,1).*XY(:,4);
%     XY(:,8) = XY(:,2).*XY(:,3);
%     XY(:,9) = XY(:,2).*XY(:,4);
% end
    XY_eff(:,4) = XY_eff(:,1).*XY_eff(:,2);
    XY_eff(:,5) = XY_eff(:,1).*XY_eff(:,3);
    XY_eff(:,6) = XY_eff(:,2).*XY_eff(:,3);
    XY_eff(:,7) = XY_eff(:,1).*XY_eff(:,2).*XY_eff(:,3);
%     XY_reprog(:,4) = XY_reprog(:,1).*XY_reprog(:,2);
%     XY_reprog(:,5) = XY_reprog(:,1).*XY_reprog(:,3);
%  
    elec.channeighbstructmat = 0;
%     [datac(s).B,datac(s).Bt,datac(s).STATS,datac(s).T] = regntcfe(auxdata,XY,1,'dummy',elec,0);
    [datac_eff(s).B,datac_eff(s).Bt,datac_eff(s).STATS,datac_eff(s).T] = regntcfe(auxdata(1:rnval,:,:),XY_eff,1,'effect',elec,0);
    
%     [datac_reprog(s).B,datac_reprog(s).Bt,datac_reprog(s).STATS,datac_reprog(s).T] = regntcfe(auxdata_reprog,XY_reprog,1,'dummy',elec,0);
    clear XY auxdata XY_reprog auxdata_reprog XY_eff
%     save([path 'analysis/model_late'],'datac','datac_eff','saclat','fixlat','fixelap','subjects')
     save([path 'analysis/model_late'],'datac_eff','saclat','fixlat','fixelap','subjects')
end

%%
% effect coding
cmap = colormap('lines');
allb = squeeze([datac_eff.B]);
numb = size(datac_eff(1).B,2);
nsuj = length(datac_eff);
figure,
hold on
for e = 1:size(allb,1)/nsuj
% plot(allb(e:numb:end,:)','Color',cmap(e,:))
   h(e) = plot(mean(allb(e:numb:end,:))','Color',cmap(e,:),'LineWidt',3);
end
legend(h,{'constant','Cross','Info','LR','CrossxInfo','CrossxLR','InfoxLR','CrossxInfoxLR','Covariate'})
set(gca,'XTick',0:250:750,'XTickLabel',[0:250:750]/500,'FontSize',14)
axis([0 750 -200 200])
view([90 90])
ylabel('Gaze Position (pix)','FontSize',16)
xlabel('Time (s)','FontSize',16)
%%
% mean vs mediam!
% cmap = colormap('lines');
% allb = squeeze([datac.B]);
% figure,
% hold on
% for e = 1:7
% % plot(allb(e:5:end,:)','Color',cmap(e,:))
% h2(e) = plot(mean(allb(e:7:end,:))','Color',cmap(e,:),'LineWidt',3);
% end
% legend(h2,{'constant','Cross','Left','Right','CrossxL','CrossxR','PosCov'})
% set(gca,'XTick',0:100:700,'XTickLabel',0:1:7)
% axis([0 750 -100 100])
% view([90 90])
%%
%2nd level
% clear 
% if ismac
%     path = '/Users/jossando/trabajo/touch/analysis/';
% else
%     path = '/net/store/users/jossando/touch/analysis/';
% end
% 
% rpgroup  = {'all','reprog'};
% for rp = 1
%     if rp == 1
%         load([path 'model_late'],'datac','subjects')
%         allb        = reshape([datac.B],[9,57,750]);
%     elseif rp == 2
%         load([path 'model_late'],'datac_reprog','subjects')
%         allb        = reshape([datac_reprog.B],[9,57,750]);
%     end
%         allb(1,:,:) = allb(1,:,:)-640.*ones(ones,size(allb,2),size(allb,3)); %to test if baseline is different from midline
%     for g = 4%[1,4,5]
%         if g == 1
%             subjs = 1:57;  % 1 to 22 are the first group
%             sname    = 'allsubjects_late';
%         elseif g == 2
%             subjs = 1:22;  % 1 to 22 are the first group
%             sname    = 'firstgroup_late';
%         elseif g == 3
%             subjs = 23:40;  % 1 to 22 are the first group
%             sname    = 'secondgroup_late';
%         elseif g == 4
%             subjs = 1:47;  % 1 to 22 are the first group
%             sname    = 'allright_late';
%         elseif g == 5
%             subjs    = 48:57;  % 1 to 22 are the first group
%             sname       = 'allleft_late';   
%         elseif g == 6
%              subjs1  = 1:47;  % 1 to 22 are the first group
%              subjs2  = 48:57;  % 1 to 22 are the first group
%             
%             sname       = 'indep_late'; 
%         end
% 
%         tiempos                     = 1:2:1500;
%         data2nd(1,:,:,:)            = permute(allb(:,subjs,:),[1 3 2]);
%         
%         for methods = {'signpermT'}%'bootet','bootsym','boottrimsym','boottrimet'}
%             if g<6
%                 data2nd(1,:,:,:)            = permute(allb(:,subjs,:),[1 3 2]);
%                 elec.channeighbstructmat = 0;
%                 [result]                    = regmodel2ndstat(data2nd,tiempos,elec,10000,methods{:},'tfce');
%                 result.subjects             = subjects;
% 
%             else
%                 data2nd1(1,:,:,:)           = permute(allb(:,subjs1,:),[1 3 2]);
%                 data2nd2(1,:,:,:)           = permute(allb(:,subjs2,:),[1 3 2]);
% 
%                 elec.channeighbstructmat = 0;
%                 [result]                    = regmodel2ndstat_ind(data2nd1,data2nd2,tiempos,elec,1000,methods{:},'tfce');
% 
%                 result.subjects1             = subjects(subjs1);
%                 result.subjects2             = subjects(subjs2);  
%             end
%                 result.name                 = sname;
%                 result.times                = tiempos;
%                 result.subjects             = subjects(subjs);
% 
%                 result.allb                 = allb(:,subjs,:);
% 
%                 save([path 'model2nd_' sname '_' rpgroup{rp} '_' methods{:}],'result')
%                 clear data2nd
%                 
%        end
%     end
% end
% 
