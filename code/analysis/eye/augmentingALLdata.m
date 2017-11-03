%%
% Adds symbol and target information to alleyedata structure 
%
clear
pathexp = '/Users/jossando/trabajo/E283/';
load(fullfile(pathexp,'analysis','eyedata','alleyedata.mat'))
load(fullfile(pathexp,'analysis','eyedata','tgtPos.mat'))

% remove subject with a different search array
data = struct_elim(data,data.subject<6,2,1);

%%
% assign to each fixation an element and to every saccade an starting and
% end element, and a distance to the target

data.elfix      = nan(1,length(data.posx));
data.elfixini   = nan(1,length(data.posx));
data.elfixend   = nan(1,length(data.posx));
data.DToTarg    = nan(1,length(data.posx));
data.onTarg     = zeros(1,length(data.posx));
data.nextToTarg = nan(3,length(data.posx));

for ev = 1:length(data.posx)
   
    % where in the array is each fixation
    if ~isnan(data.posx(ev))
        ixH         = find(data.posx(ev)>posVec.hLims,1,'last');
        ixV         = find(data.posy(ev)>posVec.vLims,1,'last');
        if ~isempty(ixH) && ~isempty(ixV) 
            if ixH>0 && ixH<9 && ixV>0 && ixV<7
                data.elfix(ev) = sub2ind([6 8],ixV,ixH);
            end
            [I,J]            = ind2sub([6 8],data.tpos(ev));
        
            % gaze on the target position or inmediatly next
            if ixH-J==0 && ixV-I==0              % on target
                data.onTarg(ev) = 1;
            elseif abs(ixH-J)==1 && abs(ixV-I)==0 % horizontally next 
                data.nextToTarg(1,ev) = 1;
            elseif abs(ixH-J)==0 && abs(ixV-I)==1 % vertically next  
                data.nextToTarg(2,ev) = 1;
            elseif abs(ixH-J)==1 && abs(ixV-I)==1 % diagonally next  
                data.nextToTarg(3,ev) = 1;
            end
        end
        data.DToTarg(ev) = sqrt((data.posx(ev)-posVec.scr(1,data.tpos(ev))).^2+... % euclidian distance
            (data.posy(ev)-posVec.scr(2,data.tpos(ev))).^2);
        
    end
    
    % where in the array is each saccade
    if ~isnan(data.posinix(ev))
        ixHini      = find(data.posinix(ev)>posVec.hLims,1,'last');
        ixVini      = find(data.posiniy(ev)>posVec.vLims,1,'last');
        ixHend      = find(data.posendx(ev)>posVec.hLims,1,'last');
        ixVend      = find(data.posendy(ev)>posVec.vLims,1,'last');
        if ~isempty(ixHini) && ~isempty(ixVini) && ~isempty(ixHend) && ~isempty(ixVend) 
            if ixHini>0 && ixHini<9 && ixVini>0 && ixVini<7
                data.elfixini(ev) = sub2ind([6 8],ixVini,ixHini);
            end
            if ixHend>0 && ixHend<9 && ixVend>0 && ixVend<7
                data.elfixend(ev) = sub2ind([6 8],ixVend,ixHend);
            end
        end
    end
    
end
       
%%
% add performance, repeated and revisited fixation
load(fullfile(pathexp,'analysis','eyedata','allRT.mat'))
result         = struct_elim(result,result.subject<6,2,1);
save(fullfile(pathexp,'analysis','eyedata','allRTFULL.mat'),'result')

data.tCorrect  = zeros(1,length(data.posx));
data.orderPreT = nan(1,length(data.posx));
data.refix     = zeros(1,length(data.posx));
data.revisit   = zeros(1,length(data.posx));
data.torevisit = zeros(1,length(data.posx));
data.torefix   = zeros(1,length(data.posx));
for ss = unique(data.subject)
    ss
    for tt = unique(data.trial(data.subject==ss))
        data.tCorrect(data.subject==ss & data.trial==tt) = result.perf(find(result.subject==ss & result.trial==tt));
        if result.perf(find(result.subject==ss & result.trial==tt))
            data.orderPreT(data.subject==ss & data.trial==tt & data.type==1) = fliplr(0:sum(data.subject==ss & data.trial==tt & data.type==1)-1);
            data.orderPreT(data.subject==ss & data.trial==tt & data.type==2) = fliplr(0:sum(data.subject==ss & data.trial==tt & data.type==2)-1);
        end
        auxindx = find(data.subject==ss & data.trial==tt & data.type==1);
        for ft = 2:length(auxindx)
           if ~isnan(data.elfix(auxindx(ft)))
               if data.elfix(auxindx(ft)) == data.elfix(auxindx(ft-1))
                   data.refix(auxindx(ft)) = 1;
                   data.torefix(auxindx(ft-1)) = 1;
               end
               if ft>3
                   [revisits,lag] = ismember(data.elfix(auxindx(ft)),data.elfix(auxindx(1:ft-2)));
                   if any(revisits)
                        if ~(sum(revisits)==1 && revisits(end)==1 && data.elfix(auxindx(ft)) == data.elfix(auxindx(ft-1))) % to avoid clasifying a fix as revisit when there has been thre refix in a row
                            data.revisit(auxindx(ft)) = ft-lag;
                            data.torevisit(auxindx(lag)) = 1;
                        end
                   end
               end
           end
        end
    end
end
save(fullfile(pathexp,'analysis','eyedata','alleyedataFULL.mat'),'data')
%%
% Get the relevant info into subjects eyedata files

for tk = unique(data.subject); % subject number

    if ismac    
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk),'expfolder','/Users/jossando/trabajo/E283/'); % this is just to being able to do analysis at work and with my laptop
    else
        cfg_eeg             = eeg_etParams_E283('sujid',sprintf('s%02dvs',tk));
    end

    cfg_eeg.filename                = sprintf('s%02dvs',tk);
    load([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'])
                         
    eyedata.events.elfix        = data.elfix(data.subject==tk);
    eyedata.events.DToTarg      = data.DToTarg(data.subject==tk);
    eyedata.events.onTarg       = data.onTarg(data.subject==tk);
    eyedata.events.nextToTarg   = data.nextToTarg(:,data.subject==tk);
    eyedata.events.tCorrect     = data.tCorrect(data.subject==tk);
    eyedata.events.orderPreT    = data.orderPreT(data.subject==tk);
    eyedata.events.refix        = data.refix(data.subject==tk);
    eyedata.events.revisit      = data.revisit(data.subject==tk);
    save([cfg_eeg.eyeanalysisfolder cfg_eeg.filename 'eye.mat'],'eyedata')
end
%%
% ss=1;
% figure,hold on
% for s=unique(data.subject)
%     for preT=0:20
%         durPreT(ss,preT+1) = mean(data.dur(data.subject==s & data.type==1 & data.orderPreT==preT));
%         
%     end
%     plot(1:20,durPreT(ss,2:end),'o','Color',[.7 .7 .7])
%     ss = ss+1;
% end
% errorbar(1:20,mean(durPreT(:,2:end)),std(durPreT(:,2:end),1,1)./sqrt(ss-1),std(durPreT(:,2:end),1,1)./sqrt(ss-1),'r')
% 
% %%
%         
% [N,BIN] = histc(data.DToTarg(data.type==1),50:100:1450);
% fixdurs = data.dur(data.type==1);
% for nb = 1:max(BIN)
%     mdur(nb) = median(fixdurs(BIN==nb));
% end
