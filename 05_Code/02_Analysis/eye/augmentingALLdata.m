%%
%%
% create the file with all rT data
clear
pathexp = '/Users/jossando/trabajo/E283/';
% create the file with all rT data
load(fullfile(pathexp,'07_Analysis','03_Eye','eyedata','alleyedata.mat'))
load(fullfile(pathexp,'07_Analysis','03_Eye','eyedata','tgtPos.mat'))
data    = struct_elim(data,data.subject<6,2,1);
subject = unique(data.subject);

path = '/Users/jossando/trabajo/E283/06_RawData/';
% path2 = '/Volumes/nibaldo/trabajo/E283/data/';
for s = subject
%     mkdir(sprintf('%ss%02dvs',path,s))
%     copyfile(sprintf('%ss%02dvs/s%02dvs.mat',path2,s,s),sprintf('%ss%02dvs/s%02dvs.mat',path,s,s))
    load(sprintf('%ss%02dvs/s%02dvs.mat',path,s,s))
    win.result.subject = s.*ones(1,length(win.result.rT));
    win.result.value   = zeros(1,length(win.result.rT));
    if s==54  % this is a fix for subjects that have changes edf and result file becqause the EEG recording started late, trial are then not with the same numbering
        win.result.value(stim.trial(stim.subject==s)-17) = stim.value(stim.subject==s); 
    elseif s==66
        win.result.value(stim.trial(stim.subject==s)-3) = stim.value(stim.subject==s); 
    else
        win.result.value(stim.trial(stim.subject==s)) = stim.value(stim.subject==s);
    end
    win.result.cue     = zeros(1,length(win.result.rT));
    win.result.cue(win.result.value>0 & win.result.value<7) = 1; % informative
    win.result.cue(win.result.value>7) = 2; % uninformative
    result = struct_up('result',win.result,2);
end
save('/Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/allRT','result')

%%
% infromation of tpos trial, and stimtime in result rt structure
subj    = unique(result.subject); % which subject do we have
aux20   = find(result.subject==20);
result  =  struct_elim(result,aux20(1:26),2,0);
aux22   = find(result.subject==22);
result  =  struct_elim(result,aux22(1:7),2,0);
    
result.tpos     = zeros(1,length(result.rT));
result.trial    = zeros(1,length(result.rT));
for s = 1:length(subj)
    auxindx     = find(data.subject==subj(s));
    auxindxres  = find(result.subject==subj(s));
    auxindstim  = find(stim.subject==subj(s));
    if ~isequal(sum(result.subject==subj(s)),length(unique(data.trial(auxindx))))
        error('No same amount of trial')
    end
    
    auxtops     = [];
    auxtrial    = [];
    auxtimestim = [];

    nn=1;
    auxpos  = data.tpos(auxindx);
    auxtime = stim.time(auxindstim);
    
    for nt = min(data.trial(auxindx)):max(data.trial(auxindx))
        auxt    = find(data.trial(auxindx)==nt);
        auxr    = unique(auxpos(auxt));
        if length(auxr)>1
            error('no single tpos value')
        end
        auxtops(nn)  = auxr;
        auxtrial(nn) = nt;
        
        auxt    = find(stim.trial(auxindstim)==nt);
        if ~isempty(auxt)
            auxtt   = unique(auxtime(auxt));
            if length(auxtt)>1
                error('no single time value')
            end
            auxtimestim(nn) = auxtt;
        else
            auxtimestim(nn) = NaN;
        end
        nn =nn+1;
    end
    
    result.tpos(auxindxres)     = auxtops;
    result.trial(auxindxres)    = auxtrial;
    result.stimtime(auxindxres) = auxtimestim;
end
save('//Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/allRT','result')

%%
% assign to each fixation an element and to every saccade an starting and
% end element, and a distance to the target

data.elfix      = nan(1,length(data.posx));
data.elfixini   = nan(1,length(data.posx));
data.elfixend   = nan(1,length(data.posx));
data.DToTarg    = nan(1,length(data.posx));
data.xToTarg    = nan(1,length(data.posx));
data.yToTarg    = nan(1,length(data.posx));
data.tToTarg    = nan(1,length(data.posx));
data.onTarg     = zeros(1,length(data.posx));
data.nextToTargH = nan(1,length(data.posx));
data.nextToTargV = nan(1,length(data.posx));
data.nextToTargD = nan(1,length(data.posx));
tic
for ev = 1:length(data.posx)
   
    % where in the array is each fixation
    if ~isnan(data.posx(ev))
        ixH         = find(data.posx(ev)>posVec.hLims,1,'last');
        ixV         = find(data.posy(ev)>posVec.vLims,1,'last');
        tFixTarg    = data.start(find(data.trial == data.trial(ev) & data.subject == data.subject(ev) & data.type==1,1,'last'));
        if ~isempty(ixH) && ~isempty(ixV) 
            if ixH>0 && ixH<9 && ixV>0 && ixV<7
                data.elfix(ev) = sub2ind([6 8],ixV,ixH);
            end
            [I,J]            = ind2sub([6 8],data.tpos(ev));
        
            % gaze on the target position or inmediatly next
            if ixH-J==0 && ixV-I==0              % on target
                data.onTarg(ev) = 1;
            elseif abs(ixH-J)==1 && abs(ixV-I)==0 % horizontally next 
                data.nextToTargH(ev) = 1;
            elseif abs(ixH-J)==0 && abs(ixV-I)==1 % vertically next  
                data.nextToTargV(ev) = 1;
            elseif abs(ixH-J)==1 && abs(ixV-I)==1 % diagonally next  
                data.nextToTargD(ev) = 1;
            end
                   
        end
        data.DToTarg(ev) = sqrt((data.posx(ev)-posVec.scr(1,data.tpos(ev))).^2+... % euclidian distance
            (data.posy(ev)-posVec.scr(2,data.tpos(ev))).^2);
        data.xToTarg(ev) = posVec.scr(1,data.tpos(ev))-data.posx(ev);
        data.yToTarg(ev) = posVec.scr(2,data.tpos(ev))-data.posy(ev);
        data.tToTarg(ev) = data.start(ev)-tFixTarg;
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
% sequential movements
data.seqmovH = nan(1,length(data.posx));
data.seqmovV = nan(1,length(data.posx));
data.seqmovD = nan(1,length(data.posx));
data.seqmovRF = nan(1,length(data.posx));
data.idiff    = nan(1,length(data.posx));
data.jdiff    = nan(1,length(data.posx));
for ev = 1:length(data.trial)-2
    I1 = []; J1= [];I2 = []; J2= [];
    if data.type(ev)==1 
        if data.trial(ev)==data.trial(ev+2) &&  ...
                data.type(ev+2)==1 && ~isnan(data.elfix(ev)) && ~isnan(data.elfix(ev+2))
            [I1,J1]            = ind2sub([6 8],data.elfix(ev));
            [I2,J2]            = ind2sub([6 8],data.elfix(ev+2));

            
        end
    end
    if data.type(ev)==2
        if ~isnan(data.elfixini(ev)) && ~isnan(data.elfixend(ev))
            [I1,J1]            = ind2sub([6 8],data.elfixini(ev));
            [I2,J2]            = ind2sub([6 8],data.elfixend(ev));
        end
    end
    if ~isempty(I2) && ~isempty(I1) && ~isempty(J2) && ~isempty(J1) 
        idiff   = I2-I1;
        jdiff   = J2-J1;
        data.seqmovH(ev)   = idiff == 0 & abs(jdiff) == 1;                         % movement to the next horizontal symbol
        data.seqmovV(ev)   = abs(idiff) == 1 & jdiff == 0;                         % movement to the next vertical symbol
        data.seqmovD(ev)   = abs(idiff) == 1 & abs(jdiff)== 1;
        data.seqmovRF(ev)  = idiff == 0 & jdiff == 0; 
        data.idiff(ev) = idiff; 
        data.jdiff(ev) = jdiff;
    end
    
end
toc
%%
% add performance, type of cue repeated and revisited fixation
load(fullfile(pathexp,'07_Analysis','03_Eye','eyedata','allRT.mat'))
result         = struct_elim(result,result.subject<6,2,1);
save(fullfile(pathexp,'07_Analysis','03_Eye','eyedata','allRTFULL.mat'),'result')

data.tCorrect  = zeros(1,length(data.posx));
data.value     = zeros(1,length(data.posx));
data.stimtime  = nan(1,length(data.posx));
data.orderPreT = nan(1,length(data.posx));
data.refix     = zeros(1,length(data.posx));
data.revisit   = zeros(1,length(data.posx));
data.torevisit = zeros(1,length(data.posx));
data.torefix   = zeros(1,length(data.posx));
data.latposStim = nan(1,length(data.posx));
data.orderposStim = nan(1,length(data.posx));
tic
for ss = unique(data.subject)
    ss
    for tt = unique(data.trial(data.subject==ss))
        data.tCorrect(data.subject==ss & data.trial==tt)    = result.perf(find(result.subject==ss & result.trial==tt));
        data.value(data.subject==ss & data.trial==tt)       = result.value(find(result.subject==ss & result.trial==tt));
        data.stim(data.subject==ss & data.trial==tt)        = result.stim(find(result.subject==ss & result.trial==tt));
        data.stimtime(data.subject==ss & data.trial==tt)    = result.stimtime(find(result.subject==ss & result.trial==tt));
        data.latposStim(data.subject==ss & data.trial==tt)  = data.start(data.subject==ss & data.trial==tt)-data.stimtime(data.subject==ss & data.trial==tt);
        indxFirstPost   = data.latposStim(data.subject==ss & data.trial==tt)>0;
        orderposStim    = nan(1,sum(data.subject==ss & data.trial==tt));
        typeposStim     = data.type(data.subject==ss & data.trial==tt);
        orderposStim(find(typeposStim==2 & indxFirstPost)) = 1:sum(typeposStim==2 & indxFirstPost);
        orderposStim(find(typeposStim==1 & indxFirstPost)) = 1:sum(typeposStim==1 & indxFirstPost);
        data.orderposStim(data.subject==ss & data.trial==tt) = orderposStim;
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
                            data.torevisit(auxindx(lag)) = ft-lag;
                        end
                   end
               end
           end
        end
    end
end
toc
save(fullfile(pathexp,'07_Analysis','03_Eye','eyedata','alleyedataFULL.mat'),'data')
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
    eyedata.events.nextToTargH  = data.nextToTargH(data.subject==tk);
    eyedata.events.nextToTargV  = data.nextToTargV(data.subject==tk);
    eyedata.events.nextToTargD  = data.nextToTargD(data.subject==tk);
    eyedata.events.seqmovH      = data.seqmovH(data.subject==tk);
    eyedata.events.seqmovV      = data.seqmovV(data.subject==tk);
    eyedata.events.seqmovD      = data.seqmovD(data.subject==tk);
    eyedata.events.seqmovRF     = data.seqmovRF(data.subject==tk);
    eyedata.events.tCorrect     = data.tCorrect(data.subject==tk);
    eyedata.events.orderPreT    = data.orderPreT(data.subject==tk);
    eyedata.events.refix        = data.refix(data.subject==tk);
    eyedata.events.revisit      = data.revisit(data.subject==tk);
    eyedata.events.xToTarg      = data.xToTarg(data.subject==tk);
    eyedata.events.yToTarg      = data.yToTarg(data.subject==tk);
    eyedata.events.tToTarg      = data.tToTarg(data.subject==tk);
    eyedata.events.latposStim   = data.latposStim(data.subject==tk);
    eyedata.events.orderposStim = data.orderposStim(data.subject==tk);
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
