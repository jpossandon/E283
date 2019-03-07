clear
pathexp = '/Users/jossando/trabajo/E283/';
patheye = fullfile(pathexp,'07_Analysis','03_Eye');
group   = 2;    % 1-  subject 6 to 27, 2, 28 to 72, 3 all
load(fullfile(patheye,'eyedata','alleyedataFULL.mat'),'data')    % info from subjects EDF files, augments with augmentinfALLdata.m 
load(fullfile(patheye,'eyedata','allRTFULL.mat'))                    % info from subjects matlab files 
load(fullfile(patheye,'eyedata','tgtPos.mat'))                   % position of targets on the screen

if group==1
    data    = struct_elim(data,data.subject>27,2,1);
    result  = struct_elim(result,result.subject>27,2,1);
    namegr  = 'randStim';
elseif group==2
    data    = struct_elim(data,data.subject<28,2,1);
    result  = struct_elim(result,result.subject<28,2,1);
    namegr  = 'earlyStim';
else
    namegr  = 'all';
end
subjects    = unique(data.subject);
Ns          = length(subjects);
subjectswin = unique(result.subject);
if ~all(subjects==subjectswin)
    error('subjects in data and result files do not match')
end
clear subjectswin
setAbsoluteFigureSize

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RT AND PERFORMANCE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SUMMARY

s=1;
SUMMARY.rTBINS   = 0:.5:12;
for ss = subjects
    SUMMARY.subjIndx(s) = ss;
    auxres  = struct_select(result,{'subject'},{['==' num2str(ss)]},2);
     if length(auxres.perf) ~= 592
        display(sprintf('\nSubject %d has only %d trials',ss,length(auxres.perf))) %  subject 20 and 22 the EEG recording was started late
     end
    SUMMARY.nT(s)         = length(auxres.perf);
    SUMMARY.perf(s)       = sum(auxres.perf)./length(auxres.perf);
    SUMMARY.rT(s)         = mean(auxres.rT(auxres.perf==1));
    SUMMARY.rTpbin(s,:)   = histc(auxres.rT(auxres.perf==1),SUMMARY.rTBINS);
    
    % rt per position
    for pos = 1:48
        SUMMARY.rtpPos(s,pos) = mean(auxres.rT(auxres.perf==1 & auxres.tpos==pos));
    end
    %nfixs
   
    if group==2
        % value: 1 - LUI, 2 - RUI, 5 - LCI, 6 - RCI, 9 - LUunI, 10 - RUuni
        %  13 - LCunI, 14 - RCuni
        % cue: 0/1 - informative, 1/2 - uninformative ,2 test trials
        SUMMARY.rTpCond(s,:) = [mean(auxres.rT(auxres.perf==1 & auxres.value==1)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==2)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==5)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==6)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==9)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==10)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==13)),...
            mean(auxres.rT(auxres.perf==1 & auxres.value==14))];
    end
     s = s+1;
end
SUMMARY.condLabels  = {'LUI','RUI','LCI','RCI','LUunI','RUunI','LCunI','RCunI'};
SUMMARY.condValues = [1,2,5,6,9,10,13,14];
SUMMARY.infoLabels  = {'Info','Info','Info','Info','unInfo','unInfo','unInfo','unInfo'};
SUMMARY.meanRT      = mean(SUMMARY.rT);
SUMMARY.medianRT    = median(SUMMARY.rT);
SUMMARY.sdRT        = std(SUMMARY.rT);
SUMMARY.meanRTse    = SUMMARY.sdRT./sqrt(Ns);

SUMMARY.meanRTpCond     = mean(SUMMARY.rTpCond);
SUMMARY.medianRTpCond   = median(SUMMARY.rTpCond);
SUMMARY.sdRTpCond       = std(SUMMARY.rTpCond);
SUMMARY.seRTCond        = SUMMARY.sdRTpCond./sqrt(Ns);