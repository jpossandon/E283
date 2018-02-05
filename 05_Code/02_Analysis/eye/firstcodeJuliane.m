% loads the data
load('allsampledata.mat')  % sample data
load('alleyedata.mat')     % event and stim data

%%
% First an example of exploration
% TODO: here we should plot the image of the search
% indxtrials = stim.trial(stim.trial>20 & stim.value==2 & stim.subject ==11); % find all trials with a specific condition and subject
% 
% % plot every one of this trials
% for t=indxtrials  % loop across trials
%     indxsamples = find(sample.trial==t & sample.subject==11);  % find indexes for the specifcit trials
%     auxsamples  = sample.pos(:,indxsamples);  % find the position of the sample (1st rox is x and 2nd row is y)
%     figure  % create a figure window
%     plot(auxsamples(1,:),auxsamples(2,:),'.') % plots it, only dots and not lines
%     axis ij % changes the axxis to 0,0 is top,left
%     axis([0 1920 0 1080]) % set the border of x aixs between 0 and 1920 and y between 0 and 1080
%     box off  
% end

%%
% all fixation on the different conditions
trigger = [1,2,5,6,9,10,13,14]; % different trigger
leg     = {'LUI','RUI','LCI','RCI','LUU','RUU','LCU','RCU'}; % labels of the condition
% loop through triggers
for trig = 1:length(trigger)
xpos = []; % create variables for the x and y positions of fications
ypos = [];
for subj = unique(stim.subject) % loop across subject
    indxtrials = stim.trial(stim.trial>16 & stim.value==trigger(trig) & stim.subject ==subj); % find indexes for a specific triger and subject
    for t=indxtrials
        xpos = [xpos,data.posx(data.type==1 & data.subject==subj & data.trial ==t)]; % aand get the data for fixation, for the specific subject and trial
        ypos = [ypos,data.posy(data.type==1 & data.subject==subj & data.trial ==t)];
    end
end

% plot
figure
plot(xpos,ypos,'.k')
axis([0 1920 0 1080])
axis ij
box off
vline(1920/2)  % this look for in mathworks
hline(1080/2)
title(leg{trig})
end

%%
% create the file with all rT data
load('alleyedata.mat')
subject = unique(stim.subject);

path = '/Users/jossando/trabajo/E283/data/';

for s = subject
    load(sprintf('%ss%02dvs/s%02dvs.mat',path,s,s))
    win.result.subject = s.*ones(1,length(win.result.rT));
    win.result.value   = zeros(1,length(win.result.rT));
    win.result.value(stim.trial(stim.subject==s)) = stim.value(stim.subject==s);
    win.result.cue     = zeros(1,length(win.result.rT));
    win.result.cue(win.result.value>0 & win.result.value<7) = 1; % informative
    win.result.cue(win.result.value>7) = 2; % uninformative
    result = struct_up('result',win.result,2);
end
save('/Users/jossando/trabajo/E283/analysis/eyedata/allRT','result')
%%
% plot subject reaction times different condition
% you can try to do it per side
load('allRT.mat')
cond = [0,1,2]; % 0 no cue, 1 cue nformative, 2 cue uninformative
subjects = unique(result.subject); % which subject do we have
RTresult = [];
for s=1:length(subjects)
    for c=1:length(cond)
        RTresult(s,c) = mean(result.rT(result.perf==1 & result.subject==subjects(s) & result.cue==cond(c))); % this get the mean RT for succesful trials, for every subject condition
    end
end

%%
% plot the rt result
figure,
hold on % so we can plot several things in the same figure
for e=1:3
    plot(e,RTresult(:,e),'.k','MarkerSize',12) % this plot every subject data per condition e
end
line(repmat([1;2;3],1,size(RTresult,1)),RTresult','Color',[.8 .8 .8])
plot(1:3,median(RTresult),'sb','MarkerSize',14) % this plot the median value
errorbar(1:3,mean(RTresult),std(RTresult)/sqrt(s),'LineWidth',2,'Color',[1 0 0]) % errorbar with SEM
axis([0 4 0 8]) % limits of the plot
ylabel('Reaction Time (s)','FontSize',20) % labels uh
xlabel('Condition','FontSize',20)
set(gca,'XTick',1:3,'XTickLabel',{'No Cue','Instruct.','No Instruct.'},'FontSize',16) % set further parameters of the axes

%%
clear
load('allRT.mat')
load('alleyedata.mat');
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
save('/Users/jossando/trabajo/E283/analysis/eyedata/allRT','result')