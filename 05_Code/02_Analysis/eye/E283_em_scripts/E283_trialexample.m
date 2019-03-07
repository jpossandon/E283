load(fullfile(patheye,'eyedata','allsampledata.mat')) 
%%

for examp = 1:20
    ss = randsample(data.subject,1);
    tt = randsample(data.trial(data.subject==ss),1);

    subjindx  = find(data.subject==ss);
    targetpos = unique(data.tpos(data.subject==ss & data.trial == tt));
    indxtrial = find(sample.subject==ss & sample.trial ==tt);
    cmap1 = cbrewer('qual','Set1',9);
    fh = figure;
    set(fh, 'InvertHardCopy', 'off')
    fh.Position = [100 100 posVec.stimRes/2];
    axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
    axis off
    % ha = gca;
    % ha.Color = [127 127 127]./255;
    patch([0 posVec.stimRes(1) posVec.stimRes(1) 0],[0 0 posVec.stimRes(2) posVec.stimRes(2)],[127 127 127]./255)
    for el = 1:size(posVec.scr,2)
        rectangle('Position',[posVec.scr(1,el)-posVec.pixxdeg*.75/2 posVec.scr(2,el)-posVec.pixxdeg*.75/2 posVec.pixxdeg*.75 posVec.pixxdeg*.75],...
            'Curvature',[1 1],'EdgeColor',[1 1 1],'LineWidth',.2)
        if el~=targetpos
           line([posVec.scr(1,el) posVec.scr(1,el)],[posVec.scr(2,el)-posVec.pixxdeg*.75/2-10 posVec.scr(2,el)-posVec.pixxdeg*.75/2+10],'Color',[1 1 1],'LineWidth',.2)
        hold on  
        end
    end
%         line([posVec.scr(1,targetpos) posVec.scr(1,targetpos)],[posVec.scr(2,targetpos)-posVec.pixxdeg*.75/2-10 posVec.scr(2,targetpos)-posVec.pixxdeg*.75/2+10],'Color',[1 1 1],'LineWidth',.2)
%         hold on
   
    if unique(data.cue(data.subject==ss & data.trial == tt))==0
        othertrials = unique(data.trial(find(data.subject==ss & data.tpos==targetpos & data.cue==1 & data.value==unique(data.value(data.subject==ss & data.trial == tt))-8)));
        colindx = [1 2];
    else
        othertrials = unique(data.trial(find(data.subject==ss & data.tpos==targetpos & data.cue==0 & data.value==unique(data.value(data.subject==ss & data.trial == tt))+8)));
        colindx = [2 1];
    end
     plot(sample.pos(1,indxtrial),sample.pos(2,indxtrial),'.-','MarkerSize',4,'Color',[0 0 0],'MarkerEdgeColor',cmap1(colindx(1),:),'LineWidth',.2)

    if ~isempty(othertrials)
        ot = othertrials(randsample(length(othertrials),1));
        indxtrial = find(sample.subject==ss & sample.trial ==ot);
%         plot(sample.pos(1,indxtrial),sample.pos(2,indxtrial),'.-','MarkerSize',4,'Color',[0 0 0],'MarkerEdgeColor',cmap1(colindx(2),:),'LineWidth',.2)
    end
        axis off
    figsize = [17.2/2 17.6/2*fh.Position(4)/fh.Position(3)];
    doimage(gcf,fullfile(patheye,'figures','trial_examples'),'png',sprintf('Subj_%d_trials_%d_%d',ss,tt,ot),'600','opengl',figsize,1)
end