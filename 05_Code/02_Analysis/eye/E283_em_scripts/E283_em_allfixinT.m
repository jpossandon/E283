%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple check that eyedata correspond with target positions and border
% between elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, hold on
set(gcf,'Position',[0 0 1280 700])
axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
plot(data.posx(data.type==1 & data.start>0),data.posy(data.type==1 & data.start>0),'.k','MarkerSize',3)
plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',16)
axis ij
vline(posVec.hLims)
hline(posVec.vLims,'r:')
tightfig
doimage(gcf,fullfile(patheye,'figures'),...
            'pdf',['allFixs_' namegr],[],1)
        
% this is for the different types of cue
cuevalues = [1 2 5 6;9 10 13 14];
infLabel  = {'Info','unInfo'};
for ii = 1:2

    figure, hold on
    set(gcf,'Position',[0 0 1280 700])
    for sp = [1,3,2,4]
        subplot(2,2,sp), hold on
        axis([0 posVec.stimRes(1) 0 posVec.stimRes(2)])
        plot(data.posx(data.type==1 & data.value==cuevalues(ii,sp) & data.start>0),...
            data.posy(data.type==1 & data.value==cuevalues(ii,sp) & data.start>0),'.k','MarkerSize',3)
        plot(posVec.scr(1,:),posVec.scr(2,:),'o','MarkerSize',16)
        axis ij
        vline(posVec.hLims)
        hline(posVec.vLims,'r:')
        title(SUMMARY.condLabels{sp+(ii-1)*4})
        
    end
    tightfig
    doimage(gcf,fullfile(patheye,'figures'),...
                    'pdf',['allFixs_' namegr '_' infLabel{ii}],[],1)
end