%%
load('/Users/jossando/trabajo/E283/06_RawData/s06vs/s06vs.mat')
vsPosGrid       = win.vsTrials(1).stimulus.posGrid;
[grdY, grdX]    = size(vsPosGrid);
posVec.borderH  = 10;
posVec.borderV  = 5;
posVec.stimRes  = win.res;  
posVec.borderH  = posVec.borderH/100*posVec.stimRes(1);
posVec.borderV  = posVec.borderV/100*posVec.stimRes(2);
posVec.dsX      = (posVec.stimRes(1)-2*posVec.borderH)/(grdX+1); % one more distance than dots | - o - o - o - |
posVec.dsY      = (posVec.stimRes(2)-2*posVec.borderV)/(grdY+1);
[vrtPosX, vrtPosY] = meshgrid(linspace(posVec.dsX+posVec.borderH, posVec.stimRes(1)-posVec.dsX-posVec.borderH, grdX),...
                             linspace(posVec.dsY+posVec.borderV, posVec.stimRes(2)-posVec.dsY-posVec.borderV, grdY));
posVec.vrt      = [vrtPosX(:) vrtPosY(:)]'; 
scrPosX         = floor( (win.res(1)-posVec.stimRes(1))/2 + vrtPosX ); % round to pixels
scrPosY         = floor( (win.res(2)-posVec.stimRes(2))/2 + vrtPosY ); % round to pixels
posVec.scr      = [scrPosX(:) scrPosY(:)]'; % 2-row vector, for PTB on DisplayPC
posVec.wdth     = win.wdth;
posVec.hght     = win.hght;
posVec.target_thr   = win.targ_thr;
posVec.pixxdeg  = posVec.stimRes(1)/(2*(180/pi)*atan(posVec.wdth/2/66));

[mins] = min(posVec.scr');
posVec.hLims = [round(mins(1)-posVec.dsX/2),round(unique(posVec.scr(1,:))+posVec.dsX/2)];
posVec.vLims = [round(mins(2)-posVec.dsY/2),round(unique(posVec.scr(2,:))+posVec.dsY/2)];
save('/Users/jossando/trabajo/E283/07_Analysis/03_Eye/eyedata/tgtPos','posVec')