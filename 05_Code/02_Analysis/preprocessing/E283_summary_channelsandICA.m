%%
E283_params
subjcorpath = '/Users/jossando/trabajo/E283/07_Analysis/02_Preprocessing/subjects_master_files';
icapath     = '/Users/jossando/trabajo/E283/07_Analysis/02_Preprocessing/ICAm';
chanelim    = [];
icaelim     = [];
for ss = p.subj
    filename = sprintf('%s/S%02dVS_channels_corrections.mat',subjcorpath,ss);
    if exist(filename)
       load(filename)
       if isempty(chan_cor.elim_chan)
        chanelim    = [chanelim;0];   
       else
        chanelim    = [chanelim;length(chan_cor.elim_chan)];
       end
    else
        chanelim    = [chanelim;0];  
    end
    filename = sprintf('%s/s%02dvs/s%02dvs_ICA.mat',icapath,ss,ss);
    load(filename)
    icaelim = [icaelim;length(union(cfg_ica.comptoremove,cfg_ica.comptoremove_m))];
end

fprintf('max number chan removed: %d\n',max(chanelim))
fprintf('mean number chan removed: %1.2f\n',mean(chanelim))
fprintf('0/%d,(chanremov/#subj)\n1/%d,(chanremov/#subj)\n2/%d,(chanremov/#subj)\n3/%d,(chanremov/#subj)\n4/%d,(chanremov/#subj)\n',sum(chanelim==0),sum(chanelim==1),sum(chanelim==2),sum(chanelim==3),sum(chanelim==3))

fprintf('\n\naverage number ICAcomp removed: %1.2f(SD:%1.2f)\n',mean(icaelim),std(icaelim))
