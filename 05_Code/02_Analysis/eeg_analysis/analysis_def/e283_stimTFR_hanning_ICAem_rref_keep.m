p.subj                  = p.subj;
p.times                 = [1000 1500];
p.npermute              = 1;
p.data                  = 'stimTFR';
p.rref                  = 'yes';
p.keep                  = 'yes';
p.plot                  = 1;
p.analysis_type         = {'ICAem'};

p.interval              = [-.15 .6 .025]; % [start end step]
% p.colorlim              = [-10 10];

p.model                 = 1;
% p.cov                   = '';
% p.model_cov             = '[]';
p.cov                   = 'cross_inf';   
p.model_cov             = 'cat(2,covs.covariate_cross,covs.covariate_Inf)'; 
p.interact              = {[1 2];[1 3];[2 3];[1 2 3]};
p.coeff                 = {'const','LR','Cross','Inf','LRxCross','LRxInf','CrossxInf','LRxCrossxInf'};

p.mirror                = [];

p.analysisname          = [p.data '_' p.analysis_type{1} '_rref' p.rref '_cov_' p.cov];   

p.cfgTFR.channel            = 'all';	
p.cfgTFR.keeptrials         = 'yes';	                
p.cfgTFR.method             = 'mtmconvol';
p.cfgTFR.taper              = 'hanning';
% p.cfgTFR.width              = 5; 
p.cfgTFR.output             = 'pow';	
p.cfgTFR.foi                = 4:1:35;
p.cfgTFR.winsize            = 0.300;
p.cfgTFR.t_ftimwin          = p.cfgTFR.winsize*ones(1,length(p.cfgTFR.foi));
p.cfgTFR.tapsmofrq          = ones(1,length(p.cfgTFR.foi));
p.cfgTFR.toi                = (-p.times(1)+p.cfgTFR.winsize*1000/2+10:10:-10+p.times(2)-p.cfgTFR.winsize*1000/2)/1000;	
p.cfgTFR.pad                = 2*sum(p.times)/1000;    
p.bsl                       = [-(p.times(1)/1000)+p.cfgTFR.winsize/2 -p.cfgTFR.winsize/2];
%