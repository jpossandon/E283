[trls.LU_I,events]              = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==1'},p.times);            
[trls.RU_I,events]              = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==2'},p.times);
[trls.LC_I,events]              = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==5'},p.times);            
[trls.RC_I,events]              = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==6'},p.times);            
[trls.LU_unI,events]            = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==9'},p.times);            
[trls.RU_unI,events]            = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==10'},p.times);
[trls.LC_unI,events]            = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==13'},p.times);            
[trls.RC_unI,events]            = define_event(cfg_eeg,eyedata,'ETtrigger',{'value','==14'},p.times);            

trls.left                       = [trls.LU_I;trls.LC_I;trls.LU_unI;trls.LC_unI];
trls.right                     = [trls.RU_I;trls.RC_I;trls.RU_unI;trls.RC_unI];
% trls.all                       = [trls.left;trls.right];
p.trls_stim                        = {'LU_I','LC_I','RU_I','RC_I','LU_unI','LC_unI','RU_unI','RC_unI','left','right'};
p.trls_glm                     = '{trls.left;trls.right}'; 

% covs.covariate_cross                = {{[-1*ones(1,size(trls.LU_I,1)),ones(1,size(trls.LC_I,1)),-1*ones(1,size(trls.LU_unI,1)),ones(1,size(trls.LC_unI,1))]};...
%                                             {[-1*ones(1,size(trls.RU_I,1)),ones(1,size(trls.RC_I,1)),-1*ones(1,size(trls.RU_unI,1)),ones(1,size(trls.RC_unI,1))]}};
%    
% covs.covariate_Inf                = {{[ones(1,size(trls.LU_I,1)),ones(1,size(trls.LC_I,1)),-1.*ones(1,size(trls.LU_unI,1)),-1.*ones(1,size(trls.LC_unI,1))]};...
%                                             {[ones(1,size(trls.RU_I,1)),ones(1,size(trls.RC_I,1)),-1*ones(1,size(trls.RU_unI,1)),-1.*ones(1,size(trls.RC_unI,1))]}};
% p.covariates                   = eval(p.model_cov);
p.trls                         = eval(p.trls_glm);