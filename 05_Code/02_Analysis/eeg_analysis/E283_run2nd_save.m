load(cfg_eeg.chanfile)
if ~isempty(stimBspl)
    % this one was for 1d splines
       % this one works for 2d not sure if for 1d
    if any(strfind(model,'2D'))
     coeffs_spl             = strcat([ufpredict.param.event]','_',...
                             strrep(strrep(strrep({ufpredict.param.name},':','XX'),'(',''),')','')','_',...
                             regexprep(mat2cell(num2str(cell2mat({ufpredict.param.value}')),ones(1,length(ufpredict.param)),17), '\s+', '_'));
    else
      coeffs_spl      = strcat([ufpredict.param.event]','_',strrep(strrep(strrep({ufpredict.param.name},':','XX'),'(',''),')','')','_',num2str([ufpredict.param.value]'));
    end  
resultSplines.B       = stimBspl;
    resultSplines.coeffs  = coeffs_spl;
    resultSplines.splines = ufpredict.unfold.splines;
    resultSplines.times   = ufpredict.times;
    
     stimBdiffSpl = [];
%     sdiff = 1;
%     splnames = cellfun(@(x) x{1},regexp(resultSplines.coeffs,'^([.*_]?.*_)','match'), 'UniformOutput', false);
%     for splType=unique(splnames)'
%         splb        = strmatch(splType,splnames);
%         for splVal = 1:floor(length(splb)/2)
%             stimBdiffSpl(:,sdiff,:,:) = resultSplines.B(:,splb(end-splVal+1),:,:)-resultSplines.B(:,splb(splVal),:,:);
%             Spldiff_coeffs{sdiff} = [resultSplines.coeffs{splb(end-splVal+1)} '_minus_' resultSplines.coeffs{splb(splVal)}];
%             sdiff = sdiff+1;
%         end
%     end
%     resultSpldiff        = regmodel2ndstat(stimBdiffSpl,resultSplines.times,elec,2000,stattype,mc);
%     resultSpldiff.coeffs = Spldiff_coeffs; 
end

result        = regmodel2ndstat(stimB,unfold.times,elec,np,stattype,mc);
coeffs_nospl  = strcat([unfold.param.event]','_',strrep(strrep(strrep({unfold.param.name},':','XX'),'(',''),')','')');
result.coeffs = coeffs_nospl;
mkdir(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm'))


if ~isempty(bslcor) 
    if exist('freqbands','var')
        glmname = ['glmALLbslcorr_' freqbands{fb} stattype '_' mc];
    else
        glmname = ['glmALLbslcorr_' stattype '_' mc];
    end
    if ~isempty(stimBspl)
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',glmname),'result','resultSplines','cfg_eeg')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',glmname),'result','cfg_eeg')
    end
    if ~isempty(stimBdiffSpl)
         save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',[glmname '_spldiff']),'resultSpldiff','resultSplines','cfg_eeg')
    end
else
    if exist('freqbands','var')
        glmname = ['glmALL_' freqbands{fb} stattype '_' mc];
    else
        glmname = ['glmALL_' stattype '_' mc];
    end
    if ~isempty(stimBspl)
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',glmname),'result','resultSplines','cfg_eeg')
    else
        save(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',glmname),'result','cfg_eeg')
    end
end%


fileID = fopen(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',[glmname '.txt']),'w');

fprintf(fileID,'%s %s\n%s\n',p.analysisname,model,result.statlabel);
for cb = 1:length(result.clusters)
      fprintf(fileID,'\n%s\n',result.coeffs{cb});
    if ~isempty(result.clusters(cb).texto)
        for ll = 1:size(result.clusters(cb).texto,1)
            fprintf(fileID,'%s\n',result.clusters(cb).texto(ll,:));
        end
    else
        fprintf(fileID,'*no significant cluster*\n',result.coeffs{cb});
    end
end
fclose(fileID);

fileID = fopen(fullfile(cfg_eeg.eeganalysisfolder,p.analysisname ,model,'glm',[glmname '_spldiff.txt']),'w');

fprintf(fileID,'%s %s\n%s\n',p.analysisname,model,resultSpldiff.statlabel);
for cb = 1:length(resultSpldiff.clusters)
      fprintf(fileID,'\n%s\n',resultSpldiff.coeffs{cb});
    if ~isempty(resultSpldiff.clusters(cb).texto)
        for ll = 1:size(resultSpldiff.clusters(cb).texto,1)
            fprintf(fileID,'%s\n',resultSpldiff.clusters(cb).texto(ll,:));
        end
    else
        fprintf(fileID,'*no significant cluster*\n',resultSpldiff.coeffs{cb});
    end
end
fclose(fileID);