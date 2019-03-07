%%
% add splines
splinestocomb = {'fixpre_pxdiff_','fixpre_pydiff_'};
comblab       = 'fixpre_pxdiffpydiffcomb_';
fstspl       = strmatch(splinestocomb{1},resultSplines.coeffs);
scndspl       = strmatch(splinestocomb{2},resultSplines.coeffs);

splvals = cellfun(@(x) x{1},regexp(resultSplines.coeffs,'[^_]+$','match'), 'UniformOutput', false);

ii =1;
for fst = fstspl'
    for snd = scndspl'
        resultSplinescomb.B(:,ii,:,:) = resultSplines.B(:,fst,:,:)+resultSplines.B(:,snd,:,:);
        resultSplinescomb.coeffs{ii,1} = [comblab splvals{fst} '_' splvals{snd}];
        ii = ii+1;
    end
end
resultSplinescomb.times = resultSplines.times;
