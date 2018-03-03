 function saveScreenshotAsPNG(winH, fname) %#ok<DEFNU>
        if nargin < 2, fname = ['saveScreen' datestr(now(),30)]; end
        [ign, fname] = fileparts(fname); %#ok<ASGLU> % strip file extension
        imgarray = Screen('GetImage', winH);
        savepath = [fullfile(pwd(),fname) '.png'];
        fprintf('Saving window content to "%s" ... ', savepath)
        imwrite(imgarray, savepath, 'png')
        fprintf('DONE.\n')