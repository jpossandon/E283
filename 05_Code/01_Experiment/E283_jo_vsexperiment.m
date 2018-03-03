%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E275 - visual search task
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all                                                                 % we clear parameters?
% this is for debugging
win.DoDummyMode             = 0;                                            % (1) is for debugging without an eye-tracker, (0) is for running the experiment
win.stim_test               = 1;                                            % (1) for testing the stimulators (always when doing the experiment), (0) to skip
% PsychDebugWindowConfiguration(0,0.7);                                       % this is for debugging with a single screen

% Screen parameters
win.whichScreen             = 0;                                             % (CHANGE?) here we define the screen to use for the experiment, it depend on which computer we are using and how the screens are conected so it might need to be changed if the experiment starts in the wrong screen
win.FontSZ                  = 20;                                           % font size
win.bkgcolor                = 127;                                          % screen background color, 127 gray
win.Vdst                    = 66;                                           % (!CHANGE!) viewer's distance from screen [cm]         
win.res                     = [1920 1080];                                  %  horizontal x vertical resolution [pixels]
win.wdth                    = 51;                                           %  51X28.7 cms is teh size of Samsung Syncmaster P2370 in BPN lab EEG rechts
win.hght                    = 28.7;                                         % 
win.pixxdeg                 = win.res(1)/(2*180/pi*atan(win.wdth/2/win.Vdst));% 
win.trial_max_length        = 12;                                           % this is the max length of search,
win.center_thr              = 3*win.pixxdeg;                                % this is overriden below to set it to the center two columns
win.targ_thr                = win.pixxdeg;                                  % this is overriden below accoring to the spacing of the targets
win.tfix_target             = .5;

% Tactile stimulation settings (using the box stimulator)
win.tact_freq               = 200;                                          % frequency of stimulation in Hz
win.stim_dur                = .025;                                         % duration of tactile stimulation. The vibrator takes some time to stop its motion so for around 50 ms we use 25 ms of stimulation time (ask Tobias for the exact latencies they have measured)
win.stim_min_latency        = .325;                                         % minimum time from trial start (new image appearance) or previous stimulation for a tactile stimulation to occur
win.halflife                = win.trial_max_length/4;                       % we use an exponential distribution for a flat hazard function. Here the denominator set the duration in which half of the times will occur an stimulation
win.stim_lambda             = log(2)./win.halflife;                                 
win.rT_ratio                = 2;
% Blocks and trials
win.ncols                   = 8;
win.nrows                   = 6;
win.rep_item                = 12;
win.exp_trials              = win.ncols*win.nrows*win.rep_item;
win.test_trials             = 16;
win.t_perblock              = 24;
if mod(win.exp_trials,win.t_perblock),error('Number of trials per block do not match with total amount of trials'),end
win.calib_every             = 2; 
win.nBlocks                 = win.exp_trials/win.t_perblock;
win.catchp                  = .1;
% Audio, white noise parameters
win.wn_vol                  = .2;                                           % (CHANGE?) adjust to subject comfort

% Device input during the experiment
win.in_dev                  = 1;                                            % (1) - keyboard  (2) - mouse  (3) - pedal (?)    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENT START
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (!CHANGE!) adjust this to the appropiate screen
display(sprintf('\n\n\n\n\n\nPlease check that the display screen resolution is set to %dx%d and the screen borders are OK. \nIf screen resolution is changes matlab (and the terminal) need to be re-started.\n',win.res(1),win.res(2)))
if ismac                                                                    % this bit is just so I can run the experiment in my mac without a problem
    exp_path                = '/Users/jossando/trabajo/E275/';              % path in my mac
else
    exp_path                = 'C:\Users\bpn\Documents\jpossandon\E275\';
end

win.s_n                     = input('Subject number: ','s');                % subject id number, this number is used to open the randomization file
win.fnameEDF                = sprintf('s%02dvs.EDF',str2num(win.s_n));       % EDF name can be only 8 letters long, so we can have numbers only between 01 and 99
pathEDF                     = [exp_path 'data/' sprintf('s%02dvs/',str2num(win.s_n))];                           % where the EDF files are going to be saved
if exist([pathEDF win.fnameEDF],'file')                                         % checks whether there is a file with the same name
    rp = input(sprintf('!Filename %s already exist, do you want to overwrite it (y/n)?',win.fnameEDF),'s');
    if (strcmpi(rp,'n') || strcmpi(rp,'no'))
        error('filename already exist')
    end
else
    mkdir(pathEDF)
end

win.s_age                   = input('Subject age: ','s');
win.s_hand                  = input('Subject handedness for writing (l/r): ','s');
win.s_gender                = input('Subject gender (m/f): ','s');
setStr                      = sprintf('Subject %d\nAge %s\nHandedness %s\nGender %s\n',win.s_n,win.s_age,win.s_hand,win.s_gender); % setting summary
fprintf(setStr); 

AssertOpenGL();                                                             % check if Psychtoolbox is working (with OpenGL) TODO: is this needed?
ClockRandSeed();                                                            % this changes the random seed

[IsConnected, IsDummy] = EyelinkInit(win.DoDummyMode);                      % open the link with the eyetracker
assert(IsConnected==1, 'Failed to initialize EyeLink!')
Eyelink('command', '!*write_ioport 0x378 0');                               % set lpt to 0
                                                          % disable key listening by MATLAB windows(CTRL+C overridable)

prevVerbos = Screen('Preference','Verbosity', 2);                           % this two lines it to set how much we want the PTB to output in the command and display window 
prevVisDbg = Screen('Preference','VisualDebugLevel',3);                     % verbosity-1 (default 3); vdbg-2 (default 4)
Screen('Preference', 'SkipSyncTests', 2)                                    % for maximum accuracy and reliability
win.start_time = clock;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START PTB SCREEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[win.hndl, win.rect]        = Screen('OpenWindow',win.whichScreen,win.bkgcolor);   % starts PTB screen

[win.cntr(1), win.cntr(2)] = WindowCenter(win.hndl);                        % get where is the display screen center
Screen('BlendFunction',win.hndl, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     % enable alpha blending for smooth drawing
HideCursor(win.hndl);                                                       % this to hide the mouse
Screen('TextSize', win.hndl, win.FontSZ);                                   % sets teh font size of the text to be diplayed
KbName('UnifyKeyNames');                                                    % recommended, called again in EyelinkInitDefaults


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START AUDIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq                        = 48000;                                        % noise frequency in Hz
dur                         = 16;                                           % noise duration in seconds
wavedata1                   = rand(1,freq*dur)*2-1;                         % the noise
wavedata2                   = sin(2.*pi.*win.tact_freq.*[0:1/freq:dur-1/freq]);      % the stimulus
wavedata                    = [wavedata1 ; wavedata2];                      % double for stereo
nrchannels                  = 2;                                            % ditto
% 
 InitializePsychSound;                                                       % Perform basic initialization of the sound driver:
% 
 try                                                                         % Try with the 'freq'uency we wanted:
     pahandle = PsychPortAudio('Open', 9, [], 0, freq, nrchannels);
 catch                                                                       % Failed. Retry with default frequency as suggested by device:
     fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
     fprintf('Sound may sound a bit out of tune, ...\n\n');
     psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
 end
 PsychPortAudio('FillBuffer', pahandle, wavedata);                           % Fill the audio playback buffer with the audio data 'wavedata':
%  PsychPortAudio('Volume', pahandle,win.wn_vol);                              % Sets the volume (between 0 and 1)
 s = PsychPortAudio('GetStatus', pahandle);                                  % Status of the port, necessary later to defice wheter to start or stop audio


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTRUCTIONS IN GERMAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if win.in_dev == 1
    txtdev = ['Leertaste dr' 252 'cken (press Space)'];
elseif win.in_dev == 2
    txtdev = 'Maustate klicken';
end

% This is the safest way so umlaut are corretly displayed UTF. 
% code for umlauts 228(a),252(u)
txt1     = double(['Die Funktionsf' 228 'higkeit der Ger' 228 'te muss '...    
            252 'berpr' 252 'ft werden\n' txtdev ' zum Fortfahren..']);             
txt2     = double(['Stimulatoren werden auf den Handr' 252 'cken besfestigt ...\n' txtdev]);
txt3     = double(['Der linke Stimulator wird getestet (vibriert drei mal), ...\n danach ' txtdev]);
txt4     = double(['Der rechte Stimulator wird getestet (vibriert drei mal), ...\n danach ' txtdev]);
txt5     = double(['Beide Stimulatoren vibrieren gleichzeitig (drei mal), ...\n danach '  txtdev]); 
txt6     = double(['Beginn des Experiments \n Test Block \n Die H' 228 ...
            'nde bitte parallel positionieren (parallel). \n Zum Beginnen die ' txtdev]);
txt7     = double(['Beginn des Experiments\n Test Block \n  F' 252 'r den n' 228 ...
            'chsten Block bitte die H' 228 'nde ' 252 'berkreuzen (crossed). \n Zum Fortfahren die ' txtdev]);      
% txt 9&10 are contingent to the block so they are defined later
txt10    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde parallel positionieren. Bitte ignoriere die Vibration. Zum Fortfahren die ' txtdev]);
txt11    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde ' 252 'berkreuzen. Bitte ignoriere die Vibration. Zum Fortfahren die '  txtdev]);
txt12    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde parallel positionieren. Die Vibration gibt an, auf welcher Seite des Bildschirms das Target (Kreis ohne Strich) erscheint. Zum Fortfahren die ' txtdev]);
txt13    = double(['F' 252 'r den n' 228 'chsten Block die H' 228 ...
        'nde ' 252 'berkreuzen. Die Vibration gibt an, auf welcher Seite des Bildschirms das Target (Kreis ohne Strich) erscheint. Zum Fortfahren die '  txtdev]);
txt14    = double(['Das Experiment ist nun beendet. Vielen Dank f' 252 'r die Teilnahme! '  txtdev]);

%these are for debugging
handstr  = {'Left','Right','','','Left','Right','','','Left','Right','','','Left','Right'};
crossstr = {'Uncrossed','Crossed'};
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE-TRACKER SETUP, OPEN THE EDF FILE AND ADDS INFO TO ITS HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[win.el, win.elcmds] = setup_eyetracker( win, 1);                            % Setups the eye-tracker with the before mentioned parameters
 win.el.eye_used = 0;
DrawFormattedText(win.hndl,txt1,'center','center',255,55);                  % This 'draws' the text, nothing is displayed until Screen flip
Screen('Flip', win.hndl);                                                   % This is the command that changes the PTB display screen. Here it present the first instructions.

if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
    waitForKB_linux({'space'});                                             % press the space key in the keyboard
elseif win.in_dev == 2
    [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
end

OpenError = Eyelink('OpenFile', win.fnameEDF);                                  % opens the eye-tracking file. It can only be done after setting-up the eye-tracker 
if OpenError                                                                % error in case it is not possible, never happened that I know, but maybe if the small hard-drive aprtition of the eye=tracker is full
    error('EyeLink OpenFile failed (Error: %d)!', OpenError), 
end
Eyelink('Command', sprintf('add_file_preamble_text ''%s''', setStr));       % this adds the information about subject to the end of the header  
wrect = Screen('Rect', win.hndl);                                           
Eyelink('Message','DISPLAY_COORDS %d %d %d %d', 0, 0, wrect(1), wrect(2));  % write display resolution to EDF file

ListenChar(2)                                                               % disable MATLAB windows' keyboard listen (no unwanted edits)
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMULATOR TEST
% ALWAYS: left-hand is # 1 and right-hand is # 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if win.stim_test
    DrawFormattedText(win.hndl,txt2,'center','center',255,55);
    Screen('Flip', win.hndl);
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,txt3,'center','center',255,55);                  % TEST LEFT STIMULATOR (three times)
    Screen('Flip', win.hndl);
     t1              = PsychPortAudio('Start', pahandle, 0, 0, 0);               % Start the sounds
    WaitSecs(1);
    for e=1:3
        Eyelink('command', '!*write_ioport 0x378 0');                  
        WaitSecs(1+rand(1));      
        Eyelink('command', '!*write_ioport 0x378 1');                           % start stimulation by sending a signal through the parallel port (a number that was set by E275_define_tact_states)
        WaitSecs(win.stim_dur);                                                 % for the specified duration
        Eyelink('command', '!*write_ioport 0x378 0');                          % stop stimulation
        WaitSecs(win.stim_dur);
    end
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,txt4,'center','center',255,55);                  % TEST RIGHT STIMULATOR (three times)
    Screen('Flip', win.hndl);
    WaitSecs(1);
    for e=1:3
        Eyelink('command', '!*write_ioport 0x378 0');
        WaitSecs(1+rand(1));      
        Eyelink('command', '!*write_ioport 0x378 2');                           % start stimulation by sending a signal through the parallel port (a number that was set by E275_define_tact_states)
        WaitSecs(win.stim_dur);       
        Eyelink('command', '!*write_ioport 0x378 0');                          % stop stimulation
        WaitSecs(win.stim_dur);
    end
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end

    DrawFormattedText(win.hndl,txt5,'center','center',255,55);                  % TEST BOTH STIMULATORs (three times), this is not needed for this experiment but anyways
    Screen('Flip', win.hndl);
    WaitSecs(1);
    for e=1:3
        Eyelink('command', '!*write_ioport 0x378 0');
        WaitSecs(1+rand(1));      
        Eyelink('command', '!*write_ioport 0x378 3');                    
        WaitSecs(win.stim_dur);       
        Eyelink('command', '!*write_ioport 0x378 0');    
        WaitSecs(win.stim_dur);
    end
    Eyelink('command', '!*write_ioport 0x378 0');
    if win.in_dev == 1                                                          % Waiting for input according to decided device to continue
        waitForKB_linux({'space'});                                             % press the space key in the keyboard
    elseif win.in_dev == 2
        [clicks,x,y,whichButton] = GetClicks(win.hndl,0);                       % mouse clik
    end
 PsychPortAudio('Stop', pahandle);
    fprintf('\nTest ready!');
else
    fprintf('\nSkipping stimulators Test.\n Only during experiment debugging');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EYE-TRACKER CALIBRATION AND VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EyelinkDoTrackerSetup(win.el);                                              % calibration/validation (keyboard control)
Eyelink('WaitForModeReady', 500);

[image,map,alpha]   = imread([exp_path 'stimuli/blackongrt.jpg']);          % drift correction dot image
fixIndex            = Screen('MakeTexture', win.hndl, image);               % this is one of the way PTB deals with images, the image matrix is transformed in a texture with a handle that can be user later to draw the image in theb PRB screen


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Randomization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn                  = str2num(win.s_n);
win.vsTrials        = vsCreateTrials(win.trial_max_length,win.rep_item,win.ncols,win.nrows);
win.vsTrials(win.test_trials+1:end) = [];
experimentrials     = vsCreateTrials(win.trial_max_length,win.rep_item,win.ncols,win.nrows);                                    % (!!TODO: set-up number of trials)
for e = 1:length(experimentrials)
    win.vsTrials(e+win.test_trials) = experimentrials(e);
end
if rem(sn,2)                                                           % we balance across subjects on which position they start(according to their subject number) whether they start the experiment with the hand crossed or uncrossed, again this should not matter that much
    win.blockcond   = [zeros(1,win.test_trials),...                         % blockcond refer to the crossing (0-uncross; 1-cross)
        repmat([zeros(1,win.t_perblock),ones(1,win.t_perblock)],1,win.nBlocks/2)];
else
    win.blockcond   = [zeros(1,win.test_trials),...
        repmat([ones(1,win.t_perblock),zeros(1,win.t_perblock)],1,win.nBlocks/2)];
end

    win.blockcue    = [2*ones(1,win.test_trials),...                        % we balance the blocking of cue informativeness across subjces and first block hand position (I guess this is completly irrelevant anyways)
     ismember(sn,[1:4:200,2:4:200])*ones(1,win.exp_trials/2),...            % blockcue refer to the informativenes of the cue
     1-ismember(sn,[1:4:200,2:4:200])*ones(1,win.exp_trials/2)];            % 2 = test trials, 0 - uninformative, 1 - informative 

nTrials             = win.exp_trials+win.test_trials;
win.block_start     = [[1,zeros(1,win.test_trials-1)],[repmat([1,zeros(1,win.t_perblock-1)],1,win.nBlocks)]];                                         % (!!TODO: is this anymore necessary?)

%%% position of the two center column for center threhsold and target threshold
posVec              = vsCreateStimulus(win.vsTrials(1).stimulus.posGrid, ...
                    win, win.hndl, true,10,5);
Screen('FillRect', win.hndl, win.bkgcolor);  

win.center_thr      = ceil(max(unique(diff(unique(posVec.scr(1,:)))))+...
    .5*max((unique(diff(unique(posVec.scr(1,:)))))));        
win.targ_thr        = [ceil(.75*(max(unique(diff(unique(posVec.scr(1,:))))))) ...
    ceil(.75*(max(unique(diff(unique(posVec.scr(2,:)))))))];

win.result.perf     = zeros(1,nTrials);
win.result.rT       = zeros(1,nTrials);
win.result.stim     = zeros(1,nTrials);
            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE ACTUAL EXPERIMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b                   = 0;                                                    % block flag
for nT = 1:nTrials                                                          % loop throught the experiment trials
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BLOCK START
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if  win.block_start(nT) == 1                                            % if it is a trial that starts a block   
         if nT>1
              PsychPortAudio('Stop', pahandle);                                       % Stop the white noise after the last trial
         end
        % Hand position block intructions
        if nT ==1 && win.blockcond(nT) == 0                                  % practice trials and uncrossed
            draw_instructions_and_wait(txt6,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif nT ==1 && win.blockcond(nT) == 1                              % practice trials amd crossed, this according to randomization should be unnecesary
            draw_instructions_and_wait(txt7,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif win.blockcond(nT) == 0                % uncrossed
            txt8    = double(['Block ' num2str(b) '/' num2str(win.nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde parallel positionieren (parallel). \n Zum Fortfahren die ' txtdev]);
            draw_instructions_and_wait(txt8,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif win.blockcond(nT) == 1                % crossed
            txt9    = double(['Block ' num2str(b) '/' num2str(win.nBlocks+1) ' beendet \n Pause \n  F' 252 'r den n' 228 ... 
            'chsten Block bitte die H' 228 'nde ' 252 'berkreuzen (crossed). \n Zum Fortfahren die ' txtdev]);
            draw_instructions_and_wait(txt9,win.bkgcolor,win.hndl,win.in_dev,1)
        end      
        b = b+1;
        if nT>win.test_trials
            win.halflife(b)                = median(win.result.rT(find([ones(1,nT), zeros(1,nTrials-nT)])))/win.rT_ratio;    % half-life set to median of RT/1.5
            win.stim_lambda(b)             = log(2)./win.halflife(b); 
        end 
        if nT>1 %&& ismember(nT, win.t_perblock+win.test_trials+1:win.calib_every*win.t_perblock:nTrials)                              % we calibrate every two small blocks
            EyelinkDoTrackerSetup(win.el);
        end
        
        if win.blockcond(nT) == 0 && win.blockcue(nT)==0                                % this gives the instructions of crossing and stimulation relevance all again
            draw_instructions_and_wait(txt10,win.bkgcolor,win.hndl,win.in_dev,1)
%             DrawFormattedText(win.hndl, txt10,'center','center',255,55);
        elseif win.blockcond(nT) == 1 && win.blockcue(nT)==0
            draw_instructions_and_wait(txt11,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif win.blockcond(nT) == 0 && win.blockcue(nT)==1
            draw_instructions_and_wait(txt12,win.bkgcolor,win.hndl,win.in_dev,1)
        elseif win.blockcond(nT) == 1 && win.blockcue(nT)==1
           draw_instructions_and_wait(txt13,win.bkgcolor,win.hndl,win.in_dev,1)
        end
 
          t1          = PsychPortAudio('Start', pahandle, 0, 0, 0);
          tt = 1;
%         Screen('Flip', win.hndl);
%         if win.in_dev == 1                                                              
%             waitForKB_linux({'space'});                                           
%         elseif win.in_dev == 2
%             GetClicks(win.hndl,0);                                                      
%         end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TRIALS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Eyelink('message','TRIALID %d', nT);
    if win.blockcond(nT) == 0 
        Eyelink('Command','record_status_message ''Block Uncross %d/%d Trial %d/%d %d/%d''',b,win.nBlocks+1,tt,win.t_perblock,nT,win.exp_trials+win.test_trials);
    elseif win.blockcond(nT) == 1
        Eyelink('Command','record_status_message ''Block Cross %d/%d Trial %d/%d %d/%d''',b,win.nBlocks+1,tt,win.t_perblock,nT,win.exp_trials+win.test_trials);
    end
    Eyelink('WaitForModeReady', 50);
    tt=tt+1;
    Screen('FillRect', win.hndl, win.bkgcolor);                         % remove what was writte or displayed
    Screen('DrawTexture', win.hndl, fixIndex);                          % drift correction image
    Screen('Flip', win.hndl);
   
    EyelinkDoDriftCorrect2(win.el,win.res(1)/2,win.res(2)/2,0)          % drift correction 
    Screen('FillRect', win.hndl, win.bkgcolor);
 
        posVec      = vsCreateStimulus(win.vsTrials(nT).stimulus.posGrid, ...
                            win, win.hndl, true,10,5);
  
    Eyelink('StartRecording');  
    Eyelink('WaitForModeReady', 50);      
    
    targetpos   = posVec.scr(:,find(win.vsTrials(nT).stimulus.posGrid));
    Eyelink('message','METATR block %d',win.blockcond(nT));                 % block condition
    Eyelink('WaitForModeReady', 50);
    Eyelink('message','METATR cue %d',win.blockcue(nT));                    % cue condition
    Eyelink('WaitForModeReady', 50);
    Eyelink('message','METATR block_start %d',win.block_start(nT));         % if it was the first image in the block
    Eyelink('WaitForModeReady', 50);
    Eyelink('message','METATR tpos %d',win.vsTrials(nT).stimulus.tgtIndx);  % which target
%     
    Screen('Flip', win.hndl);                                               % actual image change, we message it to the eye-tracke and set the timer
    Eyelink('message','SYNCTIME');                                          % so trial zero time is just after image change
    Eyelink('command', '!*write_ioport 0x378 %d',96);                       % image appearance trigger, same as in my other free-viewing data
    tstart      = GetSecs;                                                  % timer trial start
%     last_stim   = tstart;                                                   % in this case, this mean image appearace
    Eyelink('WaitForModeReady', 50);
%    fname = sprintf('vsStim_Target%02d', win.vsTrials(nT).stimulus.tgtIndx)
% saveScreenshotAsPNG(win.hndl, fname)
    Eyelink('command', '!*write_ioport 0x378 %d',0);                        % flush the parallel port
%     Eyelink('WaitForModeReady', 50);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TACTILE STIMULATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % need area IN WHICH STIMULATION COULD OCCUUR
    % % of catch trials
    % timing of stimulation
    % condition of crossinf
    % condition of informativeness
    % information to the eye-ttracker about a tactile stimulus occuring
    
%     rvals     =  win.stim_min_latency + (win.stim_max_latency-win.stim_min_latency).*rand(100,1); % interstimulation intervals defined by an uniform distribution    
%     rvals       = win.stim_min_latency + ...                                % interstimulation intervals defined by an exponential distribution with half-life = 1/lambda
%                 (-1./win.stim_lambda(b) .* log(rand([1 100])));                % rvals = win.stim_min_latency + exprnd(1./win.stim_lambda,1,100);    by pass stat toolbox
    rvals = .150+rand(1)/10;
%     if rand(1)<win.catchp                                                   % catch trials in which stimulation occurs anywhere
        catcht = 1;
%     else
%         catcht = 0;
%     end
    stim_idx    = 1;
    % here we check when the subjects find the target and whether we
    % stimulate or not
    center      = 1; 
    prevst      = 0; fixt = 0;
    
    while GetSecs<tstart+win.trial_max_length                                       % lopp until trials finishes
         if fixt && GetSecs>targettime+win.tfix_target                              % if time runs off, then perf and rT = 0 
            win.result.perf(nT)     = 1;
            win.result.rT(nT)       = targettime-tstart; 
            break
        end
        if GetSecs>tstart+rvals(stim_idx) && nT>win.test_trials && ~prevst && ~fixt
            if  center || catcht                                                    % Stimulation infromative: 1 - left uncross; 2 - right uncross; 5 - left cross; 6 - right cross
                if win.blockcue(nT)==1                                              % when cue is instructive
                    if targetpos(1)<win.res(1)/2 && win.blockcond(nT) == 0          % left target / uncrossed / LH stimulation
                        stim = 1;
                    elseif targetpos(1)>win.res(1)/2 && win.blockcond(nT) == 0      % right target / uncrossed / RH stimulation
                        stim = 2;
                    elseif targetpos(1)<win.res(1)/2 && win.blockcond(nT) == 1      % left target / crossed / RH stimulation
                        stim = 6;
                    elseif targetpos(1)>win.res(1)/2 && win.blockcond(nT) == 1      % right target / crossed / LH stimulation
                        stim = 5;
                    end
                else                                                                % if stimulation is not instructive, we choose randomly
                    stim      = round(1+rand(1))+8;                                 % Stimulation uninformative: 9 - left uncross; 10 - right uncross; 13 - left cross; 14 - right cross      
                     if win.blockcond(nT) == 1
                        stim      = stim+4;
                    end
                end
                Eyelink('command', '!*write_ioport 0x378 %d',stim);               % start stimulation by sending a signal through the parallel port 
                WaitSecs(win.stim_dur);       
                prevst =1;                                                       % only one stim per trial
                Eyelink('command', '!*write_ioport 0x378 %d',0);                 % stop stimulation
    %           last_stim = GetSecs;
              %for testing purporses
                display(sprintf('\n Delivered stimulation at the %s hand %s  %4.2f',handstr{stim},crossstr{win.blockcond(nT)+1},GetSecs-tstart))
            end
        end
        % here we check if the target was found and whether the gaze is on
        % the central columns
        [data,type] = get_ETdata;
        if type ==6 % start fixation
            if abs(data.genx(win.el.eye_used+1)-win.res(1)./2)<win.center_thr
                center = 1;
            else
                center = 0;
            end
            if data.genx(win.el.eye_used+1)>targetpos(1)-win.targ_thr(1) && ...
                data.genx(win.el.eye_used+1)<targetpos(1)+win.targ_thr(1) && ...
                data.geny(win.el.eye_used+1)>targetpos(2)-win.targ_thr(2) && ...
                data.geny(win.el.eye_used+1)<targetpos(2)+win.targ_thr(2)
                if fixt == 0
                    fixt        = 1;
                    targettime  = GetSecs;                   % (!!TODO)maybe datatime?
                end
            else
                fixt = 0;
            end
        end
    end
    if prevst
        win.result.stim(nT)     = 1;
    end
    if ~fixt
        win.result.rT(nT)       = win.trial_max_length;
    end
        
    Eyelink('message','METATR rt %d',round(win.result.rT(nT)*1000));       % reaction time in the eye=tracker file
    Eyelink('WaitForModeReady', 50);
    Eyelink('message','METATR stim %d',prevst);       % reaction time in the eye=tracker file
    Eyelink('WaitForModeReady', 50);
    
    Screen('FillRect', win.hndl, win.bkgcolor);                             % remove what was writte or displayed
    Eyelink('command', '!*write_ioport 0x378 %d',0);                        % flush the parallel port for last time
    Eyelink('WaitForModeReady', 50);
    Eyelink('StopRecording');
    Eyelink('WaitForModeReady', 50);
    save([pathEDF,win.fnameEDF(1:end-3),'mat'],'win')
end
WaitSecs(1)
 win.end_time = clock;
 save([pathEDF,win.fnameEDF(1:end-3),'mat'],'win')
 Screen('FillRect', win.hndl, win.bkgcolor);                         % remove what was writte or displayed
 Screen('Flip', win.hndl);
 draw_instructions_and_wait(txt14,win.bkgcolor,win.hndl,win.in_dev,1)
 
 PsychPortAudio('Stop', pahandle);                                       % Stop the white noise after the last trial
                                          

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finishing EDF file and transfering info (in case experiment is interrupted
% this can be run to save the eye-tracking data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Task Iteration done; save files, restore stuff, DON'T clear vars %%%%
Eyelink('CloseFile');
Eyelink('WaitForModeReady', 500); % make sure mode switching is ok
if ~win.DoDummyMode
    % get EDF->DispPC: file size [bytes] if OK; 0 if cancelled; <0 if error
    rcvStat = Eyelink('ReceiveFile', win.fnameEDF, pathEDF,1);
    if rcvStat > 0 % only sensible if real connect
        fprintf('EDF received to %s (%.1f MiB).\n',pathEDF,rcvStat/1024^2);
    else
        fprintf(2,'EDF file reception error: %d.\n', rcvStat);
    end
end

% save([pathEDF,win.fnameEDF(1:end-3),'mat'],'win')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSING ALL DEVICES, PORTS, ETC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Most of the commands below are important for being able to restart the
% experiment, in the case the experiment crushs they should be entered
% manually

% A known issue: Eyelink('Shutdown') crashing Matlab in 64-bit Linux
% cf. http://tech.groups.yahoo.com/group/psychtoolbox/message/12732
% not anymore it seems
%if ~IsLinux(true), Eyelink('Shutdown'); end
PsychPortAudio('Close', pahandle);                                          % close the audio device
Eyelink('Shutdown');                                                        % close the link to the eye-tracker

Screen('CloseAll');                                                         % close the PTB screen
Screen('Preference','Verbosity', prevVerbos);                               % restore previous verbosity
Screen('Preference','VisualDebugLevel', prevVisDbg);                        % restore prev vis dbg
ListenChar(1)                                                               % restore MATLAB keyboard listening (on command window)
