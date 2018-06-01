%%
PsychDebugWindowConfiguration(0,0.2);  
win.DoDummyMode             = 0;       
 win.whichScreen             = 0;        % (CHANGE?) here we define the screen to use for the experiment, it depend on which computer we are using and how the screens are conected so it might need to be changed if the experiment starts in the wrong screen
win.bkgcolor                = 127;   
win.tact_freq               = 200;                                          % frequency of stimulation in Hz
win.Vdst                    = 66;                                           % (!CHANGE!) viewer's distance from screen [cm]         
win.res                     = [1920 1080];                                  %  horizontal x vertical resolution [pixels]
win.wdth                    = 51;                                           %  51X28.7 cms is teh size of Samsung Syncmaster P2370 in BPN lab EEG rechts
win.hght                    = 28.7;                                         % 

[win.hndl, win.rect]        = Screen('OpenWindow',win.whichScreen,win.bkgcolor);   % starts PTB screen
[IsConnected, IsDummy] = EyelinkInit(win.DoDummyMode);                      % open the link with the eyetracker



%%
freq                        = 48000;                                        % noise frequency in Hz
dur                         = 15;                                           % noise duration in seconds
wavedata1                   = rand(1,freq*dur)*2-1;                         % the noise
wavedata2                   = sin(2.*pi.*win.tact_freq.*[0:1/freq:dur-1/freq]);      % the stimulus
wavedata                    = [wavedata1 ; wavedata2];     
nrchannels                  = 2;                                            % ditto

InitializePsychSound;                                                       % Perform basic initialization of the sound driver:

try                                                                         % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
catch                                                                       % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
    fprintf('Sound may sound a bit out of tune, ...\n\n');
    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
end
PsychPortAudio('FillBuffer', pahandle, wavedata);                           % Fill the audio playback buffer with the audio data 'wavedata':
 s = PsychPortAudio('GetStatus', pahandle);                                  % Status of the port, necessary later to defice wheter to start or stop audio

%%
[win.el, win.elcmds] = setup_eyetracker( win, 1);                            % Setups the eye-tracker with the before mentioned parameters
win.el.eye_used = 0;

OpenError = Eyelink('OpenFile', 'ERASEme.EDF'); 
%%
EyelinkDoTrackerSetup(win.el);                                              % calibration/validation (keyboard control)
Eyelink('WaitForModeReady', 500);
Screen('FillRect', win.hndl, win.bkgcolor); 

Screen('Flip', win.hndl);   
%%
tic
Eyelink('message','TRIALID %d', 1);
    pause(.1)
 toc
% pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
% PsychPortAudio('FillBuffer', pahandle, wavedata);                           % Fill the audio playback buffer with the audio data 'wavedata':
toc
    pause(.1)
t1          = PsychPortAudio('Start', pahandle, 0, 0, 0);
 EyelinkDoDriftCorrect2(win.el,win.res(1)/2,win.res(2)/2,0)          % drift correction 

    Eyelink('StartRecording');  
    Screen('FillRect', win.hndl, 200); 

    pause(3)
Screen('Flip', win.hndl);   


toc
     Eyelink('StopRecording')

      PsychPortAudio('Stop', pahandle);                                       % Stop the white noise after the last trial
 toc
 
%       PsychPortAudio('Close', pahandle);  
