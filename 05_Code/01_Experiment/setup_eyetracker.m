function [el cmds] = setup_eyetracker(win, ver)
% [el cmds] = setup_eyetracker(run)
% Sets up custom settings for a Desktop Remote EyeLink 1000 gaze tracker,
% for the configuration defined in a "win" struct, e.g. by experiment.m.
%
% INPUT
%   win - struct from experiment.m, defining an experiment run with fields:
%      .scr : struct, defining screen settings, output of setup_geometry.m
%      .win : struct, defining settings of the Psychtoolbox main window
%      .DoDummyMode : logical, whether to initialize EyeLink in Dummy Mode
%   ver - [0,1,2] verbosity level; what to do errors: nothing/warning/error
%
% OUTPUT
%   el - struct from EyelinkInitDefaults(), with some settings modified
% cmds - command strings and return values of all sent Eyelink('Command')s
%
% Most comments are copied from the EyeLink Programmer's Guide (ELPG) and
% the respective INI file on the Host PC.
% Any command CMD found in the INI-files (the Manuals are incomplete!) can
% be sent as an EyelinkCmd( 'CMD') over the network link.
% The INI-files are provided by SR Research in their Host PC application
% (downloadable from their support forum¹), the NBP-svn source tree², or
% may be copied from the EyeLink Host PC to a USB stick (in Windows mode).
% ¹https://www.sr-support.com/
% ²TRAC url is https://ikw.uni-osnabrueck.de/trac/nbp/, search "INI files"
%
% Created 11/2012 by Johannes Keyser (jkeyser@uos.de), with various help
%                    from Jose Ossandon, Hannah Knepper, Benedict Ehinger.
% Rev. 11/2012 jkeyser: +comments, -cleaning up, +explicit parameters
% Rev. 12/2012 jkeyser: +comments, +commands, +automatic cal/val, +comments
% Rev. 01/2013 jkeyser: made run-struct the only parameter: geometry is now
%                       set dynamically, overwriting PHYSICAL.INI settings;
%                       +function EyelinkCmd() wrapping Eyelink('Command')
% Rev. 02/2013 jossando: changefor using in experiment touch,

if nargin < 1, error('Need a "run" struct as first argument!'), end
if nargin < 2, ver = 1; end % per default, switch on warnings
assert(~isempty(intersect(ver, [0 1 2])), 'valid verbosity levels: 0,1,2!')
if ~isfield(win,'bkgcolor'),win.bkgcolor=127;end %normally set by caller
% Connect to the Eyetracker in the requested mode
[IsConnected IsDummy] = EyelinkInit(win.DoDummyMode); % this is already in the experiment code?
if IsDummy, warning('SetupEL:dummy','EyeLink is in Dummy Mode!'), end
if ~IsConnected, error('Failed to initialize EyeLink!'), end
% initialize with defaults, tell the Eyetracking on what screen to draw
el = EyelinkInitDefaults(win.hndl);
% accumulate sent strings + respective return values of Eyelink('Command'):
cmds = struct('str','','ret',0);

    function cmd = EyelinkCmd(varargin)
        % Instead of using the Eyelink toolbox function directly, let's use
        % Eyelink('Command') through this shallow wrapper. By doing so:
        % 1) you may use the allmighty, familiar options of sprintf(); thus
        % 2) you can inspect exactly what is sent to the eyetracker Host PC
        % which 2a) makes it easier to debug, and 2b) enables complete
        % logging of commands and their return values!
        % Also it detects basic errors and reacts more MATLAB-like on them.
        cmd.str = sprintf(varargin{:}); % let sprintf() handle all varargin
        cmd.ret = Eyelink('Command', cmd.str); % send str, get return value
        if ver && cmd.ret % only for verbosity level >0 and error code ~= 0
            if ver==1, reactfun = @warning; else reactfun = @error; end
            reactfun(sprintf('%s returned with code %d!',cmd.str, cmd.ret))
        end
    end

%%% ELPG 25.4.1 / PHYSICAL.INI
% screen_pixel_coords = <left> <top> <right> <bottom>
% Sets the gaze-position coordinate system, which is used for all
% calibration target locations and drawing commands. Usually set to
% correspond to the pixel mapping of the subject display. Issue the
% calibration-type command after changing this to recompute fixation target
% positions. You should also write a DISPLAY_COORDS message to the start of
% the EDF file to record the display resolution.
% Arguments:
% <left>   X coordinate of  left  of display area
% <top>    Y coordinate of   top  of display area
% <right>  X coordinate of right  of display area
% <bottom> Y coordinate of bottom of display area
% Example: screen_pixel_coords = 0.0 0.0 1024.0 768.0
% NOTE: My idea to make the eyetracker aware of only a centrally placed,
%       "virtual" stimulus area ("coordrect = [0 0win.vrt.res];") does
% b     not work easily, for the following reason: The Eyelink toolbox
% y     uses the MEX-callback PsychEyelinkDispatchCallback() by default,
%       which mainly enables drawing of the Camera Image on the Display PC,
% j     but also making it quite hard to tell all the functions about such
% k     a windows's sub-area we care about. My attempt lead to correct
% e     drawings to the HostPC's gaze curser (my goal), but then all the
% y     Eyelink toolbox's functions about drift-check and calibr/validation
% s     use wrong coordinates for drawing the targets!
% e     I can imagine 3 ways around this, neither of which I have time for:
% r     1? Extending Eyelink's MEX code (but this is very complex).
%       2? Defining your own callback dispatcher (probably best idea),
%       3? Not using callback-code (no Camera Image, needs >1 m-file edits)
coordrect = Screen('Rect',win.hndl); % Eyelink is aware of whole window
cmds(end+1) = EyelinkCmd('screen_pixel_coords = %d %d %d %d', coordrect);
% NOTE: this was also already done in EyelinkInitDefaults()

% screen_distance = <mm to center> | <mm to top> <mm to bottom>
% Used for visual angle and velocity calculations.
% Providing <mm to top> <mm to bottom> parameters will give better
% estimates than <mm to center>
% <mm to center> = distance from display center to subject, in MILLIMETERS
% <mm to top>    = distance from display top    to subject, in MILLIMETERS
% <mm to bottom> = distance from display bottom to subject, in MILLIMETERS
% Example from PHYSICAL.INI: screen_distance = 600 660
% TODO: use top + bottom? -> measure and set in setup_geometry.m ..?
cmds(end+1)=EyelinkCmd('screen_distance = %d',win.Vdst*10);% cm->mm

% Give warning in case of questionable viewer distance setting (supposing
% that the camera-eye distance is the same as the screen-eye distance..!).
% From EyeLink 1000 User Manual version (9/13/2007, 1.3.0),
% Operational / Functional Specifications for Desktop Mount: Remote Option
% Optimal Camera-Eye Distance: Between 40-70 [cm]
if win.Vdst < 40 ||win.Vdst > 70
    warning('Viewing distance (%.1fcm) outside of recommended range!', ...
           win.Vdst) %#ok<WNTAG>
end

% screen_phys_coords = <left>, <top>, <right>, <bottom>
% NOTE: (jkeyser) commas don't matter
% Meaure the distance of the visible part of the display screen edge
% relative to the center of the screen (measured in in MILLIMETERS).
% Position of display area corners relative to display center
% Example (PHYSICAL.INI): screen_phys_coords = -188.0, 146.0, 188.0, -146.0
% UP is positive, DOWN is negative; RIGHT is positive, LEFT is negative
% NOTE: scr.wdth and scr.hght come in [cm] while function assumes [mm]!
cmds(end+1) = EyelinkCmd('screen_phys_coords = %.1f, %.1f, %.1f, %.1f', ...
                         -.5*win.wdth*10, ...
                         +.5*win.hght*10, ...
                         +.5*win.wdth*10, ...
                         -.5*win.hght*10);
% sanity check; factor 10 should be detectable
if   win.wdth < 10 ||win.wdth > 100 ...
   ||win.hght < 10 ||win.hght > 100
    warning('Sure your virtual area is %gcm x %gcm (CENTIMETERS)?', ...
           win.wdth,win.hght) %#ok<WNTAG>
end

%%% Make sure that we get gaze data from the Eyelink and that we have the
%%% all information in the EDF file (cf. User Manual 4.4 - File Data Types
%%% and ELPG pages 179-182, File/Link Data Control).
% Sample data:
%   TODO: ELPG 75 lower example defines also LEFT,RIGHT?!?
%   GAZE    - screen xy (gaze) position
%   GAZERES - units-per-degree screen resolution
%   HREF    - head-referenced gaze
%   PUPIL   - raw eye camera pupil coordinates
%   AREA    - pupil area
%   STATUS  - warning and error flags
%   BUTTON  - button state and change flags
%   INPUT   - input port data lines
% Event data:
%   GAZE     - screen xy (gaze) position
%   GAZERES  - units-per-degree angular resolution
%   HREF     - HREF gaze position
%   AREA     - pupil area or diameter
%   VELOCITY - velocity of eye motion (avg, peak)
%   STATUS   - warning and error flags for event
%   FIXAVG   - include ONLY average data in ENDFIX events
%   NOSTART  - start events have no data, just time stamp
% Event filter:
%   LEFT,RIGHT - events for one or both eyes
%   FIXATION   - fixation start and end events
%   FIXUPDATE  - fixation (pursuit) state updates
%   SACCADE    - saccade start and end
%   BLINK      - blink start an end
%   MESSAGE    - messages (user notes in file)
%   BUTTON     - button 1..8 press or release
%   INPUT      - changes in input port lines
%
%%% which data should be saved to EDF file? %TODOs here!!
cmds(end+1) = EyelinkCmd( ...
    'file_sample_data = GAZE,GAZERES,HREF,AREA,STATUS,BUTTON');
cmds(end+1) = EyelinkCmd( ...
	'file_event_data = GAZE,GAZERES,HREF,AREA,VELOCITY,STATUS');
cmds(end+1) = EyelinkCmd( ...
    'file_event_filter=LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
%%% Which data should be sent through the link?
cmds(end+1) = EyelinkCmd('link_sample_data = GAZE,STATUS,BUTTON');
cmds(end+1) = EyelinkCmd('link_event_data = GAZE,STATUS');
cmds(end+1) = EyelinkCmd(...
    'link_event_filter = LEFT,RIGHT,SACCADE,FIXATION,BLINK,MESSAGE,BUTTON');
%%% From ANALOG.INI: Select type of data for analog output
% 	OFF    turns off analog output
% 	PUPIL is raw pupil x,y
% 	HREF  is headref-calibrated x,y
%   GAZE  is screen gaze x,y
cmds(end+1) = EyelinkCmd('analog_out_data_type = OFF');
%%% von hannah: (TODO: consider)
% remote mode possible add HTARGET (head target)
% cmds(end+1) = EyelinkCmd( ...
%  'file_sample_data  = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS');
% cmds(end+1) = EyelinkCmd( ...
%  'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
% % set link data (used for gaze cursor)
% cmds(end+1) = EyelinkCmd( ...
%     'link_event_filter = LEFT,RIGHT,FIXATION,BUTTON');
% cmds(end+1) = EyelinkCmd( ...
%     'link_sample_data  = LEFT,RIGHT,GAZE,AREA');

%%% ELPG 25.6.4 heuristic_filter (EyeLink II, EyeLink 1000)
% heuristic_filter = <ON or OFF>
% heuristic_filter = <linkfilter> <filefilter>
% Can be used to set level of filtering on the link and analog output, and
% on file data. An additional delay of 1 sample is added to link or analog
% data for each filter level. If an argument of <on> is used, link filter
% level is set to 1 to match EyeLink I delays. The file filter level is not
% changed unless two arguments are supplied. The default file filter level
% is 2. Arguments:
% 0 or OFF - disables link filter
% 1 or ON  - sets filter to 1 (moderate filtering, 1 sample delay)
% 2        - applies an extra level of filtering (2 sample delay)
% Example from DEFAULTS.INI: heuristic_filter 1 2
cmds(end+1) = EyelinkCmd('heuristic_filter 1 2');
% NOTE: gets set again with 'heuristic_filter ON' in EyelinkDoTrackerSetup

%%% Parser Configuration
% ELPG 25.10.8  Typical Parser Configurations
% The Cognitive configuration is optimal for visual search and reading, and
% and ignores most saccades smaller than 0.5°. It's less sensitive to setup
% problems. The Pursuit configuration is designed to detect small (<0.25°)
% saccades, but may produce false saccades if subject setup is poor.
% ELPG 25.10.1  select_parser_configuration = <set>
% EyeLink II and EyeLink1000 ONLY! Selects the preset standard parser setup
% (0) or more sensitive saccade detector (1). These are equivalent to the
% cognitive and psychophysical configurations listed below. Arguments:
% <set> - 0: standard, 1:high sensitivity saccade detector configuration.
cmds(end+1) = EyelinkCmd('select_parser_configuration = 0');
% The remaining commands are documented for use with EyeLink I only.
% Cognitive Configuration:
% recording_parse_type           = GAZE
% saccade_velocity_threshold     = 30
% saccade_acceleration_threshold = 9500
% saccade_motion_threshold       = 0.15
% saccade_pursuit_fixup          = 60
% fixation_update_intecmdsal       = 0
% Pursuit and Neurological Configuration:
% recording_parse_type           = GAZE
% saccade_velocity_threshold     = 22
% saccade_acceleration_threshold = 5000
% saccade_motion_threshold       = 0.0
% saccade_pursuit_fixup          = 60
% fixation_update_intecmdsal       = 0

% set pupil Tracking model [in Remote mode, only "ellipse" exists!]
% no = centroid, yes = ellipse
cmds(end+1) = EyelinkCmd( 'use_ellipse_fitter = YES');
% set sample rate in camera setup screen
cmds(end+1) = EyelinkCmd( 'sample_rate = 500');

%%% ELPG 25.2.1 / CALIBR.INI  calibration_type = <type>
% Sets calibration type and recomputed the calibration targets after a
% display resolution change. Arguments:
% <type>: one of these calibration type codes:
% H3   - horizontal 3-point calibration, horizontal-only 3-point quadratic
% HV3  -  3-point calibration, poor linearization, bilinear
% HV5  -  5-point, poor at corners, biquadratic
% HV9  -  9-point, best overall, biquadratic with corner correction
% HV13 - 13-point, for large calibration region (ELII >=ver2.0 or EL 1000)
%%% from CALIBR.INI:
% HV13 = 13-point, bicubic calibration
% HV13 works best with larger angular displays (> +/-20°). It should NOT be
% used when accurate data is needed from corners of calibrated area!
cmds(end+1) = EyelinkCmd('calibration_type = HV13');
% calibration sequencing "YES": auto, "NO": manual. With PD patients: NO!
cmds(end+1) = EyelinkCmd('enable_automatic_calibration = NO');
% should the Eyetracker generate its default target placements
cmds(end+1) = EyelinkCmd('generate_default_targets = NO');
% should the calibration/validation points appear in random order?
cmds(end+1) = EyelinkCmd('randomize_calibration_order = YES');
cmds(end+1) = EyelinkCmd('randomize_validation_order = YES');
% number of points to resample after validation [default 2]
cmds(end+1) = EyelinkCmd('validation_resample_worst = 2');
% error required for resampling in degrees [default 1]
cmds(end+1) = EyelinkCmd('validation_worst_error = 1.5');
% which error is unacceptable during drift correction? [default 2]
cmds(end+1) = EyelinkCmd('drift_correction_rpt_error = 1.5');
% do we show drift correction during validation?
cmds(end+1) = EyelinkCmd('validation_online_fixup = YES');
% Should the calibration messages be printed in the EDF file?
cmds(end+1) = EyelinkCmd('auto_calibration_messages = YES');
% Should the 1st point of the calibration and validation be repeated?
cmds(end+1) = EyelinkCmd('cal_repeat_first_target = YES');
cmds(end+1) = EyelinkCmd('val_repeat_first_target = YES');
%%% forum post:
% https://www.sr-support.com/showthread.php?2468-Custom-calibration-points/
% page2&highlight=hv13+positions
% If you simply want to shrink/compress the location of the targets towards
% the center of the screen (with the center-most point remaining where it
% is by default) then you can adjust the following parameters:
% calibration_area_proportion = <x,y display proportion>
% validation_area_proportion  = <x,y display proportion> (from CALIBR.INI)
% For auto generated calibration point positions, these set the part of the
% width or height of the display to be bounded by. Targets each may have a
% single proportion, or a horizontal followed by a vertical proportion.
% Default for both calibration and validation are: 0.88 and 0.83
% NOTE: setting for calibration also sets validation
%%% set according to preset defined in setup_geometry
% prprts =win.vrt.res ./win.res; % maximally yields [1 1]
% % don't enlarge the proportions w.r.t. Default (too close to screen border)
% prprts = min([0.88 0.83; prprts]);
% cmds(end+1) = EyelinkCmd('calibration_area_proportion = %.2f %.2f',prprts);
% cmds(end+1) = EyelinkCmd( 'validation_area_proportion = %.2f %.2f',prprts);

%%% Setup Cedrus Response Box buttons (test this with EyeLinkPPortCapture)
% EL Programmer's Guide (ELPG), p171: 25.3.5  write_ioport <ioport> <data>
% Writes data to I/O port. Useful to configure I/O cards.
% For this case of reading the parallel port, we need bi-directional!
%
% From FINAL.INI example (from SR support forum, link above):
% write_ioport 0x37A 0x20  ;enables bidirectional mode [what we need here]
% input_data_ports = 0x379 ;use the status register to read input data 
% input_data_masks = 0xFF  ;mask port bits, so changes don't trigger events
%
% EyelinkCmd( 'write_ioport 0x37A 0x0')   % UNI-directional
% EyelinkCmd( 'input_data_ports = 0x379') % UNI-directional
% cmds(end+1) = EyelinkCmd('write_ioport 0x37A 0x20'); % BI-directional
% cmds(end+1) = EyelinkCmd('input_data_ports = 0x378');% BI-directional
% cmds(end+1) = EyelinkCmd('input_data_masks = 0xFF'); % TODO: is it needed?
% 
% ELPG p178 25.3.2  create_button <button> <ioport> <bitmask> <inverted>
% Defines a button to a bit in a hardware port. Arguments:
% <button>   - button number, 1 to 8
% <ioport>   - address of hardware port
% <bitmask>  - 8-bit mask ANDed with port to test button line
% <inverted> - 1 if active-low, 0 if active-high
% Example: EyelinkCmd( 'create_button 4 0x378 0x08 0')
% CREATE_BUTTON = 'create_button %d 0x378 0x%s 0';
% for bt = 1:8
%     bitmask = dec2hex(2^(bt-1)); % [0x] 01, 02, 04, 08, 10, 20, 40, 80
%     cmds(end+1) = EyelinkCmd(CREATE_BUTTON, bt, bitmask); %#ok<AGROW>
% end

% ELPG 25.3.3  button_function <button> <presscmd> <relcmd>
% Assigns a command to a button. This can be used to control recording with
% a digital input, or to let a button be used instead of the spacebar
% during calibration. Arguments:
% <button>   - hardware button 1 to 8, keybutton 8 to 31
% <presscmd> - command to execute when button pressed
% <relcmd>   - command to execute when button released
%
% From FINAL.INI example (from SR support forum, link above):
% Buttons may be assigned to generate commands
% Give the button number (1..31), then the command string in quotes.
% The first command is executed when the button is pressed.
% The optional second is executed when button is released.
% Giving NO strings deletes the button function.
% BUTTON_FUNCTION = ['button_function %d ' ...
%                    '"data_message ''BUTTON_%d_PRESSED''" '...
%                    '"data_message ''BUTTON_%d_RELEASED''"'];
% for bt = 1:8
%     cmds(end+1) = EyelinkCmd(BUTTON_FUNCTION, bt,bt,bt); %#ok<AGROW>
% end

% ELPG p. 63, 12.1 Calibration Colors
% "The entire display is cleared to target_background_color before
% calibration and drift correction, and this is also the background for the
% camera images. The background color should match the average brightness
% of your experimental display, as this will prevent rapid changes in the
% subject’s pupil size at the start of the trial. This will provide the
% best eye-tracking accuracy as well.
% [...]
% A background color of black with white targets is used by many 
% experimenters, especially for saccadic tasks. A full white (255,255,255)
% or a medium white (200,200,200) background with black targets is 
% preferable for text. Using white or gray backgrounds rather than black
% helps reduce pupil size and increase eye-tracking range, and may reduce
% retinal afterimages."
% Also https://www.sr-support.com/showthread.php?19-EyeLink-MATLAB-Toolbox
el.foregroundcolour = 255; % CLUT color idx (for Psychtoolbox functions)
el.backgroundcolour =win.bkgcolor;
el.msgfontcolour    = 200;
% ELPG 12.2  Calibration Target Appearance
% The standard calibration and drift correction target is a filled circle
% (for peripheral delectability) with a central "hole" target (for accurate
% fixation). The sizes of these features may be set with
% set_target_size(diameter, holesize). If holesize is 0, no central feature
% will be drawn. The disk is drawn in the calibration foreground color, and
% the hole is drawn in the calibration background color.
%%% Size and color of the overall calibration target
% If you are later using EyelinkDoTrackerSetup(), which calls
% EyelinkTargetModeDisplay(), which in turn calls EyelinkDrawCalTarget(),
% the units of calibrationtargetsize and calibrationtargetwidth are
% percentage of the PTB-window's width; just set to taste?
% Example in ELPG suggests size: SCRWIDTH/60, hole: SCRWIDTH/300 [pixels]
el.calibrationtargetsize   = 1; % [% of window width]
% size of white dot in the center of the calibration target 
el.calibrationtargetwidth  = .3; % [% of window width]
el.calibrationtargetcolour = 255; % CLUT color index
% beep when a target is presented? default true
el.targetbeep = false; 
el.feedbackbeep = false;
% show calibration results on Display Computer? default true
el.displayCalResults = true;
% el.eyeimgsize = 50; % [percentage of screen], default 30
el.allowlocaltrigger = true; % allow user to trigger him or herself
el.allowlocalcontrol = true; % allow control from subject-computer
EyelinkUpdateDefaults(el) % pass changes to the MEX-callback function

%%% If you set 'generate_default_targets = NO', you have to specify the
%%% coordinates of calibration and validation targets manually:
% from CALIBR.INI:
% Specified target locations
% Orderings of points for 9-point and 13-point calibration / validations:
%     HV9                    HV13
%  6   2   7              6   2   7
%                           10  11
%  4   1   5              4   1   5
%                           12  13
%  8   3   9              8   3   9
%
% some manual pixel coordinates for HV13 on a 1920x1080 monitor (x,y pairs)
% scr_w = 1920;
% scr_h = 1080;
 scr_w = win.res(1);
 scr_h = win.res(2);
% % left, left-middle, middle, middle-right, right
% ll = 510;
% lm = 735;
% mm = scr_w/2;
% mr = 1185;
% rr = 1410;
% % top, top-center, center, center-bottom, bottom
% tt = 140;
% tc = 340;
% cc = scr_h/2;
% cb = 740;
% bb = 940;
%%% Bene's (behinger@uos.de) way to calculate HV13 positions:
 scale = 6;
% % left, left-middle, middle, middle right, right
ll = scr_w/2 - 2*scr_w/scale;
lm = scr_w/2 - 1*scr_w/scale;
mm = scr_w/2 + 0*scr_w/scale;
mr = scr_w/2 + 1*scr_w/scale;
rr = scr_w/2 + 2*scr_w/scale;
% % top, top-center, center, center-bottom, bottom
tt = scr_h/2 - 2*scr_h/scale;
tc = scr_h/2 - 1*scr_h/scale;
cc = scr_h/2 + 0*scr_h/scale;
cb = scr_h/2 + 1*scr_h/scale;
bb = scr_h/2 + 2*scr_h/scale;
%
%%% test plot of calibration/validation targets (just comment-in to run)
% %% [code-cell start]
% coords = [mm,cc; mm,tt; mm,bb; ll,cc; rr,cc; ...
%           ll,tt; rr,tt; ll,bb; rr,bb; lm,tc; lm,cb; mr,tc; mr,cb];
% figure(), ax = axes();
% scatter(ax, coords(:,1), coords(:,2), 50, 'filled',...
%         'MarkerFacecolor', repmat(200, 1,3)/255)
% set(ax, 'Color', repmat(el.backgroundcolour, 1,3)/255)
% set(ax, 'YDir', 'reverse')
% axis image, axis([0 scr_w 0 scr_h])
% %% [code-cell end]
%
%%% calibration coordinates
cmds(end+1) = EyelinkCmd( ...
                ['calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d'...
                 ' %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d'], ...
                 mm,cc, mm,tt, mm,bb, ll,cc, rr,cc, ...
                 ll,tt, rr,tt, ll,bb, rr,bb, lm,tc, lm,cb, mr,tc, mr,cb);
%% validation coordinates 
cmds(end+1) = EyelinkCmd( ...
                 ['validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d'...
                 ' %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d'], ...
                 mm,cc, mm,tt, mm,bb, ll,cc, rr,cc, ...
                 ll,tt, rr,tt, ll,bb, rr,bb, lm,tc, lm,cb, mr,tc, mr,cb);
end