function posVec = vsCreateStimulus(vsPosGrid, scr, winH, draw2EL,borderH,borderV)
% Create stimulus array for visual search task; E.g. used by vsRunTrials.m.
%
% INPUT
% vsPosGrid - target/distractor position grid of trials by vsCreateTrials.m
%       scr - setup_geometry.m's output (+field 'bkgcolor', background color)
%      winH - Psychtoolbox window handle to draw to
%   draw2EL - true/false; whether to also draw on EyeLink Host PC
%   borderH - screen horizontal margin not used in % of the screen (10 - 10% to the left and 10% to the right)
%
% Created 11/2012 by Johannes Keyser (jkeyser@uos.de)
%     Rev 01/2013 jkeyser: +stimuli defined using "scr", +drawing on HostPC
%
% TODO: Why is antialiasing not working properly, but e.g. it does in
%       EyelinkDrawCalTarget() with the same method? -> directly draw to
%       win.hndl instead of intermediate texture?

%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

AssertOpenGL() % error immediately if problems with Screen()
if nargin < 1, vsPosGrid = []; end
if nargin < 2 || isempty(scr), scr = setup_geometry(); end % good 4 testing
if nargin < 3, winH = []; end
if nargin < 4, draw2EL = false; end
if ~isfield(scr,'bkgcolor'), scr.bkgcolor = 80; end % normally added by caller
if isempty(vsPosGrid)
    % get trial's stimulus descriptions
    vsTrials  = vsCreateTrials();
    vsPosGrid = vsTrials(1).stimulus.posGrid;
    % OR manual control of target position
%     vsPosGrid = false(size(vsTrials(1).stimulus.posGrid));
%     vsPosGrid(3,5) = true;
end
if isempty(winH) % no window handle given by caller: create a new one
    scID = 0; % where to put the window? [0 <=> main screen]
    rect = [];%[0 0 scr.vrt.res]; % empty for fullscreen
    [win.hndl, win.rect] = Screen('OpenWindow', scID, scr.bkgcolor, rect);
    [win.res(1), win.res(2)] = Screen('WindowSize', win.hndl);
    Screen('BlendFunction', win.hndl, GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
else % window handle specified by caller: draw stimulus to that window
    win.hndl = winH;
    [win.res(1), win.res(2)] = Screen('WindowSize', winH);
end
if any(win.res < scr.res), error('Window too small for geometry!'), end
if any(win.res ~= scr.res)
	warning('Geometry differs from expectation!'); %#ok<WNTAG>
end

% get size of stim grid
[grdY, grdX] = size(vsPosGrid);
% restrict stimulus extent to requested area (in pixels [x y] "resolution")
stimRes = min([scr.res; win.res]);
borderH = borderH/100*stimRes(1);
borderV = borderV/100*stimRes(2);
% evenly distribute positions on grid, with same distances to area borders
dsX = (stimRes(1)-2*borderH)/(grdX+1); % one more distance than dots | - o - o - o - |
dsY = (stimRes(2)-2*borderV)/(grdY+1);
% create symbol positions w.r.t. virtual stimulus area
[vrtPosX, vrtPosY] = meshgrid(linspace(dsX+borderH, stimRes(1)-dsX-borderH, grdX),...
                             linspace(dsY+borderV, stimRes(2)-dsY-borderV, grdY));
posVec.vrt = [vrtPosX(:) vrtPosY(:)]'; % 2-row vector, for logging
% add offsets necessary for centered presentation on actual screen
scrPosX = floor( (win.res(1)-stimRes(1))/2 + vrtPosX ); % round to pixels
scrPosY = floor( (win.res(2)-stimRes(2))/2 + vrtPosY ); % round to pixels
posVec.scr = [scrPosX(:) scrPosY(:)]'; % 2-row vector, for PTB on DisplayPC
%%% position grid controls positions of target and distractors
tgt_scrPos = posVec.scr(:,  vsPosGrid); % target
dst_scrPos = posVec.scr(:, ~vsPosGrid); % distractors
tgt_tex = drawElement(win.hndl, 'target');
dst_tex = drawElement(win.hndl, 'distractor');
% draw the textures to (multiple) locations to window handle
PsychDrawSprites2D(win.hndl, dst_tex, dst_scrPos);
PsychDrawSprites2D(win.hndl, tgt_tex, tgt_scrPos);
Screen('Close', [tgt_tex dst_tex]) % free memory from drawn texture patches
if isempty(winH)
	% since caller wants to test stimulus, draw box around stimulus area
    pW = 1; % pen width of box border
    if all(win.res > pW+scr.vrt.res) % only if it fits
        lcolor = 0; % line color
        stimRect = [win.res-stimRes-pW win.res+stimRes+pW]/2;
        Screen('FrameRect', win.hndl, lcolor, stimRect, pW);
    end
    Screen('Flip', win.hndl); % only if drawing to new win
    ListenChar(1); % keyboard-listen on Command Window for typing "sca"
end
% on request, draw a simplified version of the stimulus to the Eyelink Host
if draw2EL, vsShowStimOnHost(), end % NOTE: May change the Eyelink's mode!
%%% calls to test- and documentation functions
% disp(vsTrials(1).stimulus.posGrid)
%
% for saving example images:
%  fname = sprintf('vsStim_Target%02d', find(vsPosGrid));
%  saveScreenshotAsPNG(win.hndl, fname)
%
% testGridPositions() % draw all grid positions for inspection

    function tex = drawElement(winH, elemStr)
        % Draws any vsearch-element to a texture, returns it for caller.
        % Targets are circles ("O"s), while distractors have an additional
        % line (TODO: cite Behrmann et al correctly), i.e. are "Q"s.
        %
        % INPUT
        %   winH    - Psychtoolbox window handle
        %   elemStr - 'target'/'distractor', element to draw
        % OUTPUT
        %   elem_tex - texture of element, later drawn to winH elsewhere
        %
        validElems = {'target','distractor'};
        assert(ismember(elemStr, validElems), ...
               ['Valid : ''' validElems{1} ''' or ''' validElems{2}])
        %%% SETTINGS (coming from "def" parameter)
        DIAM    = scr.pixxdeg.*.75;  % "O/Q" diameter (target/distractors)
        LLEN    = scr.pixxdeg.*.3; % length of line in "Q" (distractors)
        LWIDTH  = 4;   % width of any drawn line %TODO: not for extremes!
        COLOR   = 255; % CLUT index or [r g b a]
        BGCOLOR = [0 0 0 0]; % patch background color, CLUT or [r g b a]
        RECT    = [0 0 DIAM round(DIAM+LLEN/2)];
        CENTER  = round([DIAM/2 DIAM/2]);
        %%% DRAWING
        tex = Screen('OpenOffScreenWindow', winH, BGCOLOR, RECT);
        Screen('BlendFunction', tex, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        Screen('FrameOval', tex, COLOR, [0 0 DIAM DIAM], LWIDTH); % circle
% 		  Screen('FillOval', tex, COLOR, [0 0 DIAM DIAM]);
% 		  insetRect = [0 0 DIAM DIAM]+[LWIDTH LWIDTH -LWIDTH -LWIDTH];
% 		  Screen('FillOval', tex, scr.bkgcolor, insetRect);
        if strcmp(elemStr, 'distractor')
            % draw additional line for "Q" (however with line H-centered)
            fromH = CENTER(1);
            fromV = round(DIAM-LLEN/2);
            toH   = CENTER(1);
            toV   = RECT(4);
            Screen('DrawLine', tex, COLOR, fromH, fromV, toH, toV, LWIDTH);
        end
    end
 

    function vsShowStimOnHost()
        % Draw background graphics on the EyeLink EyeLink's real-time gaze
        % cursor display on the Host PC, to allow for evaluation of subject
        % performance and tracking errors.
        % Simple boxes etc around important details are sufficient.
        %
        % ELPG 25.7 Drawing Commands
        % draw_box        <x1> <y1> <x2> <y2> <color>
        % draw_filled_box <x1> <y1> <x2> <y2> <color>
        % Draws empty/filled box, in gaze-position display coordinates.
        % <x1>,<y1> corner of box
        % <x2>,<y2> opposite corner of box
        % <color>   0 to 15
        %
        % draw_cross <x> <y> <color>
        % <x1>,<y1> center point of cross
        % <color>   0 to 15
        %
        BoxS = scr.pixxdeg.*1; % should match drawElement()'s DIAM setting
        if ~Eyelink('IsConnected')
            warning('No EyeLink: Not drawing to HostPC!') %#ok<WNTAG>
            return
        end
        % 'CheckRecording' indicates a "recording in progress" by return 0,
        % ABORT_EXPT==disconnected, TRIAL_ERROR==other non-recording state.
        % However, ALL commands return 0 during Dummy Mode!
        % => Need to add check for connection status to never destroy data!
        if Eyelink('Isconnected')>0 && Eyelink('CheckRecording')==0
            warning('Record in progress: Won''t draw to HostPC')%#ok<WNTAG>
            return
        end
        Eyelink('SetOfflineMode') % need "Offline Mode" to draw to screen!
        Eyelink('Command', 'clear_screen 1'); % 0==black, 1==dark blue, ...
        % draw target as filled box
        tX = tgt_scrPos(1);
        tY = tgt_scrPos(2);
        Eyelink('Command', sprintf('draw_filled_box %d %d %d %d 7', ...
                        round([tX-BoxS/2 tY-BoxS/2 tX+BoxS/2 tY+BoxS/2])));
        % draw all the targets as empty boxes
        for dp = 1:length(dst_scrPos)
            dX = dst_scrPos(1, dp);
            dY = dst_scrPos(2, dp);
            Eyelink('Command', sprintf('draw_box %d %d %d %d 7', ...
                        round([dX-BoxS/2 dY-BoxS/2 dX+BoxS/2 dY+BoxS/2])));
        end
    end

    function saveScreenshotAsPNG(winH, fname) %#ok<DEFNU>
        if nargin < 2, fname = ['saveScreen' datestr(now(),30)]; end
        [ign, fname] = fileparts(fname); %#ok<ASGLU> % strip file extension
        imgarray = Screen('GetImage', winH);
        savepath = [fullfile(pwd(),fname) '.png'];
        fprintf('Saving window content to "%s" ... ', savepath)
        imwrite(imgarray, savepath, 'png')
        fprintf('DONE.\n')
    end

    function testGridPositions() %#ok<DEFNU>
        % function to inspect the grid positions even without Psychtoolbox
        figure(), ax = axes();
        scatter(vrtPosX(:), vrtPosY(:), 'MarkerEdgeColor', 'k')
        axis(ax,'image')
        axis(ax, [0 win.wdth 0 win.hght])
        set(ax, 'YDir', 'reverse')
    end
end