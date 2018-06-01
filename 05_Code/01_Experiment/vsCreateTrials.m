function vsTrials = vsCreateTrials(initDur, Nreps,Ncols,Nrows)
% vsTrials = vsCreateTrials() creates randomized trials for visual search.
% They are used by vsRunTrials.m to draw stimuli and control flow.
% The fields of vsTrials are declared general enough to only be filled with
% content (e.g. subject responses); no new fields should be defined later.
%
% Created 11/2012 by Johannes Keyser (jkeyser@uos.de)
%     Rev 01/2013 jkeyser: -response (RT is enough) -screen struct +initDur
%     Rev 11/2013 jkeyser: +fixed bug indicating wrong target repetition

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

if nargin < 1, initDur = 6; end % default is 6 seconds
if nargin < 2, Nreps = 3; end
N.grdX = Ncols; % horizontal amount of (evenly spaced) positions on grid
N.grdY = Nrows; %   vertical amount of (evenly spaced) positions on grid
N.reps = Nreps; % number of repetitions of any condition
N.trls = N.reps*N.grdX*N.grdY; % how many trials will be generated in total
%%% initialize a strctuct with all the fields a vsTrial should ever need
init = struct(...
 'task',    'vsearch',... % task ID the trial is part of; read-only
 'number',  nan,      ... % n-th in presentation order; 1st = 1, read-only
 'repttion',nan,      ... % n-th repetition of trial;   1st = 0, read-only
 'duration',nan,      ... % duration [sec]; read, then set by vsRunTrials.m
 'clockS',  nan,      ... % MATLAB clock now() of START, incl. any setup
 'clockE',  nan,      ... % MATLAB clock now() of END,   incl. any setup
 'stimulus',struct(), ... % vsearch grid parameters; read by vsRunTrials.m
 'RT',      nan,      ... % subject's reaction time; set by vsRunTrials.m
 'valid',   nan,      ... % true / false; set later by analysis scripts
 'fixmat',  struct());    % run's fixmat; set later by analysis scripts
vsTrials = repmat(init, N.trls, 1);
%%% make actual settings
% desired DURATIONS of all trials; will be overwritten with actual duration
[vsTrials.duration] = deal(initDur); % in seconds
% create STIMULUS description for each trial, and set REPETITION number
for pp = 1:(N.grdX*N.grdY) % each position in grid...
    for rep = 1:N.reps     % ...repeated N.reps times
        t = (pp-1)*N.reps +rep; % linear trial index
        vsTrials(t).stimulus.posGrid = stimDescription([N.grdY N.grdX],pp);
        vsTrials(t).stimulus.tgtIndx = pp; % target index, for convenience
    end
end
%%% randomize the trials (TODO: how? chunk-wise?)
% TODO: set rand stream? -> currently done in experiment.m...
vsTrials = vsTrials(randperm(length(vsTrials)));
% set NUMBERS for all trials: ascending integers, beginning from 1
numbercell = num2cell(int16(1:N.trls));
[vsTrials.number] = deal(numbercell{:});
% set REPETITIONS for all trials: 0=1st presentation; 1=1st repetition
stimuli = [vsTrials.stimulus];
tgtInds = [stimuli.tgtIndx];
[bla,idx] = sort(tgtInds); % index aligning same target indices
tgtreps = num2cell(repmat(0:N.reps-1, 1, N.grdX*N.grdY));
[vsTrials(idx).repttion] = deal(tgtreps{:});

    function posGrid = stimDescription(gridDims, tgtIndx)
        % A search array is sufficiently described by a 2D-logical array
        % indicating 'true' in the place where the TARGET is, and 'false'
        % otherwise (where DISTRACTORS are placed).
        % While the visual appearance of target and distractors is defined
        % elsewhere (in vsRunTrials.m), the mapping from indices to
        % positions shall be "1:1". This means that the layout of 
        % disp(posGrid) shall be "identical" to what is seen as actual
        % stimulus, so posGrid denotes locations in an "image coordinate
        % system", where the y-direction is inverted, and the origin in the
        % upper left. Of course, grid spacing is done in vsRunTrials.m.
        %
        % INPUT
        %   gridDims - [i j] grid size in vertical and horizontal dimension
        %   tgtIndx  - target index in grid of singular(!) target position
        %
        % OUTPUT
        %   posGrid - logical array of size gridDims, true only at tgtIndx
        
        % try to convert between all possible indexing schemes
        if isscalar(tgtIndx) % single to subscript indexing
            [tgt_i tgt_j] = ind2sub(gridDims, tgtIndx);
        elseif numel(tgtIndx)==2 && isnumeric(tgtIndx) % subscript indexing
            tgt_i = tgtIndx(1);
            tgt_j = tgtIndx(2);
        elseif islogical(tgtIndx) % logical to subscript indexing
            [tgt_i tgt_j] = ind2sub(gridDims, find(tgtIndx));
        else
            error('Unknown way of indexing a 2D-array!')
        end
        posGrid = false(gridDims);
        posGrid(tgt_i, tgt_j) = true;
    end
end