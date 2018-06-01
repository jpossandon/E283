function [obj] = stim_test_ET(tostim,wT,tact_freq,stim_dur)
% simple test, will stimulate the stimulator in variable tostim in the
% order they have in teh vector tostim, very wT seconds and at the end all together

% Tester stimulators
[IsConnected IsDummy] = EyelinkInit(0); % open links to eye-tracker
assert(IsConnected==1, 'Failed to initialize EyeLink!')

% tactile settings and setup
%tact_freq = 200; % 200 in Hz
% calculate from this the settings for the stimulators that are handed over
% to define_tactile_states function. cycle_time sets frequency in microseconds
cycle_time              = 1000/tact_freq*1000
% intensity is specified with on_time. if stimulator is on for half of the
% cycle time, intensity is maximal
on_time = cycle_time-1/2*cycle_time            % TODO: this war wrong before?
%stim_dur = .025% before .025

% this opens the serial port
device = '/dev/ttyS0';          % this is the serial port
obj = serial(device,'BaudRate', 115200,'Parity','none','StopBits',1,'FlowControl','none'); % and specs, this is D/P/S = 8/N/1 which it seems to be kind of standard
try                 % try to open the serial port if is nottouch_randomization.mat
    if strcmpi(obj.Status, 'closed')
        fopen(obj);
    end
catch
    error('MATLAB:serial:fwrite:openserialfailed', lasterr);
    fclose(obj);
end
fprintf('done');%exp_path

% this tells the box what to do through the serial port and according to
% what is written in define_tact_states
% 1 - stimulator 1
% 2 - stimulator 2
% 3 - stimulator 1 and 2
fprintf('\nwriting to serial port...');
ok = define_tact_states(obj, on_time, cycle_time);      % here I am trying to do bilateral stimulation with only one code, do not know if possible, this needs to be done once to setup the stimulator to a given code for the different stimulation channels and stimulation frequency
if ok ~= 1
    error('\n failed writing stimulus codes to serial port!\n');
end
fprintf('done');

% simple test, will stimulate the stimulator in variable tostim in the
% order described and at the end all together
%tostim = [1,2];
%wT     = 2; %time between stimulations
for e = 1:3
for t = 1:length(tostim)
    WaitSecs(wT);
    sprintf('Stimulating #%d',tostim(t))
    Eyelink('command', sprintf('!*write_ioport 0x378 %d',tostim(t)));                    % start stimulation by sending a signal through the parallel port (a number that was set by js_E174_define_tact_states)
    WaitSecs(stim_dur);       
    Eyelink('command', '!*write_ioport 0x378 10');    % stop stimulation
end
end
WaitSecs(wT);
sprintf('All Stimulators at the same time')
Eyelink('command', '!*write_ioport 0x378 9');                    % start stimulation by sending a signal through the parallel port (a number that was set by js_E174_define_tact_states)
WaitSecs(stim_dur);       
Eyelink('command', '!*write_ioport 0x378 10'); 
Eyelink('Shutdown') 
fclose(obj)
fprintf('\nTest ended\n');