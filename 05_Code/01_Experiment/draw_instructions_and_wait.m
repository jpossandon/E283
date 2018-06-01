function draw_instructions_and_wait(text,bkg,win,in_dev,ws)

Screen('FillRect', win, bkg);
Screen('Flip', win);
   
DrawFormattedText(win, text, 'center', 'center',255, 60);
Screen('Flip', win);
if in_dev == 1                                                              % Waiting for input according to decided device to continue
    waitForKB_linux({'space'});                                             % press the space key in the keyboard
elseif in_dev == 2
    GetClicks(win,0);                                                       % mouse clik
end
Screen('FillRect', win, bkg);
Screen('Flip', win);
WaitSecs(ws);