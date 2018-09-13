function RunBat( filename )
%RUNBAT      runs a batch file via Windows start menu using Java Robot class
%Input:
% <filename> batch file (.bat)
%CL

% Get Windows version
winvers = system_dependent('getos');

% Create java robot to generate native system input events
import java.awt.Robot;
import java.awt.event.*;
robot = Robot();

% Open the windows start menu
scrnsize = get(0,'screensize');
robot.mouseMove(0,scrnsize(4));
robot.mousePress(InputEvent.BUTTON1_MASK);
robot.mouseRelease(InputEvent.BUTTON1_MASK);

% Execute the batch file by typing the filename
if ~isempty(strfind(winvers,'10'))
    
    % For Windows 10, batch execution via 'run'
    pause(1)
    lang = get(0,'Language');
    if strcmp('de',lang(1:2))
        typestr('Ausf�hren')
    else
        typestr('run')
    end
    robot.keyPress(KeyEvent.VK_ENTER)
    pause(1)
    % Type filename
    typestr(filename)
    robot.keyPress(KeyEvent.VK_ENTER)
    
else
    
    % For Windows 7 and 8, batch execution via 'search'
    typestr(filename)
    robot.keyPress(KeyEvent.VK_ENTER)

end

end


function typestr( str )
%TYPESTRING  types a given string on the keyboard using the java robot

% Create java robot to generate native system input events
import java.awt.Robot;
import java.awt.event.*;
robot = Robot();

for i = 1:length(str)
    c = str(i);
    if isstrprop(c,'punct') % Punctuation character
        switch c
            case '-'
                robot.keyPress(KeyEvent.VK_MINUS)
            case '.'
                robot.keyPress(KeyEvent.VK_PERIOD)
            case ':'
                robot.keyPress(KeyEvent.VK_SHIFT)
                robot.keyPress(KeyEvent.VK_PERIOD)
                robot.keyRelease(KeyEvent.VK_SHIFT)
            case '\' % ASCII value
                robot.keyPress(KeyEvent.VK_ALT);
                robot.keyPress(KeyEvent.VK_NUMPAD9);
                robot.keyPress(KeyEvent.VK_NUMPAD2);
                robot.keyRelease(KeyEvent.VK_ALT);               
        end
    elseif isstrprop(c,'wspace') % Space
        robot.keyPress(KeyEvent.VK_SPACE);
    elseif isstrprop(c,'digit') % Number
        eval(strcat('robot.keyPress(KeyEvent.VK_NUMPAD',c,')'));
    elseif strcmp(c,'�')
        robot.keyPress(KeyEvent.VK_ALT);
        robot.keyPress(KeyEvent.VK_NUMPAD1);
        robot.keyPress(KeyEvent.VK_NUMPAD5);
        robot.keyPress(KeyEvent.VK_NUMPAD4);
        robot.keyRelease(KeyEvent.VK_ALT);
    elseif strcmp(c,'�')
         robot.keyPress(KeyEvent.VK_ALT);
        robot.keyPress(KeyEvent.VK_NUMPAD1);
        robot.keyPress(KeyEvent.VK_NUMPAD5);
        robot.keyPress(KeyEvent.VK_NUMPAD3);
        robot.keyRelease(KeyEvent.VK_ALT);
    elseif strcmp(c,'�')
        robot.keyPress(KeyEvent.VK_ALT);
        robot.keyPress(KeyEvent.VK_NUMPAD1);
        robot.keyPress(KeyEvent.VK_NUMPAD4);
        robot.keyPress(KeyEvent.VK_NUMPAD2);
        robot.keyRelease(KeyEvent.VK_ALT);
    else
        c = upper(c); % Upper case character
        eval(strcat('robot.keyPress(KeyEvent.VK_',c,')'));
    end   
    pause(0.1)
end

end