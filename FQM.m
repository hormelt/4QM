function [A,B,C,D] = FQM(subData)

% FQM returns the sum of the 4 quadrants from a stack of frames.
%
% INPUTS:
%   subData: Stack of images.
%
% OUTPUTS:
%   [A,B,C,D]: Sum over the four quadrants UL,UR,LL and LR.

% Find the middle of the frame.
Middle = size(subData,1)/2;

if (mod(Middle,1) >= 0.1) && (mod(Middle,1) <= 0.9)
    disp([char(9) 'Warning: The MiddleCoordinate is placed in a pixel.'])
end

% Define quadrants.
QLR = subData(Middle+1:end,Middle+1:end,:); 
QUR = subData(1:Middle,Middle+1:end,:);
QLL = subData(Middle+1:end,1:Middle,:); 
QUL = subData(1:Middle,1:Middle,:);

% Perform sums to obtain bias.
A = squeeze(sum(sum(QUL))); 
B = squeeze(sum(sum(QUR)));
C = squeeze(sum(sum(QLL))); 
D = squeeze(sum(sum(QLR)));
end
