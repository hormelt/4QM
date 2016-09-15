function [A,B,C,D] = FQM(subData)

% FQM returns the 3-dimensional subpixel locations of particles in subdata
% using the Four Quadrant Method.
%
% INPUTS:
%   subData: image containing the particle to be tracked.
%   fake_dx:
%   fake_dy:
%   calibrate: TRUE to have this function perform calibration
%   p_coef:
%
% OUTPUTS:
%   res: linear corrections in position in the form p1(1)*x+p1(2)...
%       p2(1)*y+p2(2) are output as res = [p1 errx p2 erry]
%

cutoff = round(size(subData,1)/2); %are we requiring even, odd input?

% Define quadrants
QLR = subData(cutoff+1:end,cutoff+1:end,:); QUR = subData(1:cutoff,cutoff+1:end,:);
QLL = subData(cutoff+1:end,1:cutoff,:); QUL = subData(1:cutoff,1:cutoff,:);

% Perform sums to obtain bias
A = squeeze(sum(sum(QUL))); B = squeeze(sum(sum(QUR)));
C = squeeze(sum(sum(QLL))); D = squeeze(sum(sum(QLR)));
end