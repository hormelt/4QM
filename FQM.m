function [A,B,C,D] = FQM(subData)

% FQM returns the 3-dimensional subpixel locations of particles in subdata
% using the Four Quadrant Method.
%
% INPUTS:
%   subData: image containing the particle to be tracked.
%
% OUTPUTS:
%   A,B,C,D: sums over the four quadrants [A B; C D]

cutoff = round(size(subData,1)/2); %are we requiring even, odd input?

% Define quadrants
QLR = subData(cutoff+1:end,cutoff+1:end,:); QUR = subData(1:cutoff,cutoff+1:end,:);
QLL = subData(cutoff+1:end,1:cutoff,:); QUL = subData(1:cutoff,1:cutoff,:);

% Perform sums to obtain bias
A = squeeze(sum(sum(QUL))); B = squeeze(sum(sum(QUR)));
C = squeeze(sum(sum(QLL))); D = squeeze(sum(sum(QLR)));
end