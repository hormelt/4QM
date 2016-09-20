function [correctedMSDs, MSDs, trialCenters] = tracker_caller_4QM(FileStub,varargin)

% Segmentation and Tracking of Particles via the 4QM method in 2D.
%
% INPUTS:
%   FileStub: Path to the image_stack to be analyzed.
%   NmPerPixel: [optional] Actual pixel width (nm).
%   SecsPerFrame: [optional] Time between frames (sec).
%   NoiseSz: [optional] (pixels).
%   FeatSize: [optional] Full optical diameter of particle (pixels).
%   DeltaFit: [optional] Narrows analysis region around particle (pixels).
%   ThreshFact: [optional] maximum intensity devided by the thresfact gives
%   the threshold value.
%   TrackMem: [optional] Number of steps disconnected tracks can be
%   reconnected, in case a particle is lost.
%   Dim: [optional] Dimension of the system.
%   MinTrackLength: [optional] minimum length of track; throw away tracks
%   shorter than this.
%   PrintTrackProgress: [optional] Turns on or off printing progress to screen.
%   MaxDisp: [optional] MaxDisp should be set to a value somewhat less than
%   the mean spacing between the particles.
%   NTests: [optional] Number of trial shifts per particle.
%   FrmStart: [optional] Frame to start data analysis from.
%   PlotOpt: [optional ]Options for plotting data {'simple','bandpassed','none'}
%
% NOTES:
%   The imwrite() function is unstable when windows file explorer is opened.
%
% DEPENDENCIES:
%   This program depends on the particle_tracking toolbox from the TA-lab
%   for the following functions: bpass2D_TA() and calcMSD().
%   This program depends on the SPtrack1.0 toolbox by Eric Dufresne from
%   Yale University for the following functions: pkfnd() and cntrd().
%   This program depends on the trackin() function from Crocker
%   (http://glinda.lrsm.upenn.edu/~weeks/idl).

tic

%% Set up particle, intensity and duration parameters.

% Default Values.

info = imfinfo([FileStub '.tif']);

defNFrames = numel(info);
defNmPerPixel = 16.25; % zyla at 400x
defSecsPerFrame = 0.011179722;
defNoiseSz = 1;
defFeatSize = 15;
defDeltaFit = 3;
defThreshFact = 3.5;
defTrackMem = 0;
defDim = 2;
defMinTrackLength = 100; % This was sufficient statistics for my dissertation research -TH
defPrintTrackProgress = 1;
defMaxDisp = defFeatSize/2;
defp.FrameStart = 1;
defPlotOpt = 'none';
defErrorThresh = 0.1;
defNTests = 100;
defStepAmplitude = 1;
validPlotOpt = {'bandpass','simple','none'};
checkPlotOpt = @(x) any(validatestring(x,validPlotOpt));

% Set up required and optional inputs.
f = inputParser;
f.CaseSensitive = 0;

addRequired(f,'FileStub',@ischar)
addOptional(f,'NFrames',defNFrames,@isnumeric)
addOptional(f,'NmPerPixel',defNmPerPixel,@isnumeric)
addOptional(f,'SecsPerFrame',defSecsPerFrame,@isnumeric)
addOptional(f,'NoiseSz',defNoiseSz,@isnumeric)
addOptional(f,'FeatSize',defFeatSize,@isnumeric)
addOptional(f,'DeltaFit',defDeltaFit,@isnumeric)
addOptional(f,'ThreshFact',defThreshFact,@isnumeric)
addOptional(f,'TrackMem',defTrackMem,@isnumeric)
addOptional(f,'Dim',defDim,@isnumeric)
addOptional(f,'MinTrackLength',defMinTrackLength,@isnumeric)
addOptional(f,'PrintTrackProgress',defPrintTrackProgress,@isnumeric)
addOptional(f,'MaxDisp',defMaxDisp,@isnumeric)
addOptional(f,'FrameStart',defp.FrameStart,@isnumeric)
addOptional(f,'PlotOpt',defPlotOpt,checkPlotOpt)
addOptional(f,'ErrorThresh',defErrorThresh,@isnumeric)
addOptional(f,'NTests',defNTests,@isnumeric)
addOptional(f,'StepAmplitude',defStepAmplitude,@isnumeric)
% Parse the values from f and put results in p.
parse(f,FileStub,varargin{:})
p = f.Results;

% Get info about data.
info = imfinfo([FileStub '.tif']);
FrameWidth = info.Width;
FrameHeight = info.Height;
NFrames = numel(info);

% Set up parameters for pre-tracking.
param.mem = p.TrackMem;
param.dim = p.Dim;
param.good = p.MinTrackLength;
param.quiet = p.PrintTrackProgress;

% Set parameter for MSD calculation
CollectiveMotionFlag = 0; % 1 = subtract collective motion;
                          % 0 = leave collective motion HARDCODED OPTION

%% Particle Tracking

% Set up arrays
Data = zeros(FrameHeight,FrameWidth,NFrames);
bpData = Data;

% Read in data + bandpasfilter
disp([char(10) 'Loading and bandpassing frames... '])

for Frame = p.FrameStart:p.NFrames
    Data(:,:,Frame-p.FrameStart+1) = double(imread([FileStub '.tif'],Frame));
    bpData(:,:,Frame-p.FrameStart+1) = bpass2D_TA(Data(:,:,Frame-p.FrameStart+1), ...
                                                  p.NoiseSz,p.FeatSize);
end

% Do traditional tracking to determine averaged particle centers
disp([char(10) 'Pretracking... '])
Threshold = max(bpData(:))/p.ThreshFact;
Centers = zeros(0,5);

for Frame = 1:size(bpData,3)
    Peaks = pkfnd(bpData(:,:,Frame),Threshold,p.FeatSize);
    temp = cntrd(bpData(:,:,Frame),Peaks,p.FeatSize,0);
    Centers = [Centers; [temp repmat(Frame,[size(temp,1) 1])]];
end

Tracks = trackin(Centers,p.MaxDisp,param);
NParticles = max(Tracks(:,6));
NTracks = size(unique(Tracks(:,6)),1);

disp([char(9) 'Found a total of ' num2str(NParticles) ' particles' ...
    ' and ' num2str(NTracks) ' tracks.' ]);

% Visually check tracks if desired.
switch p.PlotOpt
    case 'simple'
        disp([char(9) 'Visual check of tracks.'])
        PlotPretracking(Data,bpData,Tracks,p.FeatSize,p.NmPerPixel,FileStub,'simple')
    case 'bandpass'
        disp([char(9) 'Visual check of tracks.'])
        PlotPretracking(Data,bpData,Tracks,p.FeatSize,p.NmPerPixel,FileStub,'bandpass')
    otherwise
        disp([char(9) 'No visual check. If desired use PlotOpt.'])
end

% Compute averaged centers to use as reference points for rest of analysis
disp([char(9) 'Find reference points from pretracking data.'])
refCenters = zeros(size(unique(Tracks(:,6)),1),3);

for ParticleID = 1:max(Tracks(:,6))
    if sum(Tracks(:,6)==ParticleID)~=0
        refCenters(ParticleID,:) = [mean(Tracks(Tracks(:,6)==ParticleID,1:2),1) ParticleID];
    end
end

% Compute noise and estimate centroiding error
disp([char(9) 'Find single particle calibration parameters.'])
[CalibParams,trialCenters] = mserror_calculator_4QM(bpData,Tracks,p.FeatSize, ...
                                                    p.DeltaFit,p.StepAmplitude, ...
                                                    refCenters,p.PlotOpt, ...
                                                    p.ErrorThresh, NParticles, ...
                                                    p.NTests); 
                                        
rmserror = sqrt((CalibParams(:,3) + CalibParams(:,6)));

%% Use single particle calibrations with 4QM to process real data
disp([char(10) '4QM ... '])
disp([char(9) 'Processing real data.'])
QMTracks = QMtrackcorrection(Tracks,bpData,refCenters,CalibParams, ...
                             p.FeatSize,p.DeltaFit);

% For now we will just use the good tracks, can thing about if we want to
% do something different later.
QMGood = QMTracks(QMTracks(:,5)==1,:);
QMGood = QMGood(:,1:4); %resize for use with calcMSD

% Calculate MSDs and errors

disp([char(9) 'Calculating MSDs.'])
MSDs = calcMSD(QMGood,p.NmPerPixel,CollectiveMotionFlag);
disp([char(10) 'Error correction of MSDs ... '])

correctedMSDs = MSDs(:,3:end)-2*repmat(rmserror',size(MSDs,1),1).^2*p.NmPerPixel^2;
averagecorrectedMSDs = nansum(correctedMSDs,2)./size(correctedMSDs,2);
stdcorrectedMSDs = std(correctedMSDs,[],2);
correctedMSDs = [averagecorrectedMSDs stdcorrectedMSDs correctedMSDs];
correctedMSDs(1,:) = zeros(1,size(correctedMSDs,2));
%% output

% Showing MSD graph
switch p.PlotOpt
    case {'simple','bandpass'}
        fig3 = figure();
        whitebg(fig3,[1,1,1])
        loglog(0:size(MSDs,1)-1,correctedMSDs(:,1),'.','color','k')
        hold on
        ylim([min(correctedMSDs(2:end,1))*0.1,max(correctedMSDs(:,1))*10])
        title('Averaged MSD');
        xlabel('\tau (s)');
        ylabel('corrected \langledR^{2}\rangle (nm^{2})');
        savefig([FileStub '_correctedmsd'])
end

% Writing files.
disp([char(10) 'Writing output files...'])
disp([char(9) 'Writing MSD file.'])
csvwrite([FileStub '_msd.csv'],MSDs);
disp([char(9) 'Write rms error file.'])
csvwrite([FileStub '_rmserror.csv'],rmserror);
disp([char(9) 'Writing data as .mat file.'])
save([FileStub '4QMData'],'Tracks','QMTracks','MSDs','correctedMSDs')
disp([char(9) 'Writing input as .mat file.'])
save([FileStub '4QMInput'],'p')

disp(char(9)); toc
end