function [] = tracker_caller_4QM(FileStub,NFrames,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
% SEGMENTATION AND TRACKING OF PARTICLES VIA THE 4QM METHOD IN 2D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% FileStub -- Path to the image_stack to be analyzed. 
% NFrames -- Number of frames in the image stack.
% set_length -- Number of images in a set. 
% NmPerPixel -- [optional] Actual pixel width (nm).
% SecsPerFrame -- [optional] Time between frames (sec).
% NoiseSz -- [optional] (pixels).
% FeatSize -- [optional] Full optical diameter of particle (pixels).
% DeltaFit -- [optional] Narrows analysis region around particle (pixels).
% treshfact -- [optional] maximum intensity devided by the thresfact gives 
% the threshold value.
% TrackMem -- [optional] Number of steps disconnected tracks can be 
% reconnected, in case a particle is lost.
% Dim -- [optional] Dimension of the system.
% MinTrackLength -- [optional] minimum length of track; throw away tracks 
% shorter than this.
% PrintTrackProgress -- [optional] Turns on or off printing progress to screen.
% MaxDisp -- [optional] MaxDisp should be set to a value somewhat less than 
% the mean spacing between the particles.
% PlotOpt -- Options for plotting data ['simple','bandpassed',0] 
%
% NOTES
%
% The imwrite() function is unstable when windows file explorer is opened.
%
% DEPENDENCIES
%
% This program depends on the particle_tracking toolbox from the TA-lab for
% the following functions: bpass2D_TA() and msd_manual2().
% This program depends on the SPtrack1.0 toolbox by Eric Dufresne from Yale 
% University for the following functions: pkfnd() and cntrd(). 
% This program depends on the trackin() function from Crocker 
% (http://glinda.lrsm.upenn.edu/~weeks/idl). 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
f = inputParser;
f.CaseSensitive = 0;
%% Set up particle, intensity and duration parameters.

% Default Values
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
defp.FrmStart = 1;
defPlotOpt = 'none';
validPlotOpt = {'bandpass','simple','none'};
checkPlotOpt = @(x) any(validatestring(x,validPlotOpt));

% Set up required and optional inputs
addRequired(f,'FileStub',@ischar)
addRequired(f,'NFrames',@isnumeric)
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
addOptional(f,'FrmStart',defp.FrmStart,@isnumeric)
addOptional(f,'PlotOpt',defPlotOpt,checkPlotOpt)

% Parse the values to p
parse(f,FileStub,NFrames,varargin{:})
p = f.Results;

% Set up parameters for pre-tracking.
param.mem = p.TrackMem;
param.dim = p.Dim;
param.good = p.MinTrackLength;
param.quiet = p.PrintTrackProgress;

% Set parameters for error calculation.
step_amplitude = 1;

%% Particle Tracking 
    
    % Set up arrays   
    temp = double(imread([FileStub '.tif'],p.FrmStart));
    data = zeros(size(temp,1),size(temp,2),NFrames);
    b = data;
    
    % Read in data + bandpasfilter
    disp([char(10) 'Loading and bandpassing frames... '])
        for frame = p.FrmStart:NFrames
        data(:,:,frame-p.FrmStart+1) = double(imread([FileStub '.tif'],frame));
        b(:,:,frame-p.FrmStart+1) = bpass2D_TA(data(:,:,frame-p.FrmStart+1), ...
                                             p.NoiseSz,p.FeatSize);
        end
    
    % Do traditional tracking to determine averaged particle centers
    disp([char(10) 'Pretracking... '])    
    thresh = max(b(:))/p.ThreshFact;
    cnt = zeros(0,5);
    
    for frame = 1:size(b,3)
        pk = pkfnd(b(:,:,frame),thresh,p.FeatSize);  
        temp = cntrd(b(:,:,frame),pk,p.FeatSize,0);
        cnt = [cnt; [temp repmat(frame,[size(temp,1) 1])]];
    end
  
    tracks = trackin(cnt,p.MaxDisp,param);
    Ntracks = size(unique(tracks(:,6)));
    disp([char(9) 'Found a total of ' num2str(Ntracks(1)) ' tracks.'])
    clear cnt

    % Visually check tracks if desired.
    switch p.PlotOpt
        case 'simple'
            disp([char(9) 'Visual check of tracks.'])
            PlotPretracking(data,b,tracks,p.FeatSize,p.NmPerPixel,FileStub,'simple')
        case 'bandpass'
            disp([char(9) 'Visual check of tracks.'])
            PlotPretracking(data,b,tracks,p.FeatSize,p.NmPerPixel,FileStub,'bandpass')            
        otherwise
            disp([char(9) 'No visual check. If desired use PlotOpt.'])              
    end
                
    % Compute averaged centers to use as reference points for rest of analysis
    disp([char(9) 'Find reference points from pretracking data.'])  
    ptclecnt = 0;
    
    for ptcle = 1:max(tracks(:,6))
        ptclecnt = ptclecnt + 1;
        
        if sum(tracks(:,6)==ptcle)~=0
            ref_cnts(ptclecnt,:) = [mean(tracks(tracks(:,6)==ptcle,1:2),1) ptcle];
        end
    end
    
    % Compute noise and estimate centroiding error
    disp([char(9) 'Find single particle calibration parameters.'])
    calibration_params = mserror_calculator_4QM(b,tracks,p.FeatSize, ...
                                                p.DeltaFit,step_amplitude, ...
                                                p.ThreshFact,ref_cnts); 
    rmserror = sqrt((calibration_params(:,3) + calibration_params(:,6)));
    mean(rmserror);
    
    %% Now use single particle calibrations with 4QM to process real data
    disp([char(10) '4QM ... '])
    disp([char(9) 'Processing real data.'])
    tracks_4QM = zeros(0,4);
    
    for particle = 1:max(tracks(:,6))       
        frames = tracks(tracks(:,6)==particle,5);
        x_coarse = ref_cnts(particle,1);
        y_coarse = ref_cnts(particle,2);      
        cols = SetAxisSubdata(x_coarse,p.FeatSize,p.DeltaFit);
        rows = SetAxisSubdata(y_coarse,p.FeatSize,p.DeltaFit);      
        subdata = b(rows,cols,frames);
        
        tracks_4QM = [tracks_4QM; [x_coarse*ones(numel(frames),1), ...
                      y_coarse*ones(numel(frames),1), zeros(numel(frames),1)] ...
                      + FQM(subdata,[],[],0,calibration_params(particle,:),1) ...
                      particle*ones(numel(frames),1)];
    
    end
    
    % Calculate MSD's and write to CSV files per set
    disp([char(9) 'Calculating MSDs.'])
    collective_motion_flag = 0; % 1 = subtract collective motion; 0 = leave collective motion HARDCODED OPTION
    msd_temp = calcMSD(tracks_4QM,p.NmPerPixel,collective_motion_flag);

    disp([char(9) 'Writing MSD file.'])
    csvwrite([FileStub '_msd.csv'],msd_temp);
    disp([char(9) 'Write rms error file.'])
    csvwrite([FileStub '_rmserror.csv'],rmserror);
    
    if ~strcmp(p.PlotOpt,'none')
        loglog((0:size(msd_temp,1)-1),msd_temp(:,1)-2*mean(rmserror)^2,'.') % These are the blue shapes
        hold on
        getframe;
    end
    
    disp(char(9))
    toc

corrected_SPmsds = zeros(size(msd_temp,1),0);
disp([char(10) 'Error correction of MSDs ... '])
    
% msds = csvread([FileStub '_msd.csv']);
% rmserrors = csvread([FileStub '_rmserror.csv']);
% 
% full_error = repmat(rmserrors',[size(msds,1) 1]);

corrected_SPmsds = msd_temp(:,3:end)-2*repmat(rmserror',size(msd_temp,1),1).^2*p.NmPerPixel^2;

% corrected_AVEmsds = [mean(corrected_SPmsds(2:end,:),1)' ...
%                      std(corrected_SPmsds(2:end,:),[],1)'];

% disp([char(9) 'Write corrected MSD file.'])
% csvwrite([FileStub '_corrected_msd.csv'],msd_temp);
% disp([char(9) 'Write corrected rms error file.'])
% csvwrite([FileStub '_corrected_rmserror.csv'],rmserror);
% 
% final_AVEmsds = nanmean(corrected_AVEmsds,1);
% final_AVEmsds(:,2) = final_AVEmsds(:,2)/sqrt(size(corrected_AVEmsds,1));
% 
% disp([char(9) 'Write error corrected msd file.'])    
% csvwrite([FileStub '_error_corrected_msd.csv'], ...
%          [mean(corrected_SPmsds,2) std(corrected_SPmsds,[],2) corrected_SPmsds]);
% 
% %corrected_SPmsds(:,1)=[];
% 
% corrected_SPmsds = [(0:(size(corrected_SPmsds,1)-1))'*SecsPerFrame corrected_SPmsds*p.NmPerPixel^2];