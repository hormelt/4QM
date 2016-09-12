function [] = tracker_caller_4QM(filestub,nframes,set_length, ...
                                 nm_per_pixel,secs_per_frame,noise_sz, ...
                                 feat_size,delta_fit,threshfact, ...
                                 TrackMem,Dim,MinTrackLength, ...
                                 PrintTrackProgress, maxdisp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                             
% SEGMENTATION AND TRACKING OF PARTICLES VIA THE 4QM METHOD IN 2D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% filestub -- Path to the image_stack to be analyzed. 
% nframes -- Number of frames in the image stack.
% set_length -- Number of images in a set. 
% nm_per_pixel -- [optional] Actual pixel width (nm).
% secs_per_frame -- [optional] Time between frames (sec).
% noise_sz -- [optional] (pixels).
% feat_size -- [optional] Full optical diameter of particle (pixels).
% delta_fit -- [optional] Narrows analysis region around particle (pixels).
% treshfact -- [optional] maximum intensity devided by the thresfact gives 
% the threshold value.
% TrackMem -- [optional] Number of steps disconnected tracks can be 
% reconnected, in case a particle is lost.
% Dim -- [optional] Dimension of the system.
% MinTrackLength -- [optional] minimum length of track; throw away tracks 
% shorter than this.
% PrintTrackProgress -- [optional] Turns on or off printing progress to screen.
% maxdisp -- [optional] maxdisp should be set to a value somewhat less than 
% the mean spacing between the particles.
 
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

%% Set up particle, intensity and duration parameters.
if nargin < 4, nm_per_pixel = 16.25; end  % zyla at 400x
if nargin < 5, secs_per_frame = 0.011179722; end
if nargin < 6, noise_sz = 1; end
if nargin < 7, feat_size = 15; end
if nargin < 8, delta_fit = 3; end
if nargin < 9, threshfact = 3.5; end
if nargin < 10, TrackMem = 0; end
if nargin < 11, Dim = 2; end
if nargin < 12, MinTrackLength = set_length; end
if nargin < 13, PrintTrackProgress = 1; end
if nargin < 14, maxdisp = feat_size/2; end

nsets = floor(nframes/set_length);

% Set up parameters for pre-tracking.
param.mem = TrackMem;
param.dim = Dim;
param.good = MinTrackLength;
param.quiet = PrintTrackProgress;

% Set parameters for error calculation.
step_amplitude = 1;
ntests = 100;

%% Step through sets

for set = 1:nsets   
    frmstart = (set-1)*set_length + 1;
    frmend = frmstart + set_length - 1;
    
    % Set up arrays   
    temp = double(imread([filestub '.tif'],frmstart));
    data = zeros(size(temp,1),size(temp,2),frmend-frmstart+1);
    b = data;
    
    % Read in data + bandpasfilter
    disp([char(10) 'Loading and bandpassing frames... '])
        for frame = frmstart:frmend
        data(:,:,frame-frmstart+1) = double(imread([filestub '.tif'],frame));
        b(:,:,frame-frmstart+1) = bpass2D_TA(data(:,:,frame-frmstart+1) ...
                                             ,noise_sz,feat_size);
        
        % Put out progress of reading process                                                                          
        if mod(frame-set_length*(set-1),50) == 0 || frame == set_length*set 
           disp([char(9) num2str(frame-set_length*(set-1), '%3d') ' of ' ...
                 num2str(set_length) ' done.']);
        end
        
    end
    
    % Do traditional tracking to determine averaged particle centers
    disp([char(10) 'Pretracking... '])    
    thresh = max(b(:))/threshfact;
    cnt = zeros(0,5);
    
    for frame = 1:size(b,3)
        pk = pkfnd(b(:,:,frame),thresh,feat_size);  
        temp = cntrd(b(:,:,frame),pk,feat_size,0);
        cnt = [cnt; [temp repmat(frame,[size(temp,1) 1])]];
    end
  
    tracks = trackin(cnt,maxdisp,param);
    Ntracks = size(unique(tracks(:,6)));
    disp([char(9) 'Found a total of ' num2str(Ntracks(1)) ' tracks for set ' num2str(set) '.'])
    clear cnt
    
                % Visually check tracks
                disp([char(9) 'Visual check of tracks for set ' num2str(set) '. '])
                for frame = 1:size(data,3)
    
                    tempx = tracks(tracks(:,5)==frame,1);
                    tempy = tracks(tracks(:,5)==frame,2);
    
                    hold off
                    imagesc(data(:,:,frame))
                    colormap gray
                    hold on
                    scatter(tempx,tempy,'r')
                    truesize
                    f = getframe;
                    imwrite(frame2im(f),[filestub 'tracking_movie.tif'], ...
                            'tiff','compression','none','writemode','append');
    
                end
                close
                
    % Compute averaged centers to use as reference points for rest of analysis
    disp([char(9) 'Find reference points from pretracking data for set ' num2str(set) '.'])  
    ptclecnt = 0;
    
    for ptcle = 1:max(tracks(:,6))
        ptclecnt = ptclecnt + 1;
        
        if sum(tracks(:,6)==ptcle)~=0
            ref_cnts(ptclecnt,:) = [mean(tracks(tracks(:,6)==ptcle,1:2),1) ptcle];
        end
    end
    
    % Compute noise and estimate centroiding error
    disp([char(9) 'Find single particle calibration parameters for set ' num2str(set) '.'])
    calibration_params = mserror_calculator_4QM(b,tracks,feat_size, ...
                                                delta_fit,step_amplitude, ...
                                                ntests,threshfact,ref_cnts);
    rmserror = sqrt((calibration_params(:,3) + calibration_params(:,6)));
    mean(rmserror);
    
    %% Now use single particle calibrations with 4QM to process real data
    disp([char(10) '4QM ... '])
    disp([char(9) 'Process real data for set ' num2str(set) '.'])
    tracks_4QM = zeros(0,4);
    
    for particle = 1:max(tracks(:,6))
        
        frames = tracks(tracks(:,6)==ptcle,5);
        x_coarse = ref_cnts(ptclecnt,1);
        y_coarse = ref_cnts(ptclecnt,2);
        
        if round(x_coarse) > x_coarse
            cols = (round(x_coarse)-(feat_size-delta_fit)): ...
                   (round(x_coarse)+(feat_size-delta_fit))-1;
        else
            cols = (round(x_coarse)-(feat_size-delta_fit))+1: ...
                   (round(x_coarse)+(feat_size-delta_fit));
        end
        
        if round(y_coarse) > y_coarse
            rows = (round(y_coarse)-(feat_size-delta_fit)): ...
                   (round(y_coarse)+(feat_size-delta_fit))-1;
        else
           rows = (round(y_coarse)-(feat_size-delta_fit))+1: ...
                  (round(y_coarse)+(feat_size-delta_fit));
        end
        
        subdata = b(rows,cols,frames);
        
        tracks_4QM = [tracks_4QM; [[[x_coarse*ones(frames(end),1), ...
                      y_coarse*ones(frames(end),1), zeros(frames(end),1)] ...
                      + FQM(subdata,[],[],0,calibration_params(particle,:),1)] ...
                      particle*ones(frames(end),1)]];
    
    end
    
    % Calculate MSD's and write to CSV files per set
    disp([char(9) 'Calculate MSDs for set ' num2str(set) '.'])
    collective_motion_flag = 0; % 1 = subtract collective motion; 0 = leave collective motion
    msd_temp = msd_manual2(tracks_4QM,nm_per_pixel,collective_motion_flag);%-2*(rmserror^2)*nm_per_pixel^2);

    disp([char(9) 'Write MSD file for set ' num2str(set) '.'])
    csvwrite([filestub '_set' num2str(set) '_msd.csv'],msd_temp);
    disp([char(9) 'Write rms error file for set ' num2str(set) '.'])
    csvwrite([filestub '_set' num2str(set) '_rmserror.csv'],rmserror);
    
    loglog((0:size(msd_temp,1)-1),msd_temp(:,1)-2*mean(rmserror)^2,'.') % These are the blue shapes
    hold on
    getframe;
    
    disp(char(9))
    toc
    
end

corrected_SPmsds = zeros(size(msd_temp,1),0);
disp([char(10) 'Error corecction of MSDs ... '])
for set = 1:nsets
    
    msds = csvread([filestub '_set' num2str(set) '_msd.csv']);
    rmserrors = csvread([filestub '_set' num2str(set) '_rmserror.csv']);
    
    full_error = repmat(rmserrors',[size(msds,1) 1]);
    
    corrected_SPmsds = [corrected_SPmsds msds(:,3:end)-4*full_error];
    corrected_AVEmsds = [mean(corrected_SPmsds(2:end,:),1)' ...
                         std(corrected_SPmsds(2:end,:),[],1)'];

    disp([char(9) 'Write corrected MSD file for set ' num2str(set) '.'])
    csvwrite([filestub '_set' num2str(set) '_corrected_msd.csv'],msd_temp);
    disp([char(9) 'Write corrected rms error file for set ' num2str(set) '.'])
    csvwrite([filestub '_set' num2str(set) '_corrected_rmserror.csv'],rmserror);
    
    final_AVEmsds(set,:) = mean(corrected_AVEmsds,1);
    final_AVEmsds(set,2) = final_AVEmsds(set,2)/sqrt(size(corrected_AVEmsds,1));
    
    disp([char(9) 'Write error corrected msd file for set ' num2str(set) '.'])    
    csvwrite([filestub '_error_corrected_msd.csv'], ...
             [mean(corrected_SPmsds,2) std(corrected_SPmsds,[],2) corrected_SPmsds]);
    
end

%corrected_SPmsds(:,1)=[];

corrected_SPmsds = [(0:(size(corrected_SPmsds,1)-1))'*secs_per_frame corrected_SPmsds*nm_per_pixel^2];