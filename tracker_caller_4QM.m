function [] = tracker_caller_4QM(filestub,nframes,set_length, ...
                                 nm_per_pixel,secs_per_frame,noise_sz, ...
                                 feat_size,delta_fit,threshfact)

% SEGMENTATION AND TRACKING OF PARTICLES VIA THE 4QM METHOD IN 2D.
%
% INPUTS
%
% filestub -- Path to the image_stack to be analyzed. 
% nframes -- Number of frames in the image stack.
% set_length -- Number of images in a set. 
% nm_per_pixel -- [optional] Actual pixel width (nm).
% secs_per_frame -- [optional] Time between frames (sec).
% noise_sz -- [optional] (pixels).
% feat_size -- [optional] Full optical radius of particle (pixels).
% delta_fit -- [optional] Narrows analysis region around particle (pixels).
% treshfact -- [optional] maximum intensity devided by the thresfact gives
%              the threshold value

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up particle, intensity and duration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9, threshfact = 3.5; end
if nargin < 8, delta_fit = 3; end
if nargin < 7, feat_size = 15; end
if nargin < 6, noise_sz = 1; end
if nargin < 5, secs_per_frame = 0.011179722; end
if nargin < 4, nm_per_pixel = 16.25; end  % zyla at 400x

nsets = floor(nframes/set_length);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step through sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for set = 1:nsets   
    frmstart = (set-1)*set_length + 1;
    frmend = frmstart + set_length - 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up arrays
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    temp = double(imread([filestub '.tif'],frmstart));
    data = zeros(size(temp,1),size(temp,2),frmend-frmstart+1);
    b = data;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read in data + bandpasfilter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for frame = frmstart:frmend
        data(:,:,frame-frmstart+1) = double(imread([filestub '.tif'],frame));
        b(:,:,frame-frmstart+1) = bpass2D_TA(data(:,:,frame-frmstart+1) ...
                                                  ,noise_sz,feat_size)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do traditional tracking to determine averaged particle centers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    thresh = max(b(:))/threshfact;
    cnt = zeros(1,5);
    
    for frame = 1:size(b,3)
        pk = pkfnd(b(:,:,frame),thresh,2*feat_size);
        size(pk)
        temp = cntrd(b(:,:,frame),pk,2*feat_size,0);
        size(temp)
        cnt = [cnt; [temp repmat(frame,[size(temp,1) 1])]];
    end
    
    cnt(1,:)=[];
    
    param.mem = 0; %number of steps disconnected tracks can be reconnected,in case a particle is lost
    param.dim = 2; %dimension of the system
    param.good = size(data,3); %minimum length of track; throw away tracks shorter than this
    param.quiet = 0; %turns on or off printing progress to screen
    maxdisp = feat_size/2; %maxdisp should be set to a value somewhat less than the mean spacing between the particles.
    
    tracks = trackin(cnt,maxdisp,param);
    clear cnt
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Visually check tracks
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute averaged centers to use a reference points for rest of
    % analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ptclecnt = 0;
    
    for ptcle = 1:max(tracks(:,6))
        ptclecnt = ptclecnt + 1;
        
        if sum(tracks(:,6)==ptcle)~=0
            ref_cnts(ptclecnt,:) = [mean(tracks(tracks(:,6)==ptcle,1:2),1) ptcle];
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute noise and estimate centroiding error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    param.feat_size = feat_size;
    param.delta_fit = delta_fit;
    param.step_amplitude = 1;
    param.ntests = 100;
    param.threshfact = threshfact;
    param.ref_cnts = ref_cnts;
    
    calibration_params = mserror_calculator_4QM(b,tracks,param);
    rmserror = sqrt((calibration_params(:,3) + calibration_params(:,6)));
    mean(rmserror);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now use single particle calibrations with 4QM to process real data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tracks_4QM = zeros(1,4);
    
    for particle = 1:max(tracks(:,6))
        
        frames = tracks(tracks(:,6)==ptcle,5);
        x_coarse = ref_cnts(ptclecnt,1);
        y_coarse = ref_cnts(ptclecnt,2);
        
        if round(x_coarse) > x_coarse
            cols = (round(x_coarse)-(feat_size-delta_fit)):(round(x_coarse)+(feat_size-delta_fit))-1;
        else
            cols = (round(x_coarse)-(feat_size-delta_fit))+1:(round(x_coarse)+(feat_size-delta_fit));
        end
        
        if round(y_coarse) > y_coarse
            rows = (round(y_coarse)-(feat_size-delta_fit)):(round(y_coarse)+(feat_size-delta_fit))-1;
        else
           rows = (round(y_coarse)-(feat_size-delta_fit))+1:(round(y_coarse)+(feat_size-delta_fit));
        end
        
        subdata = b(rows,cols,frames);
        
        tracks_4QM = [tracks_4QM; [[[x_coarse*ones(frames(end),1), y_coarse*ones(frames(end),1), zeros(frames(end),1)]...
            + FQM(subdata,[],[],0,calibration_params(particle,:))] particle*ones(frames(end),1)]];
    
    end
    
    tracks_4QM(1,:) = [];
    
    collective_motion_flag = 0; % 1 = subtract collective motion; 0 = leave collective motion
    msd_temp = msd_manual2(tracks_4QM,nm_per_pixel,collective_motion_flag);%-2*(rmserror^2)*nm_per_pixel^2);
    
    csvwrite([filestub '_set' num2str(set) '_msd.csv'],msd_temp);
    csvwrite([filestub '_set' num2str(set) '_rmserror.csv'],rmserror);
    
    loglog((0:size(msd_temp,1)-1),msd_temp(:,1)-2*mean(rmserror)^2,'.')
    hold on
    getframe
    
    toc
    
end

corrected_SPmsds = zeros(size(msd_temp,1),1);

for set = 1:nsets
    
    msds = csvread([filestub '_set' num2str(set) '_msd.csv']);
    rmserrors = csvread([filestub '_set' num2str(set) '_rmserror.csv']);
    
    full_error = repmat(rmserrors',[size(msds,1) 1]);
    
    corrected_SPmsds = [corrected_SPmsds msds(:,3:end)-4*full_error];
    corrected_AVEmsds = [mean(corrected_SPmsds(2:end,:),1)' std(corrected_SPmsds(2:end,:),[],1)'];
    
    csvwrite([filestub '_set' num2str(set) '_corrected_msd.csv'],msd_temp);
    csvwrite([filestub '_set' num2str(set) '_corrected_rmserror.csv'],rmserror);
    
    final_AVEmsds(set,:) = mean(corrected_AVEmsds,1);
    final_AVEmsds(set,2) = final_AVEmsds(set,2)/sqrt(size(corrected_AVEmsds,1));
    
    csvwrite([filestub '_error_corrected_msd.csv'],[mean(corrected_SPmsds,2) std(corrected_SPmsds,[],2) corrected_SPmsds]);
    
end

corrected_SPmsds(:,1)=[];

corrected_SPmsds = [(0:(size(corrected_SPmsds,1)-1))'*secs_per_frame corrected_SPmsds*nm_per_pixel^2];