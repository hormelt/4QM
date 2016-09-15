function res = mserror_calculator_4QM(Data,Tracks,FeatSize,DeltaFit, ...
                                      StepAmplitude,ThreshFact,refCenters,PlotOpt)

% Calculates the mean squared error in the particle position upon subpixel
% displacements.
%
% INPUTS:
%   Data: Collection of image frames.
%   Tracks: The collection of particle tracks in the frames.
%   FeatSize: Full optical diameter of particle (pixels).
%   DeltaFit: Narrows analysis region around particle (pixels).
%   StepAmplitude: Maximum Amplitude of shift.
%   ntests: Number of shift tests.
%   refCenters: The value used as a reference for the particle centers.
%
% OUTPUTS:
%   res: The calibration parameters. [p1 errx p2 erry].
%
% CODE FOR FUTURE:
%   Currently the amplitude of the trial shifts is determined by
%   StepAmplitude. It is best to set the StepAmplitude in the trial shifts 
%   equal to the amplitude of the shifts in the real data.
%   In the future we would like to determine the shift amplitude from the 
%   real data. The code snippets below could be used for that. Also there
%   are some ideas about adding noise?
%
%   % Make sure noise in calibration data represents noise in real data
%   mean_noise = mean(subdata(subdata(:)<(max(subdata(:))/ThreshFact))); % target mean in noise
%   std_noise = std(subdata(subdata(:)<(max(subdata(:))/ThreshFact))); % target std in noise
%   max_noise = max(abs(subdata(subdata(:)<(max(subdata(:))/ThreshFact))-mean_noise));
%   SNR = mean(subdata(round(size(subdata,1)/2),round(size(subdata,2)/2),:))/mean_noise; %approximate signal to noise ratio
%   temp_noise = 2*(rand([size(shifted_data)])-0.5)*std_noise/2 + mean_noise;

%   mean_shifted = mean(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % mean in shifted noise
%   std_shifted = std(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % std in shifted noise
%   normdata = (shifted_data - mean_shifted)/std_shifted;        
%   added_noise = 2*(rand([size(normdata)])-0.5)/(0.5*SNR);       
%   noisydata = normdata + added_noise;        
%   normnoisydata = (noisydata - mean(noisydata(noisydata(:)<(max(noisydata(:))/threshfact))))/std(noisydata(noisydata(:)<(max(noisydata(:))/threshfact)));
%   scalednoisydata = normnoisydata*std_noise + mean_noise;        
%   calib_params(ptcle,:) = [pkcnt_4QM_calibrator_0(scalednoisydata,tracks,fake_dx,fake_dy,feat_size,delta_fit) ptcle];
%   for frame = 1:NFrames
%       imagesc([subdata(:,:,frame) shifted_data(:,:,frame) scalednoisydata(:,:,frame)])
%       getframe
%   end


%% Setup pixel grid for subdata frames.
xGrid = 1:(2*(FeatSize-DeltaFit));
yGrid = 1:(2*(FeatSize-DeltaFit));
[xGrid, yGrid] = meshgrid(xGrid,yGrid);

%% Find calibration parameters for every particle with tracks.
TrackID = 0;

for ParticleID = 1:max(Tracks(:,6))  
    subTracks = Tracks(Tracks(:,6)==ParticleID,:);
    
    if ~isempty(subTracks)
        TrackID = TrackID + 1;
        NFrames = numel(subTracks(:,5));
        
        % Determine the subData frame around the particle of interest.
        xCoarse = refCenters(TrackID,1);
        yCoarse = refCenters(TrackID,2);
        Cols = SetAxisSubdata(xCoarse,FeatSize,DeltaFit);
        Rows = SetAxisSubdata(yCoarse,FeatSize,DeltaFit);
        subData = Data(Rows,Cols,subTracks(:,5));
         
        % Find the frame in which the particle is closest to its refCenter.
        metricDistance = sqrt((subTracks(:,1)-refCenters(TrackID,1)).^2 ...
                         + (subTracks(:,2)-refCenters(TrackID,2)).^2);
        [~, refStep] = min(metricDistance);
        refsubData = subData(:,:,refStep);

        % Determine the trial shifts of the refFrame for 4QM calibration.
        dxTrialShift = StepAmplitude*(randn(NFrames,1));
        dyTrialShift = StepAmplitude*(randn(NFrames,1));
        
        % Create the trial data by shifting refsubData by TrialShift.
        shiftedData = zeros([size(xGrid),NFrames]);
        
        for j = 1:NFrames
            shiftedxGrid = xGrid + dxTrialShift(j);
            shiftedyGrid = yGrid + dyTrialShift(j);
            shiftedData(:,:,j) = interp2(xGrid,yGrid,refsubData, ...
                                         shiftedxGrid,shiftedyGrid);
        end
        
        
        % Replace NaN resulting from interp2 by data from the subrefFrame
        shiftedData(isnan(shiftedData(:))) = subData(isnan(shiftedData(:)));
        
        
        % Start Calibration
        [A,B,C,D] = FQM(shiftedData);
        cnt = [(A+C-B-D)./(A+B+C+D) (A+B-C-D)./(A+B+C+D)];
        refShift = [dxTrialShift dyTrialShift];
    
        % now find the shift that gives the smallest error
    
        [p1,fvalx] = fminsearch(@(p1) squeeze(mean((p1(1)*(cnt(:,1)+p1(2))-refShift(:,1)).^2,1)),...
                     [range(refShift(:,1))/range(cnt(:,1)),mean(refShift(:,1))]);
        errx = sqrt(fvalx);
    
        [p2,fvaly] = fminsearch(@(p2) squeeze(mean((p2(1)*(cnt(:,2)+p2(2))-refShift(:,2)).^2,1)),...
                     [range(refShift(:,2))/range(cnt(:,2)),mean(refShift(:,2))]);
        erry = sqrt(fvaly);
    
        %     if (errx<=1e-1) && (erry<=1e-1)
    
        %csvwrite([num2str(round(rand*1000)) '.csv'],[p1(1)*(cnt(:,1)+p1(2)),refshft(:,1) p2(1)*(cnt(:,2)+p2(2)),refshft(:,2)]);
    
        res1 = [p1 errx p2 erry];
    
        if PlotOpt
            scatter(p1(1)*(cnt(:,1)+p1(2)),refShift(:,1),'b')
            hold on
            scatter(p2(1)*(cnt(:,2)+p2(2)),refShift(:,2),'g')
            getframe;  
        end      
        
        CalibParams(ParticleID,:) = [res1 ParticleID];
    end
end

res = CalibParams;

end