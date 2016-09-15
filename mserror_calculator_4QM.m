function res = mserror_calculator_4QM(Data,Tracks,FeatSize,DeltaFit, ...
                                      StepAmplitude,ThreshFact,refCenters)

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

%% setup the various coordinate systems and coordinate shifts used in this code
xsub = 1:(2*(FeatSize-DeltaFit));
ysub = 1:(2*(FeatSize-DeltaFit));
[xsub, ysub] = meshgrid(xsub,ysub);

%% step through particles, find best frame to use.

TrackID = 0;

for ParticleID = 1:max(Tracks(:,6))  
    subTracks = Tracks(Tracks(:,6)==ParticleID,1:2);   
    if ~isempty(subTracks)
        TrackID = TrackID + 1;
        metricDistance = sqrt((subTracks(:,1)-refCenters(TrackID,1)).^2 + (subTracks(:,2)-refCenters(TrackID,2)).^2);
        [~, refStep] = min(metricDistance);
        
        frames = Tracks(Tracks(:,6)==ParticleID,5);
        xCoarse = refCenters(TrackID,1);
        yCoarse = refCenters(TrackID,2);
        Cols = SetAxisSubdata(xCoarse,FeatSize,DeltaFit);
        Rows = SetAxisSubdata(yCoarse,FeatSize,DeltaFit);
        subdata = Data(Rows,Cols,frames);
         
        %     first do calibration using shiftedData
        ref_subData = subdata(:,:,refStep);
        shiftedData = zeros([size(xsub),numel(frames)]);
        
        fake_dx = StepAmplitude*(randn(numel(frames),1));
        fake_dy = StepAmplitude*(randn(numel(frames),1));
        
        for j = 1:numel(frames)
            shiftedx = xsub + fake_dx(j);
            shiftedy = ysub + fake_dy(j);
            shiftedData(:,:,j) = interp2(xsub,ysub,ref_subData,shiftedx,shiftedy);
        end

% Make sure noise in calibration data represents noise in real data

%         mean_noise = mean(subdata(subdata(:)<(max(subdata(:))/ThreshFact))); % target mean in noise
%         std_noise = std(subdata(subdata(:)<(max(subdata(:))/ThreshFact))); % target std in noise
%         max_noise = max(abs(subdata(subdata(:)<(max(subdata(:))/ThreshFact))-mean_noise));
%         SNR = mean(subdata(round(size(subdata,1)/2),round(size(subdata,2)/2),:))/mean_noise; %approximate signal to noise ratio
%         temp_noise = 2*(rand([size(shifted_data)])-0.5)*std_noise/2 + mean_noise;

        shiftedData(isnan(shiftedData(:))) = subdata(isnan(shiftedData(:)));
        
%         mean_shifted = mean(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % mean in shifted noise
%         std_shifted = std(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % std in shifted noise
%         normdata = (shifted_data - mean_shifted)/std_shifted;        
%         added_noise = 2*(rand([size(normdata)])-0.5)/(0.5*SNR);       
%         noisydata = normdata + added_noise;        
%         normnoisydata = (noisydata - mean(noisydata(noisydata(:)<(max(noisydata(:))/threshfact))))/std(noisydata(noisydata(:)<(max(noisydata(:))/threshfact)));
%         scalednoisydata = normnoisydata*std_noise + mean_noise;        
%         calib_params(ptcle,:) = [pkcnt_4QM_calibrator_0(scalednoisydata,tracks,fake_dx,fake_dy,feat_size,delta_fit) ptcle];
%         for frame = 1:numel(frames)
%             imagesc([subdata(:,:,frame) shifted_data(:,:,frame) scalednoisydata(:,:,frame)])
%             getframe
%         end
        calib_params(ParticleID,:) = [FQM(shiftedData,fake_dx,fake_dy,1,[],1) ParticleID];
        
        

        
    end
end

res = calib_params;

end