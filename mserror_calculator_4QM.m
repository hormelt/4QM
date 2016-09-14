function res = mserror_calculator_4QM(data,tracks,feat_size,delta_fit, ...
<<<<<<< HEAD
    step_amplitude,threshfact, ...
    ref_cnts)
=======
                                      step_amplitude,threshfact, ...
                                      ref_cnts)
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15

% Calculates the mean squared error in the paticle position upon subpixel
% displacements
%
% INPUTS:
<<<<<<< HEAD
%    data: Collection of image frames.
=======
%    data: Collection of image frames. 
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15
%    tracks: The collection of particle tracks in the frames.
%    feat_size: Full optical diameter of particle (pixels).
%    delta_fit: Narrows analysis region around particle (pixels).
%    step_amplitude: Maximum Amplitude of shift.
%    ntests: Number of shift tests.
<<<<<<< HEAD
%    ref_cnts: The value used as a reference for the particle centers.
=======
%    ref_cnts: The value used as a reference for the particle centers. 
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15
%
% OUTPUTS:
%    res: The calibration parameters. [p1 errx p2 erry]

%% setup the various coordinate systems and coordinate shifts used in this code
<<<<<<< HEAD
=======
fake_dx = step_amplitude*(randn(size(data,3),1));
fake_dy = step_amplitude*(randn(size(data,3),1));
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15

% fake_dx = 2*(rand(ntests,1)-0.5)*step_amplitude;
% fake_dy = 2*(rand(ntests,1)-0.5)*step_amplitude;

xsub = 1:(2*(feat_size-delta_fit));
ysub = 1:(2*(feat_size-delta_fit));
[xsub, ysub] = meshgrid(xsub,ysub);

%% step through particles, find best frame to use, make sure noise in
% calibration data represents noise in real data

ptclecnt = 0;

for ptcle = 1:max(tracks(:,6))
    
    subtracks = tracks(tracks(:,6)==ptcle,1:2);
    
    if ~isempty(subtracks)
        ptclecnt = ptclecnt + 1;
        
        distance_metric = sqrt((subtracks(:,1)-ref_cnts(ptclecnt,1)).^2 + (subtracks(:,2)-ref_cnts(ptclecnt,2)).^2);
        [~, ref_step] = min(distance_metric);
        
<<<<<<< HEAD
        frames = tracks(tracks(:,6)==ptcle,5);
        x_coarse = ref_cnts(ptclecnt,1);
        y_coarse = ref_cnts(ptclecnt,2);
        cols = SetAxisSubdata(x_coarse,feat_size,delta_fit);
        rows = SetAxisSubdata(y_coarse,feat_size,delta_fit);
        
        %     first do calibration using shifted_data
        subdata = data(rows,cols,frames);
=======
        frames = tracks(tracks(:,6)==ptcle,5);       
        x_coarse = ref_cnts(ptclecnt,1);
        y_coarse = ref_cnts(ptclecnt,2);        
        cols = SetAxisSubdata(x_coarse,feat_size,delta_fit);
        rows = SetAxisSubdata(y_coarse,feat_size,delta_fit);
        
        %     first do calibration using shifted_data       
        subdata = data(rows,cols,frames);    
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15
        mean_noise = mean(subdata(subdata(:)<(max(subdata(:))/threshfact))); % target mean in noise
        std_noise = std(subdata(subdata(:)<(max(subdata(:))/threshfact))); % target std in noise
        max_noise = max(abs(subdata(subdata(:)<(max(subdata(:))/threshfact))-mean_noise));
        SNR = mean(subdata(round(size(subdata,1)/2),round(size(subdata,2)/2),:))/mean_noise; %approximate signal to noise ratio
        
        ref_subdata = subdata(:,:,ref_step);
<<<<<<< HEAD
        shifted_data = zeros([size(xsub),numel(frames)]);
        
        fake_dx = step_amplitude*(randn(numel(frames),1));
        fake_dy = step_amplitude*(randn(numel(frames),1));
        
        for j = 1:numel(frames)
            shiftedx = xsub + fake_dx(j);
            shiftedy = ysub + fake_dy(j);
            shifted_data(:,:,j) = interp2(xsub,ysub,ref_subdata,shiftedx,shiftedy);
=======
        shifted_data = zeros([size(xsub),size(data,3)]);        

        for frame = 1:size(data,3)           
            shiftedx = xsub + fake_dx(frame);
            shiftedy = ysub + fake_dy(frame);            
            shifted_data(:,:,frame) = interp2(xsub,ysub,ref_subdata,shiftedx,shiftedy);
>>>>>>> 4107f2ab22e3c274933e723d7e9e1c5b30102c15
        end
        
        %         temp_noise = 2*(rand([size(shifted_data)])-0.5)*std_noise/2 + mean_noise;
        shifted_data(isnan(shifted_data(:))) = subdata(isnan(shifted_data(:)));
        
        %         mean_shifted = mean(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % mean in shifted noise
        %         std_shifted = std(shifted_data(shifted_data(:)<(max(shifted_data(:))/threshfact))); % std in shifted noise
        %         normdata = (shifted_data - mean_shifted)/std_shifted;
        
        %         added_noise = 2*(rand([size(normdata)])-0.5)/(0.5*SNR);
        
        %         noisydata = normdata + added_noise;
        
        %         normnoisydata = (noisydata - mean(noisydata(noisydata(:)<(max(noisydata(:))/threshfact))))/std(noisydata(noisydata(:)<(max(noisydata(:))/threshfact)));
        %         scalednoisydata = normnoisydata*std_noise + mean_noise;
        
        %         calib_params(ptcle,:) = [pkcnt_4QM_calibrator_0(scalednoisydata,tracks,fake_dx,fake_dy,feat_size,delta_fit) ptcle];
        calib_params(ptcle,:) = [FQM(shifted_data,fake_dx,fake_dy,1,[],1) ptcle];
        
        
        %         for frame = 1:ntests
        %
        %             imagesc([subdata(:,:,frame) shifted_data(:,:,frame) scalednoisydata(:,:,frame)])
        %             getframe
        %
        %         end
        
    end
end

res = calib_params;

end