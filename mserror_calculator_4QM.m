function res = mserror_calculator_4QM(data,tracks,param)

% Calculates the error
%
% INPUTS
%
% OUTPUTS
%

%%%%
%unpack params
%%%%%
% mean_noise = param.mean_noise;
% std_noise = param.std_noise;
feat_size = param.feat_size;
delta_fit = param.delta_fit;
step_amplitude = param.step_amplitude;
ntests = param.ntests;
threshfact = param.threshfact;
ref_cnts = param.ref_cnts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup the various coordinate systems and coordinate shifts used in this
%code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fake_dx = step_amplitude*(randn(ntests,1));
fake_dy = step_amplitude*(randn(ntests,1));

% fake_dx = 2*(rand(ntests,1)-0.5)*step_amplitude;
% fake_dy = 2*(rand(ntests,1)-0.5)*step_amplitude;

x = 1:size(data,2);
y = 1:size(data,1);
[x y] = meshgrid(x,y);

xsub = 1:(2*(feat_size-delta_fit));
ysub = 1:(2*(feat_size-delta_fit));
[xsub ysub] = meshgrid(xsub,ysub);

%%%%%%%%%%%%%%%%
%step through particles, find best frame to use, make sure noise in
%calibration data represents noise in real data
%%%%%%%%%%%%%%%%%

ptclecnt = 0;

for ptcle = 1:max(tracks(:,6))
    
    subtracks = tracks(tracks(:,6)==ptcle,1:2);
    
    if ~isempty(subtracks)
        ptclecnt = ptclecnt + 1;
        
        distance_metric = sqrt((subtracks(:,1)-ref_cnts(ptclecnt,1)).^2 + (subtracks(:,2)-ref_cnts(ptclecnt,2)).^2);
        [C ref_step] = min(distance_metric);
        
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
        
        %     first do calibration using shifted_data
        
        subdata = data(rows,cols,frames);
        
        mean_noise = mean(subdata(subdata(:)<(max(subdata(:))/threshfact))); % target mean in noise
        std_noise = std(subdata(subdata(:)<(max(subdata(:))/threshfact))); % target std in noise
        max_noise = max(abs(subdata(subdata(:)<(max(subdata(:))/threshfact))-mean_noise));
        SNR = mean(subdata(round(size(subdata,1)/2),round(size(subdata,2)/2),:))/mean_noise; %approximate signal to noise ratio
        
        ref_subdata = subdata(:,:,ref_step);
        
        for frame = 1:ntests
            
            shiftedx = xsub + fake_dx(frame);
            shiftedy = ysub + fake_dy(frame);
            
            shifted_data(:,:,frame) = interp2(xsub,ysub,ref_subdata,shiftedx,shiftedy);
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