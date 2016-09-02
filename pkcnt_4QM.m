function res = pkcnt_4QM(b,calibration_params,param);
%%%%
%unpack params
%%%%%
mean_noise = param.mean_noise;
std_noise = param.std_noise;
feat_size = param.feat_size;
delta_fit = param.delta_fit;
step_amplitude = param.step_amplitude;
ntests = param.ntests;
threshfact = param.threshfact;
noise_sz = param.noise_sz;
widthcut =  param.widthcut;

%%%%
%setup coords for rest of measurement
%%%
x = 1:size(b,2);
y = 1:size(b,1);
[x y] = meshgrid(x,y);

xsub = 1:(2*(feat_size-delta_fit)+1);
ysub = 1:(2*(feat_size-delta_fit)+1);
[xsub ysub] = meshgrid(xsub,ysub);

xfine = 1:0.1:(2*(feat_size-delta_fit)+1);
yfine = 1:0.1:(2*(feat_size-delta_fit)+1);
[xfine yfine] = meshgrid(xfine,yfine);

%%%%%%%%
% step through each particle, calibrate, then measure
%%%%%%%%

ptclecnt = 0;
cnt = [0 0 0 0];

for step = 1:size(calibration_params,1)
    
    ptcle = calibration_params(step,end);
    ref_cnts = calibration_params(step,7:8);
    p1(1) = calibration_params(step,1);
    p1(2) = calibration_params(step,2);
    p2(1) = calibration_params(step,4);
    p2(2) = calibration_params(step,5);
    
    frames = 1:size(b,3);
    
    x_coarse = round(ref_cnts(1));
    y_coarse = round(ref_cnts(2));
    
    %     first do calibration using shifted_data
    
    subdata = b((y_coarse-(feat_size-delta_fit)):(y_coarse+(feat_size-delta_fit)),(x_coarse-(feat_size-delta_fit)):(x_coarse+(feat_size-delta_fit)),frames);
       
    for frame = 1:size(subdata,3)
        
        finedata = interp2(xsub,ysub,subdata(:,:,frame),xfine,yfine,'cubic');
        
        cutoffx = xfine(1,(size(xfine,1)-1)/2+1);
        cutoffy = yfine((size(xfine,2)-1)/2+1,1);
        
        QLR = finedata.*(xfine>cutoffx).*(yfine>cutoffy);
        QUR = finedata.*(xfine>cutoffx).*(yfine<cutoffy);
        QLL = finedata.*(xfine<cutoffx).*(yfine>cutoffy);
        QUL = finedata.*(xfine<cutoffx).*(yfine<cutoffy);
        
        A = sum(QUL(:));
        B = sum(QUR(:));
        C = sum(QLL(:));
        D = sum(QLR(:));
        
        cnt = [cnt; [x_coarse + p1(1)*(A+C-B-D)/(A+B+C+D)+p1(2) y_coarse + p2(1)*(A+B-C-D)/(A+B+C+D)+p2(2) frame ptcle]];
        
    end
    
end
cnt(1,:)=[];
res = cnt;

end
