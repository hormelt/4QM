function res = pkcnt_4QM_calibrator(b,tracks,fake_dx,fake_dy,feat_size,delta_fit)

%     get rid of edge particles
if ~isempty(tracks)
    dx1 = abs(tracks(:,1)-1);
    dy1 = abs(tracks(:,2)-1);
    dx2 = abs(tracks(:,1)-size(b,2));
    dy2 = abs(tracks(:,2)-size(b,1));
    
    cutoff = feat_size+4+feat_size/2;
    keeplist = logical((dx1>cutoff).*(dx2>cutoff).*(dy1>cutoff).*(dy2>cutoff));
    tracks(~keeplist,:)=[];
    
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
    
    for ptcle = 1:max(tracks(:,6))
        
        if sum(tracks(:,6)==ptcle)~=0
            ref_cnts = [mean(tracks(tracks(:,6)==ptcle,1:2),1)];
            
            frames = tracks(tracks(:,6)==ptcle,5);
            
            x_coarse = round(ref_cnts(1));
            y_coarse = round(ref_cnts(2));
            
            %     first do calibration using shifted_data
            
            subdata = b((y_coarse-(feat_size-delta_fit)):(y_coarse+(feat_size-delta_fit)),(x_coarse-(feat_size-delta_fit)):(x_coarse+(feat_size-delta_fit)),frames);
            
            count = 0;
            
            for frame = 1:size(subdata,3)
                
                count = count + 1;
                
                finedata(:,:,frame) = interp2(xsub,ysub,subdata(:,:,frame),xfine,yfine,'cubic');
                
                cutoffx = xfine(1,(size(xfine,1)-1)/2+1);
                cutoffy = yfine((size(xfine,2)-1)/2+1,1);
                
                QLR = finedata(:,:,frame).*(xfine>cutoffx).*(yfine>cutoffy);
                QUR = finedata(:,:,frame).*(xfine>cutoffx).*(yfine<cutoffy);
                QLL = finedata(:,:,frame).*(xfine<cutoffx).*(yfine>cutoffy);
                QUL = finedata(:,:,frame).*(xfine<cutoffx).*(yfine<cutoffy);
                
                A = sum(QUL(:));
                B = sum(QUR(:));
                C = sum(QLL(:));
                D = sum(QLR(:));
                
                cnt(count,:) = [(A+C-B-D)/(A+B+C+D) (A+B-C-D)/(A+B+C+D)];
                refcnt(count,:) = [fake_dx(frames(frame)) fake_dy(frames(frame))];
                
            end
            
            errorfun = @(p1)squeeze(mean((p1(1)*(cnt(:,1)+p1(2))-refcnt(:,1)).^2,1));
            [p1,fval] = fminsearch(errorfun,[range(refcnt(:,1))/range(cnt(:,1)),mean(refcnt(:,1))]);
            errx(ptcle) = fval;
            
            errorfun = @(p2)squeeze(mean((p2(1)*(cnt(:,2)+p2(2))-refcnt(:,2)).^2,1));
            [p2,fval] = fminsearch(errorfun,[range(refcnt(:,2))/range(cnt(:,2)),mean(refcnt(:,2))]);
            erry(ptcle) = fval;
            
            % for the future: automate error threshold
            
            if (errx(ptcle)<=1e-2).*(erry(ptcle)<=1e-2)==1
                
                ptclecnt = ptclecnt + 1;
                
                %             scatter(p1(1)*(cnt(:,1)+p1(2)),refcnt(:,1),'b')
                %             hold on
                %             scatter(p2(1)*(cnt(:,2)+p2(2)),refcnt(:,2),'g')
                %             getframe
                
                calib_params(ptclecnt,:) = [p1(1) p1(2) errx(ptcle) p2(1) p2(2) erry(ptcle) ref_cnts ptcle];
                
            end
            
            
            clear cnt
            clear refcnt
            
        end
        
    end
end

if count~=0
    res = calib_params;
else
    res = [NaN NaN NaN NaN NaN NaN NaN];
    
end

end
