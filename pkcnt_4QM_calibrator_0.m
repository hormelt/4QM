function res = pkcnt_4QM_calibrator_0(data,tracks,fake_dx,fake_dy,feat_size,delta_fit)

%%%%
%setup coords for rest of measurement
%%%
x = 1:size(data,2);
y = 1:size(data,1);
[x y] = meshgrid(x,y);

xfine = 1:0.1:(2*(feat_size-delta_fit)+1);
yfine = 1:0.1:(2*(feat_size-delta_fit)+1);
[xfine yfine] = meshgrid(xfine,yfine);

%%%%%%%%
% step through each particle, calibrate, then measure
%%%%%%%%
for frame = 1:size(data,3)
        
    finedata(:,:,frame) = interp2(x,y,data(:,:,frame),xfine,yfine,'cubic');
    
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
    
    cnt(frame,:) = [(A+C-B-D)/(A+B+C+D) (A+B-C-D)/(A+B+C+D)];
    refcnt(frame,:) = [fake_dx(frame) fake_dy(frame)];
    
end

errorfun = @(p1)squeeze(mean((p1(1)*(cnt(:,1)+p1(2))-refcnt(:,1)).^2,1));
[p1,fval] = fminsearch(errorfun,[range(refcnt(:,1))/range(cnt(:,1)),mean(refcnt(:,1))]);
errx = fval;

errorfun = @(p2)squeeze(mean((p2(1)*(cnt(:,2)+p2(2))-refcnt(:,2)).^2,1));
[p2,fval] = fminsearch(errorfun,[range(refcnt(:,2))/range(cnt(:,2)),mean(refcnt(:,2))]);
erry = fval;

% for the future: automate error threshold
if (errx<=1e-2).*(erry<=1e-2)==1
    
    scatter(p1(1)*(cnt(:,1)+p1(2)),refcnt(:,1),'b')
    hold on
    scatter(p2(1)*(cnt(:,2)+p2(2)),refcnt(:,2),'g')
    getframe
    
    csvwrite([num2str(round(rand*1000)) '.csv'],[p1(1)*(cnt(:,1)+p1(2)),refcnt(:,1) p2(1)*(cnt(:,2)+p2(2)),refcnt(:,2)]);
    
    res = [p1(1) p1(2) errx p2(1) p2(2) erry];
    
else
    
    res = [NaN NaN NaN NaN NaN NaN];
    
end

end
