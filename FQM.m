% Tommy's badass four-quadrant method. If you want to use this bad boy to
% calibrate, set calibrate = 1. If not, then set calibrate = 0.

function res = FQM(subdata,fake_dx,fake_dy,calibrate,p_coef)

%%%%%%%%
% step through each particle and frame then calibrate or measure
%%%%%%%%

if calibrate == 1

    for frame = 1:size(subdata,3)

        subdata_frame = subdata(:,:,frame);
        cutoff = size(subdata_frame,1)/2;

        QLR = subdata_frame(cutoff+1:end,cutoff+1:end);
        QUR = subdata_frame(1:cutoff,cutoff+1:end);
        QLL = subdata_frame(cutoff+1:end,1:cutoff);
        QUL = subdata_frame(1:cutoff,1:cutoff);

        A = sum(QUL(:));
        B = sum(QUR(:));
        C = sum(QLL(:));
        D = sum(QLR(:));

        cnt(frame,:) = [(A+C-B-D)/(A+B+C+D) (A+B-C-D)/(A+B+C+D)];
        refshft(frame,:) = [fake_dx(frame) fake_dy(frame)];

    end

    errorfun = @(p1)squeeze(mean((p1(1)*(cnt(:,1)+p1(2))-refshft(:,1)).^2,1));
    [p1,fval] = fminsearch(errorfun,[range(refshft(:,1))/range(cnt(:,1)),mean(refshft(:,1))]);
    errx = sqrt(fval);

    errorfun = @(p2)squeeze(mean((p2(1)*(cnt(:,2)+p2(2))-refshft(:,2)).^2,1));
    [p2,fval] = fminsearch(errorfun,[range(refshft(:,2))/range(cnt(:,2)),mean(refshft(:,2))]);
    erry = sqrt(fval);

    % for the future: automate error threshold
    if (errx<=1e-1).*(erry<=1e-1)==1

        scatter(p1(1)*(cnt(:,1)+p1(2)),refshft(:,1),'b')
        hold on
        scatter(p2(1)*(cnt(:,2)+p2(2)),refshft(:,2),'g')
        getframe;

        %csvwrite([num2str(round(rand*1000)) '.csv'],[p1(1)*(cnt(:,1)+p1(2)),refshft(:,1) p2(1)*(cnt(:,2)+p2(2)),refshft(:,2)]);

        res = [p1(1) p1(2) errx p2(1) p2(2) erry];

    else

        res = [NaN NaN NaN NaN NaN NaN];

    end
    
else
    
    cnt = zeros(1,3);
        
    for frame = 1:size(subdata,3)

    subdata_frame = subdata(:,:,frame);
    cutoff = size(subdata_frame,1)/2;

    QLR = subdata_frame(cutoff+1:end,cutoff+1:end);
    QUR = subdata_frame(1:cutoff,cutoff+1:end);
    QLL = subdata_frame(cutoff+1:end,1:cutoff);
    QUL = subdata_frame(cutoff+1:end,cutoff+1);

    A = sum(QUL(:));
    B = sum(QUR(:));
    C = sum(QLL(:));
    D = sum(QLR(:));

    cnt = [cnt; [p_coef(1)*(A+C-B-D)/(A+B+C+D)+p_coef(2) p_coef(4)*(A+B-C-D)/(A+B+C+D)+p_coef(5) frame]];

    end
    
    cnt(1,:) = [];
    res = cnt;
    
end

end

