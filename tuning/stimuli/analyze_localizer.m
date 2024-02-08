function analyze_localizer(animal,unit,expt)
    expname = sprintf('%s_%s_%s',animal,unit,expt);
    
    %Load the event file

    datapath = 'd:\2pdata\'; 
    %thepath = '\\vstim\Users\jt_admin\Documents\dario\object_view';
    thepath = 'C:\\ObjectViewLogs';

    % thepath = 'C:\\Users\\dario\\Desktop\\object_view';
    sprintf('%s\\%s\\%s',datapath,animal,expname);
    info = load(sprintf('%s\\%s\\%s\\%s',datapath,animal,expname,expname));
    info = info.info;
    
    ue = unique(info.event_id);
    if ismember(64,ue) 
        evts = info.event_id == 64;
        disp('Using Frame2TTL');
    else
        evts = info.event_id == 1;
        disp('Using Software TTL');
    end
    
    evts = info.event_id == 1;     % use sync 
    eventstarts = info.frame(evts);
    %Load the log file
    d = load(sprintf('%s\\%s.log_17',thepath,expname));
    ys = d(:,2);
    xs = d(:,3);
    
    sz = [max(ys),max(xs)];
    
    sprintf('%s\\%s\\%s_realtime',datapath,animal,expname);
    
    d = load(sprintf('%s\\%s\\%s\\%s_realtime',datapath,animal,expname,expname));
    Y = -d.rtdata;
    Y = bsxfun(@minus,Y,median(Y));
    Y = bsxfun(@times,Y,1./std(Y));
    
    rgs = bsxfun(@plus,eventstarts,(5:16));
    
    Y = [Y;zeros(100,size(Y,2))];
    
    meanY = squeeze(mean(reshape(Y(rgs(:),:),[size(rgs),size(Y,2)]),2));
    
    figure;
    
%     nx = ceil(sqrt(size(Y,2))*sqrt(1.7));
%     ny = ceil(sqrt(size(Y,2))/sqrt(1.7));
    
    nx = 1;
    ny = size(Y,2);
%     
%     ny = 3;
%     nx = 3;
    
Q = [];

    for ii = 1:size(Y,2)
        S = sparse(ys,xs,meanY(:,ii),sz(1),sz(2)); %quick way to get receptive field
        subplot(ny,nx,ii);
        if ii==1
            Z = full(S);
        else
            Z = Z+full(S);
        end
        imagesc(S,full(max(S(:)))*[-1,1]);
        axis equal off
        
        [xx,yy]= meshgrid(1:size(Z,2),1:size(Z,1));
        f2dg = @(x,xdata) x(1)+x(2)*exp(- ((xdata(:,:,1)-x(3)).^2 + (xdata(:,:,2)-x(4)).^2) / (2*x(5)^2));
        xd = zeros([size(xx) 2]);
        xd(:,:,1) = xx;
        xd(:,:,2) = yy;
        M = max(full(S(:)));
        [ii,jj] = find(M==full(S));
        xopt = lsqcurvefit(f2dg,[0 M jj ii 2], xd,full(S));
        hold on,plot(xopt(3),xopt(4),'ro','markersize',12,'linewidth',2);
        Q = [Q ; xopt(3) xopt(4)];
    end
    
    figure
    imagesc(Z,full(max(Z(:)))*[-1,1]);
    axis equal off
    [xx,yy]= meshgrid(1:size(Z,2),1:size(Z,1));
    f2dg = @(x,xdata) x(1)+x(2)*exp(- ((xdata(:,:,1)-x(3)).^2 + (xdata(:,:,2)-x(4)).^2) / (2*x(5)^2));
    xd = zeros([size(xx) 2]);
    xd(:,:,1) = xx;
    xd(:,:,2) = yy;
    M = max(full(Z(:)));
    [ii,jj] = find(M==full(Z));
    xopt = lsqcurvefit(f2dg,[0 M jj ii 2], xd,full(Z));
    hold on,plot(xopt(3),xopt(4),'ro','markersize',12,'linewidth',2);
    
    display_param;
    dx = displayWidth/size(Z,2);
    dy = displayHeight/size(Z,1);
    
    xcoord = round((xopt(3)-1)*dx + dx/2);
    ycoord = round((xopt(4)-1)*dy + dy/2);
    ycoord = displayHeight-ycoord;
    
    title(sprintf('Pop Ctr at [%d,%d]',xcoord,ycoord));
    
    for i=1:size(Q,1)
        t = text(Q(i,1),Q(i,2),sprintf('%d',i));
        t.FontSize = 14;
        t.HorizontalAlignment = 'center';
    end
end


    