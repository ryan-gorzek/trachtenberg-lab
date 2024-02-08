
function sbxeyemotion_cv2(fn)

% memory mapped files

image_out = memmapfile('h:\2pdata\sbx2dlc','Writable',true,'Format',{'uint8', [112 160], 'img'});
data_in   = memmapfile('h:\2pdata\dlc2sbx','Writable',false,'Format',{'single', [3 8], 'pose'});

% assumes eyetracker_cv2.py is running (findblob_test.py)

vr = VideoReader([fn '_eye.mj2']);

Centroid = zeros(vr.NumFrames,2);

for(k=1:vr.NumFrames)
    x = vr.readFrame;
    x(1)=1;
    image_out.Data.img = x';
    while(image_out.Data.img(1,1) ~=0)      % wait for it to finish
    end
    Centroid(k,:) = data_in.Data.pose(1:2,1)';
end

save(fn,'Centroid','-append');     % append the motion estimate data...
