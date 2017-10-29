function y_smooth = smooth(y, box_pts)

box = ones(box_pts,1)./box_pts;
y_smooth = conv2(y, box,'same');


