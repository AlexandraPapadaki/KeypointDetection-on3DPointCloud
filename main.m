clc;
close all;
clear all;

% Applied for both chinese emperor and main building

% Load and organize data
measurements = dlmread('Orangerie.ptx', '', 0, 0);
% Find number of rows and columns
cols = measurements(1,1);
rows = measurements(2,1);
% Create a matix with only point measurements
points = measurements(11:size(measurements,1), :);

% Case of having RGB values and intensity 
if (size(points,2)==7)
    image = reshape(points(:, 5:7), rows, cols, 3);
    % Flip vertically 
    image1 = image(end:-1:1,:,:); 
    % intensity with range 0-1
    image2 = image1/255;
    imshow(image2),title('Image of the scanned scene');
% Case of not having RGB values but only intensity
elseif (size(points,2)==4)
    image = reshape(points(:, 4), rows, cols);
    image2 = (image')/255;
    imshow(image2),title('Scanned scene');    
end

% Histogram stretching
image_gray = reshape(points(:, 4), rows, cols);
% Find min and max and stretch histogram
mini = min(min(image_gray));
maxi = max(max(image_gray));

image_stretched = zeros(rows, cols);
for i = 1:rows
    for j = 1:cols
        image_stretched(i,j) = (image_gray(i,j)-mini)/(maxi-mini);
    end
end
figure;
imshow(image),title('Unstretched image');
figure;
imshow(image_stretched'),title('Image with stretched Histogram');
% Observe difference
diff = image_stretched - image_gray;
figure;
imshow(diff'),title('Difference between Image and stretched image');

% Keypoint generation
% Calculate the gradients of the image
kernelx=part_der(3);
kernely=kernelx';

gx=imfilter(image_gray,kernelx);
gy=imfilter(image_gray,kernely);

gx2 = gx.^2;
gy2 = gy.^2;
gxy = gx.*gy;
gxy = abs(gxy);

% Build matrix N
N = [gx2, gxy; gxy, gy2];

% Apply Gaussian blurring to improve the quality of the results
gx2=gaussianBlur(gx2,5);
gy2=gaussianBlur(gy2,5);
gxy=gaussianBlur(gxy,5);

% Convert images to uint8 
N = im2uint8(N);
gx2 = im2uint8(gx2);
gy2 = im2uint8(gy2);
gxy = im2uint8(gxy);

% Calculate determinant and weight
det = gx2.*gy2-gxy.^2;
trace = gx2+gy2;
%Convert to double
det = im2double(det);
trace = im2double(trace);
trace = trace+0.0001;

% Calculate size and roundness of ellipse and set thresholds
w = det./trace;
q = (4.*det)./(trace.^2);
w_thres = 1.5*(sum(sum(w))/numel(w));
q_thres = 220;

% Detect the keypoints
rows = size(image_gray,1);
cols = size(image_gray,2);
keypoints = [];
for i = 1:rows
    for j = 1:cols
        if w(i,j) > w_thres && q(i,j) > q_thres
            keypoints = [keypoints; i j;];
        end
    end
end

%Show resulting keypoints in the image
figure;
imshow(image/255);
hold on;
scatter(keypoints(:,2),keypoints(:,1));

% Export keypointsinto txt
r=size(keypoints,1);
keys=zeros(r,3);

for l=1:r
    keys(l,:)=measurements((keypoints(l,1)-1)*cols+keypoints(l,2),1:3);
end
% 
% dlmwrite('keypoints_embas.txt',keys,'delimiter',' ');
% imwrite(image/255, 'embas.png');
% imwrite(image/255, 'embas_stretched.png');
% 
% dlmwrite('keypoints_TU-Main-Building.txt',keys,'delimiter',' ');
% imwrite(image, 'tu.png');
% imwrite(image, 'tu_stretched.png');

dlmwrite('keypoints_oran.txt',keys,'delimiter',' ');
imwrite(image, 'oran.png');
imwrite(image, 'oran_stretched.png');



