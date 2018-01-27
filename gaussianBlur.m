function [g] =gaussianBlur(image,kernelSize)
kernel = zeros(kernelSize,kernelSize);
sigma = kernelSize/5;
mu =(kernelSize-1)/2;

for i = (-mu):mu
    for j = (-mu):mu
        kernel(i+mu+1,j+mu+1) = 1/(2*pi*sigma^2)*exp(-0.5*(i^2/sigma^2+j^2/sigma^2));
    end
end
g = imfilter(image,kernel);
end
