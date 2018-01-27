function g = part_der(kernelSize)
g = zeros (kernelSize,kernelSize);
sigma = 0.7;
mu =(kernelSize-1)/2;

for i = (-mu):mu
    for j = (-mu):mu
        g(i+mu+1,j+mu+1) = -j*exp(-0.5*(i^2/sigma^2+j^2/sigma^2));
    end
end
end

