nsamples = 5;
npoints = 50;

for k = 1:nsamples
    currentData = rand(npoints,1);
    sampleMean(k) = mean(currentData);
end
overallMean = mean(sampleMean)
if overallMean < 0.55
    disp('mean less than expected')
else 
    disp('mean more than expected')
end


    