% Test script to perform k-means clusting
%
% Author: Colorado Reed
% load the data

colors = 'krgb';
shapes='.xo^';

npts = 5000;
mu1 = [1 -1]; Sigma1 = [.9 .4; .4 .3];
d1 = mvnrnd(mu1, Sigma1, npts);
mu2 = [0 1]; Sigma2 = [.9 .4; .4 .9];
d2 = mvnrnd(mu2, Sigma2, npts);
plot(d1(:,1),d1(:,2),'.');
data = [d1; d2];
truth = [ones(npts,1); ones(npts,1) + 1];

%% run homemade kmeans
nclasses = 2;
max_it = 10000;
tic;
[clusts,centrs] = kmeansval(data, nclasses, max_it, false);
toc
acc = sum(truth == clusts)/length(clusts);
if acc < 0.5
   acc = 1 - acc; % flip class assignments 
end
fprintf('accuracy homemade k-means: %0.5f \n\n',acc);

%% compare with Matlab's kmeans
tic;
[clusts,centrs] = kmeans(data, nclasses);
toc
acc = sum(truth == clusts)/length(clusts);
if acc < 0.5
   acc = 1 - acc; % flip class assignments 
end
fprintf('matlab k-means: %0.5f \n',acc);
