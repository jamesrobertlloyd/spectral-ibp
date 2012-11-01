function [classes,centrs]=kmeansval(indata, k, its, visualize, verbose)
% function [classes,centrs]=kmeansval(indata,k,its,visualize) 
%
% This function finds the assignments of `indata` using the k-means algorithm
%
% INPUT
% k is the number of clusters to use
% its is the max number of iterations (break if convergence)[default 1000]
% visualize indicates whether or not to visualize the clustering (only use with 2d data) [default false]
%
% Author: Colorado Reed 
if nargin < 3
    its = 1000;
end
if nargin < 4
    visualize = false;
end
if nargin < 5
    verbose = false;
end

[nrows, ncols] = size(indata);
 
% assign random points for initialization
temp=randperm(nrows);
centrs = indata(temp(1:k),:);
dists = zeros(nrows,k) + NaN;
prevclass = zeros(nrows,1) + NaN;

for ii=1:its
    % find distances 
    for j=1:k
        usecentr = centrs(j,:);
        % TODO use fileexchange distances
        dists(:,j) = sum((indata - usecentr(ones(nrows,1),:)).^2,2); % improve this via file exchange
    end
    
    % find the nearest centroid
    [~,classes] = min(dists,[],2);
    
    % break if classes converge
    if prevclass==classes
        fprintf('k-means converged in %i iterations\n',ii);
        break
    end
    
    prevclass = classes;
    
    if visualize
        colors = 'krgby';
        shapes='.xo^+';
        for j=1:k
            classdata = classes==j;
            plot(indata(classdata,1),indata(classdata,2),[shapes(mod(j,length(shapes))+1) colors(mod(j,length(colors))+1)] ,'markersize',10);
            hold on;
            plot(centrs(j,1),centrs(j,2),['p' colors(mod(j,length(colors))+1)],'markersize',20','markerfacecolor','m');
        end
        pause(0.1);
        hold off
        grid off;
        drawnow;
    end
    
     % determine new centroid
    for j=1:k
        centrs(j,:) = mean(indata(classes==j,:),1);
    end
    
    if ii == its
        frintf('WARNING: Did not converge after %i iterations\n',its);
    end
end

end

