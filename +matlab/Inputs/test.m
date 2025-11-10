clear
clc
load("CorrelationData.mat")

[idx, Z, order] = cluster_from_corr_hier(corrMatrix, 4, false);

function [idx, Z, order] = cluster_from_corr_hier(R, k, useAbs)
% Cluster variables given a correlation matrix R into k groups
% useAbs = true  -> cluster on |corr| (treat positive/negative as similar)
% useAbs = false -> cluster on signed corr (negatives are dissimilar)

if nargin < 3, useAbs = true; end
R = (R + R.')/2;              % force symmetry
R(1:size(R,1)+1:end) = 1;     % clean diagonal

% Convert to a distance. Two common choices:
%   1 - |R|  (unsigned similarity)
%   1 -  R   (signed similarity)
if useAbs
    S = abs(R);
else
    S = R;
end
D = 1 - S;
D(1:size(D,1)+1:end) = 0;     % zero self-distance

% Hierarchical clustering (average linkage is a good default)
% linkage needs a condensed distance vector:
dvec = squareform(D, 'tovector');
Z = linkage(dvec, 'average');

% Cut the tree into k clusters
idx = cluster(Z, 'maxclust', k);

% Order variables to show block structure (optimal leaf ordering)
Zord = optimalleaforder(Z, dvec);
order = optimalleaforder(Z, dvec); %#ok<NASGU> % alias for clarity

% Quick plots (optional)
figure; dendrogram(Z, 0); title('Hierarchical clustering dendrogram');
figure; imagesc(R(Zord, Zord)); axis image; colorbar
title('Correlation matrix reordered by clusters');
end
