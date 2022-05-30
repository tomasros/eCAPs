%%distance function based on spatial correlation between EEG maps

function abscorr_adjacency = pdist2_adjacency_ecaps(XI,XJ)  

% n = size(XI,2);
abscorr_adjacency = (corr(XI',XJ'));
end

