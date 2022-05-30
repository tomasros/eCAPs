%%distance function based on spatial correlation between EEG maps

function abscorr_adjacency = pdist2_adjacency(XI,XJ)  

% n = size(XI,2);
abscorr_adjacency = abs(corr(XI',XJ'));
end

