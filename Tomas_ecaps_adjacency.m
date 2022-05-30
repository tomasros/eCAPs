%%similarity function based on spatial correlation between EEG maps

function D3 = Tomas_ecaps_adjacency(XI,XJ)  

% n = size(XI,2);
D3 = (corr(XI',XJ'));

end

