%%distance function based on spatial correlation between EEG maps

function D3 = Tomas_microstatetopo_adjacency(XI,XJ)  

% n = size(XI,2);
D3 = abs(corr(XI',XJ'));

end

