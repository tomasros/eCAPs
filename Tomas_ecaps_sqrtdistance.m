%%distance function based on spatial correlation between EEG source space microstates

function D3 = Tomas_ecaps_sqrtdistance(XI,XJ)  

% n = size(XI,2);
% D3 = sqrt(1-abs(corr(XI',XJ', 'Type','Pearson')));
D3 = sqrt(1-(corr(XI',XJ', 'Type','Pearson')));
% D3 = (1-abs(corr(XI',XJ', 'Type','Pearson'))); %% can't do a square root when the correlation is negative
end

