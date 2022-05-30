function [ labels, nbMoved ] = kaverages ( simils, nbClasses, initClasses, nbIterations )
  %[ labels, nbMoved ] = kaverages ( simils, nbClasses, initClasses, nbIterations )
  %
  %Performs a clustering of the input data using the k-averages algorithm.
  %
  %PARAMETERS:
  %  - 'simils' is a square symmetrical similarity matrix,
  %  - 'nbClasses' is the desired number of output clusters,
  %  - 'initClasses' is an optional parameter specifying the initial allocation
  %     of objects to clusters,
  %  - 'nbIterations' is an optional argument limiting the maximum number of
  %     iterations before convergence (if unspecified, that value is set to 500;
  %     note that convergence is not usually an issue).
  %
  %RETURN VALUES:
  %  - 'labels' gives the result of the clustering as labels in the interval
  %     [1,nbClasses]
  %  - 'nbMoved' gives the number of objects moved from one cluster to another
  %     at each iteration of the algorithm.
  
  % This program was written by Mathias Rossignol & Mathieu Lagrange and
  % is Copyright (C) 2015 IRCAM <http://www.ircam.fr>
  % 
  % This program is free software: you can redistribute it and/or modify it
  % under the terms of the GNU General Public License as published by the Free
  % Software Foundation, either version 3 of the License, or (at your option)
  % any later version.
  % 
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  % or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  % for more details.
  % 
  % You should have received a copy of the GNU General Public License along
  % with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  if nargin<2
    error('Function needs at least 2 parameters');
  end
  
  nbObjects = length(simils);
  
  if nargin > 2
    labels = initClasses;
  else
    labels = 1+floor(nbClasses*rand(1,nbObjects));
  end
  
  selfSimilarities = arrayfun(@(i)simils(i,i), 1:nbObjects);
  
  classSizes = arrayfun(@(c) sum(labels==c), 1:nbClasses);
  
  averageSimilarityAccumulators = zeros(nbClasses, nbObjects);
  classObjectSimilarities = zeros(nbClasses, nbObjects);
  classQualities = zeros(nbClasses,1);
  for ic = 1:nbClasses
    averageSimilarityAccumulators(ic,:) = sum(simils(labels==ic,:));
    classObjectSimilarities(ic,labels==ic) = (averageSimilarityAccumulators(ic,labels==ic) - selfSimilarities(labels==ic)) / (classSizes(ic)-1);
    classObjectSimilarities(ic,labels~=ic) = averageSimilarityAccumulators(ic,labels~=ic) / classSizes(ic);
    classQualities(ic) = mean (classObjectSimilarities(ic,labels==ic));
  end
  
  % Criterion function defined below
  critFun = @objectNormalizedCriterion;
  
  if ~exist('nbIterations', 'var'), nbIterations = 500; end
  
  minIncrement = 0;
  
  nbMoved = [];

  changes = 1;
  loop = 0;
  while changes && loop < nbIterations
    loop = loop+1;
    changes = progressiveStep();
    nbMoved = [nbMoved changes]; %#ok<AGROW>
  end
  
  if (loop == nbIterations)
    disp('Warning: progressive k-averages clustering failed to converge within maximum number of iterations.');
  end
  
  % Defining the optimisation criterion as nested function here
  % Parameters: object index o, new class index c
  % Return: impact on the global objective function of moving object o to class c
  
  function criterion = objectNormalizedCriterion (c, o)
    c1 = labels(o);
    n1 = classSizes(c1);
    c2 = c;
    n2 = classSizes(c2);
    if (c1==c2)
      criterion = 0;
    else
      impactOfLeavingOld =  2 *(n1-1)*(classQualities(c1)-classObjectSimilarities(c1,o)) / (n1-2) - classQualities(c1);
      newClassQualityC2 = ((n2-1)*classQualities(c2)+2*classObjectSimilarities(c2,o)) / (n2+1);
      impactOfJoiningNew = -2 * n2   *(classQualities(c2)-classObjectSimilarities(c2,o)) / (n2+1) + newClassQualityC2;
      criterion = impactOfLeavingOld + impactOfJoiningNew;
    end
  end
  
  % Defining the base strategy as nested function here
  % One function call = one iteration step
  % Returns number of changes
  
  function changes = progressiveStep()
    changes = 0;
    % Taking objects one by one, compute updated criterion then
    % take allocation decision
    for o=1:nbObjects
      criteria = arrayfun(@(c) critFun(c, o), 1:nbClasses);
      [v,newC] = max(criteria);
      if (v>minIncrement)
        % Incremental update of class properties
        changes = changes+1;
        oldC = labels(o);
        labels(o) = newC;
        classQualities(oldC) = (classSizes(oldC)*classQualities(oldC)-2*classObjectSimilarities(oldC,o))/(classSizes(oldC)-2);
        classQualities(newC) = ((classSizes(newC)-1)*classQualities(newC)+2*classObjectSimilarities(newC,o))/(classSizes(newC)+1);
        
        % Incremental update of class-object similarities
        averageSimilarityAccumulators(oldC,:) = averageSimilarityAccumulators(oldC,:) - simils(o,:);
        averageSimilarityAccumulators(newC,:) = averageSimilarityAccumulators(newC,:) + simils(o,:);
        classSizes(oldC) = classSizes(oldC)-1;
        classSizes(newC) = classSizes(newC)+1;
        for c=[oldC,newC]
          classObjectSimilarities(c,labels==c) = (averageSimilarityAccumulators(c,labels==c) - selfSimilarities(labels==c)) / (classSizes(c)-1);
          classObjectSimilarities(c,labels~=c) = averageSimilarityAccumulators(c, labels~=c) / classSizes(c);
        end
      end
    end
  end
end
