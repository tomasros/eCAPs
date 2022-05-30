function [ranked,score] = fscmrmr(X,y,varargin)
%FSCMRMR Importance of features (predictors) for classification using MRMR algorithm.
%   IDX=FSCMRMR(TBL,Y) ranks predictors for table TBL and class label Y.
%   The ranks are learned using the Minimum Redundancy Maximum Relevance
%   algorithm.
%
%   TBL contains the predictor variables. Y can be any of the following:
%      1. An array of class labels. Y can be a categorical array, logical
%         vector, numeric vector, string array or cell array of character 
%         vectors.
%      2. The name of a variable in TBL. This variable is used as the
%         response Y, and the remaining variables in TBL are used as
%         predictors.
%      3. A formula character vector such as 'y ~ x1 + x2 + x3' specifying
%         that the variable y is to be used as the response, and the other
%         variables in the formula are predictors. Any table variables not
%         listed in the formula are not used.
%
%   IDX=FSCMRMR(X,Y) is an alternative syntax that accepts X as an
%   N-by-P matrix of predictors with one row per observation and one column
%   per predictor. Y is the response and is an array of N class labels. 
%
%   IDX is a 1-by-P vector for P predictors. IDX are indices of
%   columns in X ordered by importance, meaning IDX(1) is the index of
%   the most important predictor.  
%
%   [IDX,SCORES]=FSCMRMR(...) also returns predictor scores SCORES, a
%   1-by-P array for P predictors. SCORES have the same order as predictors
%   in the input data, meaning SCORES(1) is the score for the first
%   predictor. Large score indicates important predictor.
%
%   [...]=FSCMRMR(X,Y,'PARAM1',val1,'PARAM2',val2,...) specifies optional
%   parameter name/value pairs:
%       'CategoricalPredictors' - List of categorical predictors. Pass
%                        'CategoricalPredictors' as one of:
%                          * A numeric vector with indices between 1 and P,
%                            where P is the number of columns of X or
%                            variables in TBL.
%                          * A logical vector of length P, where a true
%                            entry means that the corresponding column of X
%                            or TBL is a categorical variable. 
%                          * 'all', meaning all predictors are categorical.
%                          * A string array or cell array of character
%                            vectors, where each element in the array is
%                            the name of a predictor variable. The names
%                            must match variable names in table TBL.
%                        Default: for a matrix input X, no categorical
%                        predictors; for a table TBL, predictors are
%                        treated as categorical if they are strings, cell
%                        arrays of character vectors, logical, or unordered
%                        of type 'categorical'.
%       'ClassNames'   - Array of class names. Use the data type that
%                        exists in Y. You can use this argument to select a
%                        subset of classes out of all classes in Y.
%                        Default: All class names in Y.
%       'Prior'        - Prior probabilities for each class. Specify as one of:
%                         * A character vector:
%                           - 'empirical' determines class probabilities
%                             from class frequencies in Y
%                           - 'uniform' sets all class probabilities equal
%                         * A vector (one scalar value for each class)
%                         * A structure S with two fields: S.ClassProbs
%                           containing a vector of class probabilities, and
%                           S.ClassNames classes containing the class names and
%                           defining the ordering of classes used for the
%                           elements of this vector.
%                        If you pass numeric values, FSCMRMR normalizes
%                        them to add up to one. Default: 'empirical'
%       'UseMissing'   - Logical value specifying whether missing values 
%                        in the predictors must be used or discarded. When
%                        'UseMissing' is set to true, FSCMRMR groups
%                        together missing values and uses them to compute
%                        mutual information between pairs of features and
%                        between each feature and the response variable.
%                        When 'UseMissing' is set to false, FSCMRMR omits
%                        missing values in each feature while computing the
%                        mutual information. Default: false
%       'Weights'      - Vector of observation weights, one weight per
%                        observation. FSCMRMR normalizes the weights to
%                        add up to the value of the prior probability in
%                        the respective class. Default: ones(size(X,1),1).
%                        For an input table TBL, the 'Weights' value can be
%                        the name of a variable in TBL.
%       'Verbose'      - Verbosity flag, a nonnegative integer:
%                         * 0  - FSCMRMR does not display any diagnostic
%                                messages (default).
%                         * 1  - FSCMRMR displays major diagnostic info.
%                         * >1 - FSCMRMR displays a lot of diagnostic info.
%
%   Example:
%       % Identify important predictors in the ionosphere dataset:
%       load ionosphere
%       [idx,scores] = fscmrmr(X,Y);
%       bar(scores(idx))
%       xlabel('Predictor rank')
%       ylabel('Predictor importance score')

%   Copyright 2019 The MathWorks, Inc.

if nargin > 1
    y = convertStringsToChars(y);
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

args = {'verbose' 'usemissing'};
defs = {        0        false};
[verbose,useNaNs,~,varargin] = internal.stats.parseArgs(args,defs,varargin{:});

validateattributes(verbose,{'numeric'},{'scalar' 'integer' 'nonnegative'},...
    'fscmrmr','Verbose');
verbose = double(verbose);

useNaNs = internal.stats.parseOnOff(useNaNs,'UseMissing');

[X,y,w,dataSummary] = ...
  classreg.learning.classif.FullClassificationModel.prepareData(...
  X,y,varargin{:},'OrdinalIsCategorical',false);

% Accept only dense floating-point.
internal.stats.checkSupportedNumeric('X',X);

D = size(X,2);

y = int32(grp2idx(y));
if max(y) > intmax('int32')-1
    error(message('stats:fscmrmr:TooManyClasses'))
end

% Edge case: one class only
if max(y)==1
    ranked = 1:D;
    score = zeros(1,D,'like',w);
    return
end

% Get indices of categorical features
catpreds = dataSummary.CategoricalPredictors;
iscat = false(1,D);
iscat(catpreds) = true;

% Use at most 256 bins to bin continuous features. These will be lumped
% into larger bins by adaptive binning.
nbin = 256*ones(D,1);

tstart = tic;

% binPredictors returns -1 for missing values; these get transformed into
% zeros.
if ~any(iscat)
    Xbinned = 1 + classreg.learning.treeutils.binPredictors(X,nbin);
else
    Xbinned = zeros(size(X),'int32');
    Xbinned(:,~iscat) = 1 + classreg.learning.treeutils.binPredictors(...
        X(:,~iscat),nbin(~iscat));
    
    Xbinned(:,iscat) = classreg.learning.fsutils.indexCategoricals(X(:,iscat));
        
    if any(max(Xbinned(:,iscat),[],1) > intmax('int32')-1)
        error(message('stats:fscmrmr:TooManyCategories'))
    end
end

% Set parameters for adaptive binning
chi2level = 0.03;
minNodeSize = 5;
miVerbosity = max(0,verbose-1);

% Find MI between features and class label
if verbose > 1
    fprintf('\nComputing mutual info between predictors and class label...\n')
end
Iy = classreg.learning.fsutils.mi2adaptive(...
    2,Xbinned,y,w,iscat,true,useNaNs,chi2level,minNodeSize,miVerbosity);

% Good predictors
predsWithPositiveMI = (Iy > 0);

% Edge case: MI between all features and class label is too small
if ~any(predsWithPositiveMI)
    ranked = 1:D;
    score = zeros(1,D,'like',w);
    return
end

% Prepare output
ranked = zeros(1,D);
score = zeros(1,D,'like',w);

% Place bad features (MI with class label is non-positive) at the end
Dbad = sum(~predsWithPositiveMI);
if Dbad > 0
    % No need to do anything with scores - they are already zeros.
    ranked(end-Dbad+1:end) = find(~predsWithPositiveMI);
    Xbinned = Xbinned(:,predsWithPositiveMI);
    Iy = Iy(predsWithPositiveMI);
    iscat = iscat(predsWithPositiveMI);
end

% Pairwise MI across features
if verbose > 1
    fprintf('\nComputing pairwise mutual info between predictors...\n')
end
Ix = classreg.learning.fsutils.mi2adaptive(...
    1,Xbinned,Xbinned,w,iscat,iscat,useNaNs,chi2level,minNodeSize,miVerbosity);
if verbose > 1
    fprintf('\n')
end

telapsed = toc(tstart);

if verbose>0
    fprintf('Mutual information computed in %10.5e s.\n',telapsed)
end

tstart = tic;

% Convert logical indices to integer indices
predsWithPositiveMI = find(predsWithPositiveMI);

% Find ranks and MIQ ratios on the subset of features with positive MI
% values
[rankedPos,scorePos] = classreg.learning.fsutils.mrmrMIQsort(Iy,Ix);

telapsed = toc(tstart);

if verbose>0
    fprintf('Predictors are ranked in       %10.5e s.\n',telapsed)
end

% Index into ranks and scores for all features
Dgood = numel(predsWithPositiveMI);
ranked(1:Dgood)        = predsWithPositiveMI(rankedPos);
score(ranked(1:Dgood)) = scorePos;

end
