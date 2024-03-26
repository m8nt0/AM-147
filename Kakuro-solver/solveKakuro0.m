
% 
% function: solveKakuro0.m
% purpose: to solve a Kakuro problem using mixed integer linear programming
%          (MILP)
% usage:
%   solMat = solveKakuro0(pMat), solve the Kakuro problem as given in the
%       problem matrix 'pMat'. The solution is returned similar to the problem
%       matrix 'pMat' except that all the solutions are filled in their proper
%       places. The default allowable values are between 1 and 9
%       inclusively, by default.
%
%   solMat = solveKakuro0(pMat, maxVal), specify the maximum allowed
%       solution value, i.e. the allowable value is between 1 and 'maxVal'
%       inclusively.
%
% author: LCC Corps
% date: March 8, 2011
% Dependencies:
%   This function requires third party Matlab package Yalmip (available at
%   http://users.isy.liu.se/johanl/yalmip/)
%
% Note:
% 1. The representation used in 'pMat' is quite different from the
%   real-life representation of a Kakuro problem. For details on the
%   representation format, please refer to the 'solveKakuro_demo.m' file.
% 2. This function is very weak especially when there are multiple
%   possible solutions to the given Kakuro problem. If this happens, this
%   function will fail. For better solver, use 'solveKakuro.m' function,
%   which uses binary model.
%
% Kakuro problem representation format:
%   -1 -> to represent an unused/blocked cell
%    0 -> to represent an empty cell which is to be filled for solution
%   x + 1i*y -> to represent the cell that specifies the hints/sums.'x',
%       which is the real part, represent the hint/sum for the vertical
%       (up-down) run, 'y' which is the imaginary part, represent the the
%       hint/sum for the horizontal (left-right) run. For example, please
%       refer to solveKakuro_demo.m
%

function solMat = solveKakuro0(pMat, maxVal)

nargchk(1, 2, nargin);

if (nargin < 2),
    maxVal = 9;
end;


% Formatting problem matrix into appropriate form
% ===============================================
% add a small imaginary value to the cell containing only the sum for the
% vertical column
pMat(imag(pMat) == 0 & pMat > 0) = pMat(imag(pMat) == 0 & pMat > 0) + 1i*eps;  

% fill the unfilled (== 0) elements
ind = find(pMat == 0);
numVars = length(ind);      %length(ind) == number of variables
pMat(ind) = 1:numVars;


% Convert the problem into MILP constraint
% ========================================
% maximum number of constraint estimated. May not be the exact number of
% constraint, but can be trimmed later
numCons = 2*sum(sum(imag(pMat) ~= 0));     

% constraints in matrix form: A*x == y
A = zeros(numCons, numVars);
y = zeros(numCons, 1);
x = intvar(numVars, 1);

[row, col] = size(pMat);

consCount = 0;  % constraint count

F = set(1 <= x <= maxVal);

[rI, cI] = find(imag(pMat) ~= 0);   % find the elements where it is complex, i.e. the given sum 

for i = 1:length(rI),   % loop through each given sum
    vSum = floor(real(pMat(rI(i), cI(i))));
    hSum = floor(imag(pMat(rI(i), cI(i))));
            
    if (vSum > 0),  % vertical sum constraint
        consCount = consCount + 1;
                
        stopInd = find(imag(pMat((rI(i)+1):end, cI(i))) ~= 0 | pMat((rI(i)+1):end, cI(i)) == -1, 1);
                
        if (isempty(stopInd)),  % if end of matrix is encountered
            stopInd = row;
        else    % adjust the end index
            stopInd = stopInd + rI(i) - 1; % stop index for the 'unfilled' cells containing the variable indices
        end;
                
        startInd = rI(i) + 1;      % start index for the 'unfilled' cells containing the variable indices
                
        A(consCount, pMat(startInd:stopInd, cI(i))) = 1;
        y(consCount) = vSum;
                
        % constraint such that no repeating number is a particular
        % sum. Number across different sums can be repeating.
        F = F + set(alldifferent(x(pMat(startInd:stopInd, cI(i)))));
    end;
            
    if (hSum > 0),  % horizontal sum constraint
        consCount = consCount + 1;
                
        stopInd = find(imag(pMat(rI(i), (cI(i)+1):end)) ~= 0 | pMat(rI(i), (cI(i)+1):end) == -1, 1);
                
        if (isempty(stopInd)),  % if end of matrix is encountered
            stopInd = col;
        else    % adjust the end index
            stopInd = stopInd + cI(i) - 1; % stop index for the 'unfilled' cells containing the variable indices
        end;
                
        startInd = cI(i) + 1;      % start index for the 'unfilled' cells containing the variable indices
                
        A(consCount, pMat(rI(i), startInd:stopInd)) = 1;
        y(consCount) = hSum;
                
        % constraint such that no repeating number is a particular
        % sum. Number across different sums can be repeating.
        F = F + set(alldifferent(x(pMat(rI(i), startInd:stopInd))));
    end;
end;

% trimming down the constraint matrix and vector
A((consCount+1):end, :) = [];
y((consCount+1):end) = [];

% Final constraints on the sum
F = F + set(A*x == y);

% Solving the Kakuro problem 
% ==========================
% using yalmip interface and branch and bound algorithm
solvesdp(F, sum(x), sdpsettings('linprog.LargeScale', 'on'));

xsol = double (x);

% Put the solution into the problem matrix and return it as the answer
solMat = pMat;
ind = find (solMat > 0 & imag(solMat) == 0);    % indices for the elements where solution need to be filled into
solMat(ind) = xsol(solMat(ind));        % solMat[r, c] = xsol[solMat[r,c]]




