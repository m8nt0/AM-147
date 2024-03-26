% 
% function: solveKakuro.m
% purpose: to solve a Kakuro problem using mixed integer linear programming
%          (MILP) with binary model
% usage:
%   solMat = solveKakuro(pMat), solve the Kakuro problem as given in the
%       problem matrix 'pMat'. The solution is returned similar to the problem
%       matrix 'pMat' except that all the solutions are filled in their proper
%       places. The default allowable values are between 1 and 9
%       inclusively, by default.
%
%   solMat = solveKakuro(pMat, maxVal), specify the maximum allowed
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
% 2. This function is much stronger (as compared to solveKakuro0.m) in a
%   sense that it will still provide valid solution even in the case where
%   there exists multiple possible solutions for the given Kakuro problem.
%   However, the formulation of the problem using binary model is much less
%   intuitive as compared to non-binary model (refer to solveKakuro0.m).
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

function [solMat, xsol] = solveKakuro(pMat, maxVal)

narginchk(1, 2);

if nargin < 2
    maxVal = 9;
end

% Formatting problem matrix into appropriate form
% ===============================================
% add a small imaginary value to the cell containing only the sum for the
% vertical column
pMat(imag(pMat) == 0 & pMat > 0) = pMat(imag(pMat) == 0 & pMat > 0) + 1i*eps;  

% fill the unfilled (== 0) elements
ind = find(pMat == 0);
numVars = length(ind);      % length(ind) == number of variables
pMat(ind) = 1:numVars;

% Convert the problem into MILP constraint
% ========================================
x = binvar(numVars, maxVal, 'full');    % only number 1 - 9 is allowed

[row, ~] = size(pMat);

% constraint such that each variable must consist of a value between 1 - 9
% inclusive
F = sum(x, 2) == 1;    

[rI, cI] = find(imag(pMat) ~= 0);   % find the elements where it is complex, i.e. the given sum 

valVec = (1:maxVal)';

for i = 1:length(rI)   % loop through each given sum
    vSum = floor(real(pMat(rI(i), cI(i))));
    hSum = floor(imag(pMat(rI(i), cI(i))));
            
    if vSum > 0  % vertical sum constraint    
        stopInd = find(imag(pMat((rI(i)+1):end, cI(i))) ~= 0 | pMat((rI(i)+1):end, cI(i)) == -1, 1);
                
        if isempty(stopInd)  % if end of matrix is encountered
            stopInd = row;
        else    % adjust the end index
            stopInd = stopInd + rI(i) - 1; % stop index for the 'unfilled' cells containing the variable indices
        end
                
        startInd = rI(i) + 1;      % start index for the 'unfilled' cells containing the variable indices
                
        if stopInd <= startInd
            warning('MATLAB:solveKakuro', 'Invalid Kakuro problem! Please check.\n');
        else
            % constraint for a given sum (from problem)
            F = [F, sum(x(pMat(startInd:stopInd, cI(i)), :)*valVec) == vSum];

            % constraint such that no repeating value in a particular  sum
            F = [F, 0 <= sum(x(pMat(startInd:stopInd, cI(i)), :), 1) <= 1];
        end
    end
            
    if hSum > 0  % horizontal sum constraint
        % fprintf('hSum: %d, startInd: %d, stopInd: %d\n', hSum, startInd, stopInd); % Debugging information
        stopInd = find(imag(pMat(rI(i), (cI(i)+1):end)) ~= 0 | pMat(rI(i), (cI(i)+1):end) == -1, 1);
        
        if isempty(stopInd)  % if end of matrix is encountered
            [row, col] = size(pMat);
            stopInd = col;
        else    % adjust the end index
            stopInd = stopInd + cI(i) - 1; % stop index for the 'unfilled' cells containing the variable indices
        end
                    
        startInd = cI(i) + 1;      % start index for the 'unfilled' cells containing the variable indices
                    
        if stopInd <= startInd
            warning('MATLAB:solveKakuro', 'Invalid Kakuro problem! Please check.\n');
        else
            % constraint for a given sum (from problem)
            F = [F, sum(x(pMat(rI(i), startInd:stopInd), :)*valVec) == hSum];
    
            % constraint such that no repeating value in a particular  sum
            F = [F, 0 <= sum(x(pMat(rI(i), startInd:stopInd), :), 1) <= 1];
        end
    end

end

% Solving the Kakuro problem 
% ==========================
% using yalmip interface and branch and bound algorithm
info = solvesdp(F, [], sdpsettings('linprog.LargeScale', 'on'));

if info.problem
    warning('MATLAB:solveKakuro - %s\n', info.info);
end

xbin = double(x);
xsol = xbin*valVec;

% Put the solution into the problem matrix and return it as
solMat = pMat;
ind = find(solMat > 0 & imag(solMat) == 0); % indices for the elements where solution need to be filled into
solMat(ind) = xsol(solMat(ind)); % solMat[r, c] = xsol[solMat[r,c]]
