
% this function is to demonstrate the use of solveKakuro0.m and solveKakuro.m

clear all;

% Kakuro problem #1
% Real-life representation of the problem
% 
% +-----+-----+-----+-----+-----+-----+-----+
% |xxxxx|xxxxx|15\  |xxxxx|xxxxx|xxxxx|xxxxx|
% |xxxxx|xxxxx|     |24\  |xxxxx|15\  | 6\  |
% |xxxxx|24\16|     |     |21\11|     |     |
% |xxxxx|     |17\22|     |     |     |     |
% |  \30|     |     |     |     |17\  |     |
% |  \16|     |     |  \17|     |     |xxxxx|
% |xxxxx|xxxxx|xxxxx|xxxxx|xxxxx|     |xxxxx|
% +-----+-----+-----+-----+-----+-----+-----+
%
% where,
%   xxxxx -> represent unused block/cell
%   \ -> (forward slash) represent the divider between the vertical and
%       horizontal hints/sums in a cell

% p1 = [-1   -1      15        -1    -1       -1 -1; ...
%       -1   -1       0        24    -1       15  6; ...
%       -1   24+1i*16  0         0     21+1i*11  0  0; ...
%       -1    0       17+1i*22   0     0        0  0; ...
%       1i*30  0       0         0     0       17  0; ...
%       1i*16  0       0        1i*17   0        0 -1; ...
%       -1   -1      -1        -1    -1        0 -1];
%  
% s1 = solveKakuro(p1);


% Kakuro problem #2
% Real-life representation of the problem
% 
% +-----+-----+-----+-----+-----+-----+-----+
% |xxxxx|xxxxx|xxxxx|22\  |36\  | 4\  |xxxxx|
% |xxxxx|17\  |14\10|     |     |     |     |
% |  \31|     |     |     |     |     |xxxxx|
% |  \26|     |     |     |     | 8\  |11\  |
% |xxxxx|xxxxx|15\10|     |     |     |     |
% |xxxxx|  \34|     |     |     |     |     |
% |  \29|     |     |     |     |xxxxx|xxxxx|
% +-----+-----+-----+-----+-----+-----+-----+
%
% where,
%   xxxxx -> represent unused block/cell
%   \ -> (forward slash) represent the divider between the vertical and
%       horizontal hints/sums in a cell

% p2 = [-1   -1   -1       22    36   4 -1; ...
%       -1   17   14+1i*10   0     0   0  0; ...
%       1i*31  0    0        0     0   0 -1; ...
%       1i*26  0    0        0     0   8 11; ...
%       -1   -1   15+1i*10   0     0   0  0; ...
%       -1   1i*34  0        0     0   0  0; ...
%       1i*29  0    0        0     0  -1 -1];
%  
% s2 = solveKakuro(p2);


% Kakuro problem #3
% Real-life representation of the problem
% 
% +-----+-----+-----+-----+-----+-----+-----+-----+
% |xxxxx|23\  |30\  |xxxxx|xxxxx|27\  |12\  |16\  |
% |  \16|     |     |xxxxx|17\24|     |     |     |
% |  \17|     |     |15\29|     |     |     |     |
% |  \35|     |     |     |     |     |12\  |xxxxx|
% |xxxxx|  \7 |     |     | 7\8 |     |     | 7\  |
% |xxxxx|11\  |10\16|     |     |     |     |     |
% |  \21|     |     |     |     |  \5 |     |     |
% |  \6 |     |     |     |xxxxx|  \3 |     |     |
% +-----+-----+-----+-----+-----+-----+-----+-----+
%
% where,
%   xxxxx -> represent unused block/cell
%   \ -> (forward slash) represent the divider between the vertical and
%       horizontal hints/sums in a cell

p3 = [-1     23     30         -1       -1       27   12  16; ...
      1i*16   0      0         -1       17+1i*24  0    0   0; ...
      1i*17   0      0         15+1i*29  0        0    0   0; ...
      1i*35   0      0          0        0        0   12  -1; ...
      -1     1i*7    0          0        7+1i*8   0    0   7; ...
      -1     11     10+1i*16    0        0        0    0   0; ...
      1i*21   0      0          0        0       1i*5  0   0; ...
      1i*6    0      0          0       -1       1i*3  0   0];
 
% F = [];
% s3 = solveKakuro(p3);
[s3, ~] = solveKakuro(p3, 9);