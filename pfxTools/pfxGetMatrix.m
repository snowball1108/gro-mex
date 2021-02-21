function [matrix, indices] = pfxGetMatrix(pfx, fmdim, varargin)
%PFXGETPCAMATRIX Reads in pfx struct and maps results to a symmetrical matrix
% Note that in order for the function to work properly the force matrix has
% to be symmetrical, i.e. in the FDA Gromacs configuration file, Group1 must
% equal Group 2
%
% Example:
%   file='pathtofile\filename.pfx'
%   fid = pfaOpen(file);
%   pfx=pfxGetFrame(fid, 'vector', 'filter', [4101 4102 4103 4104 4105 4106 4107 4108 96553 112281 45186], 'AND')
%   [matrix, indices] = pfxGetMatrix(pfx, dimension of force matrix, 'offset', offset value);
    
    if length(fmdim) == 2
    n = fmdim(1);
    m = fmdim(2);
    else
        error
    end 
    

    argStruct = struct('offset', 0) ;
    argStruct = parseArgs(varargin, argStruct, {''}); % The cell specifies optional arguments that have to be included in argStruct

    % Fill in sparse array with force data
    uniqueIndices = unique(pfx.result(:,1:2))'-argStruct.offset;    
    % normalize indices
    minIndex = min(uniqueIndices)-1;
    globalIndices = pfx.result(:,1:2)-minIndex;
    

    % Checking sense of input
    [row, column]=size(pfx.result);
    if (column == 4)
        matrix = zeros((n*m),1);
        for j=1:length(globalIndices)
            x = globalIndices(j,1);
            y = globalIndices(j,2);
            index = xy2ix(x, y, n); % matrices for mechanical coupling always symmetrical?
            matrix(index, 1) = pfx.result(j,3);
        end
    elseif (column == 6)
        matrix = zeros((n*m),3);
        for j=1:length(globalIndices)
            x = globalIndices(j,1);
            y = globalIndices(j,2);
            index = xy2ix(x, y, n);
            matrix(index, :) = pfx.result(j,3:5);
        end
    end
    % Generate indices; size: (fmdim^2)x2
    indices = zeros(n*m,2);
    [indices(:,1), indices(:,2)] = ix2xy(1:(n*m),n);

end
