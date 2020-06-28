function pfx = pfxGetFrame(fid, varargin)
%pfaGetFrame reads one frame of PFA file
%   Input parameter: 
%       fid: file identifier of open PFA file
%       vector (optional): This flag tells the function that the PFA file contains 
%       three dimensional vectors
%       filter (optional) : 1xN array containing atom indices for output
%       AND, OR (optional): The PFA file contains two columns containing atom indices
%       between atoms i and atoms j. The filter option is run for the column
%       of indices i and for the column of indices j. In the last step both
%       results are combined on a logical AND or OR. If you use AND you 
%       will get only interactins between the set of atoms you defined in 
%       filter. If you use OR you will get all interactions between the 
%       atoms you defined in filter and any other atoms. 
%   
%   Output parameter:
%       data: Nx4 array containing pairwise interactions in the following
%       form - [atom i] [atom j] [force] [interaction type]
%
%   Example for call:   
%       fid     = pfaOpen(filename);
%       frame0  = pfaGetFrame(fid);
%       frame1  = pfaGetFrame(fid);
%       ...
%       frameN  = pfaGetFrame(fid);
%       fclose(fid);
%
%   See also:           
%       fmGetFrame
%
%   Dependencies:       
%       none
%
%   Created:            Reiner J. Ribarics
%   Last modified:      $LastChangedDate: $

argStruct = struct('filter',  [], 'vector', false, 'AND', false, 'OR', false, 'single', false) ; 
argStruct = parseArgs(varargin, argStruct, {'vector', 'AND', 'OR', 'single'}); % The cell specifies optional arguments that have to be included in argStruct

if (isempty(argStruct.filter) == false)
    if (argStruct.AND == true) & (argStruct.OR == true)
        error('Both flags "AND" and "OR" are set on!')
    elseif (argStruct.AND == false) & (argStruct.OR == false)
        error('Both flags "AND" and "OR" are set off!')    
    end
end 


chunksize   = 500000; % Set chunksize above estimated number of entries per frame
blocksize   = chunksize;
data        = [];

% Get frame start
% The header consists of a single line containing "Frame N", with N being the frame number  
frameHeader = fgets(fid); 

% Skip header if present
if strfind(strtrim(frameHeader), 'pairwise')
    fgets(fid); 
end 

% When end of file is reached, fgets will return -1
if (frameHeader == -1)
    pfx.result = -1;   
    return
end

timestep        = textscan(frameHeader, 'frame %f');
pfx.timestep    = cell2mat(timestep);

% Read chunks of data 
if (argStruct.vector) 
    while (blocksize == chunksize)
        chunk = textscan(fid, '%f %f %f %f %f %f', chunksize);
        if (((numel(chunk{1}) + numel(chunk{2}) + numel(chunk{3}) + numel(chunk{4}) + numel(chunk{5}) + numel(chunk{6}))/6) ~= numel(chunk{1}))
            error('Inconsistant data! File may be corrupted');
        end 
        data        = [data; horzcat(chunk{1}, chunk{2}, chunk{3}, chunk{4}, chunk{5}, chunk{6})];
        blocksize   = numel(chunk{1});
    end
else 
    while (blocksize == chunksize)
        chunk       = textscan(fid, '%f %f %f %f', chunksize);
        data        = [data; horzcat(chunk{1}, chunk{2}, chunk{3}, chunk{4})];
        blocksize   = numel(chunk{1});
    end
end 

if (sum(sum(isnan(data))) ~= 0)
    error('Invalid data detected! Are you trying to read vector based data without the "vector" flag option?')
end

data(:,1:2) = data(:,1:2) + 1; % correct for indexing

% Filter interactions
if (isempty(argStruct.filter)) 
    if (argStruct.single)
        pfx.result  = single(data);
    else 
        pfx.result  = data;
    end
else 
    fi          = ismember(data(:,1), argStruct.filter);
    fj          = ismember(data(:,2), argStruct.filter);
    if (argStruct.OR == true)
        fij         = (fi|fj);
    elseif (argStruct.AND == true)
        fij         = (fi&fj);
    end 
    
    fdata = data(fij,:);
    if (argStruct.single)
        pfx.result  = single(fdata);
    else 
        pfx.result  = fdata;
    end
end

end 

%{
% Interaction types
                            Hex     Decimal
#define PF_INTER_NONE		0       0
#define PF_INTER_BOND		0x01    1
#define PF_INTER_ANGLE		0x02    2
#define PF_INTER_DIHEDRAL	0x04    4 
#define PF_INTER_POLAR		0x08    8
#define PF_INTER_COULOMB	0x10    16
#define PF_INTER_LJ         0x20    32
#define PF_INTER_NB14		0x40    64

Any combination of numbers, assuming that it will be used only once, can be
decomposed to the interactions listed above. 

Example 
    48 = 32 + 16   -> Lennard Jones + Coulmomb 
    7  = 2 + 4 + 1 -> Angle + Dihedral + Bond 

% defines for easy tests of bonded/nonbonded/all interactions */
#define PF_INTER_BONDED		PF_INTER_BOND + PF_INTER_ANGLE + PF_INTER_DIHEDRAL
#define PF_INTER_NONBONDED	PF_INTER_COULOMB + PF_INTER_LJ + PF_INTER_NB14
#define PF_INTER_ALL		PF_INTER_BONDED + PF_INTER_POLAR + PF_INTER_NONBONDED

%}