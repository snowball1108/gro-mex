function data = datParse(filename, precision)
%DATPARSE Parses Gromacs *.dat files
%
%   In-params           Description                                     Default value
%   precision           For single precision set to 4
%                       For double precision set to 8
%   
%   Out-params:
%   data                Square matrix 
% 
%   Example for call:   data = datParse('data.dat', 4)
%   See also:           
%   Dependencies:       none
%
%   Version:            
%   Author:             Reiner J. Ribarics, Medical University of Vienna

fid = fopen(filename,'r');

if (fid < 0) % Does the given filename exist?
    error('Could not open file "%s"', filename);
end

fseek(fid, 0, 'eof');
fileSize = ftell(fid);
frewind(fid);
    
if (precision == 4)
    nEntries=(fileSize/4);
    if (rem(sqrt(nEntries), 1) ~= 0)
        error('Data does not contain square matrix')
    end 
    data =  fread(fid, nEntries, 'single');
elseif (precision == 8)
    nEntries=(fileSize/8);
    if (rem(sqrt(nEntries), 1) ~= 0)
        error('Data does not contain square matrix')
    end
    data =  fread(fid, nEntries, 'double');
end 

fclose(fid);
data = reshape(data, sqrt(nEntries), sqrt(nEntries));

end 