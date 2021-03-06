function [fid, header] = fmOpen(filename)
%parseFM Summary of this function goes here
%   The function reads binary files containing a N x N x F 
%   force matrix generated by FDA gromacs, where N is the number of 
%   atoms in the system and F the number of frames. Detailed information 
%   about the header and file format structure is found in the code
%   within this function.
%   
%   Input parameter: 
%       filter: List of atoms to parse forces from (e.g. [1 2 3 4 5 ... N] or [1:1:N] )
%   
%   Output parameter:
%       pforce: Struct containing dimension of force matrix (fmdim), the
%       index matrix (index_matrix), timesteps and the results. The results comprise 
%       the index array, force(pforce)- and interaction matrix(interaction).
%
%   Example for call:   
%       [fid, header]=fmOpen(filename)
%       nframes = fmCountBlocks(fid, header)
%       frame1  = fmGetFrame(fid, header)
%       frame2  = fmGetFrame(fid, header)...
%       fclose(fid)
%
%   See also:           
%       parseTRR.m, parseFM.m
%
%   Dependencies:       
%       parseArgs.m, progressbar.m 
%
%   Created:            Reiner J. Ribarics
%   Last modified:      $LastChangedDate: $
%   Version:            $Id: $


%%% START Subroutines %%%
    function fid = fmOpen(filename)
        fid = fopen(filename,'r');

        if (fid < 0) % Does the given filename exist?
            error('Could not open file "%s"', filename);
        end
    end

function fm_header = fmHeader(fid)    
        % FIXME: rewrite header handling to allow for flexible headers! 
        % Read 13 lines that correspond to the header
        line                    = fgets(fid);
        forcemat_version        = cell2mat(textscan(fgets(fid), '; Forcemat version %8.1f'));
        content                 = strtrim(fgets(fid));
        index_matrix_dimension  = cell2mat(textscan(fgets(fid), '; Matrix dimension %f'));
        version                 = fgets(fid);
        groupname               = textscan(fgets(fid), 'groupname=%s');
        write_frequency         = cell2mat(textscan(fgets(fid), 'writefreq=%f'));
        nsteps                  = cell2mat(textscan(fgets(fid), 'nsteps=%f'));
        sysanr                  = cell2mat(textscan(fgets(fid), 'sysanr=%f'));
        fmdim                   = cell2mat(textscan(fgets(fid), 'fmdim=%f'));
        intsize                 = cell2mat(textscan(fgets(fid), 'intsize=%f'));
        realsize                = cell2mat(textscan(fgets(fid), 'realsize=%f'));
        line                    = fgets(fid);
        headsize                = double(ftell(fid));

        % Set int and real size for reading binary blocks            
        if (intsize == 4)
            sizeof_int = '*int';
        elseif (intsize == 8)
            sizeof_int = '*int64';
        end
        if (realsize == 4) 
            sizeof_real = '*single';
            precision   = 'single';
        elseif (realsize == 8)
            sizeof_real = '*double';
            precision   = 'double';
        end 

        % Group header information
        fm_header = struct('forcemat_version', forcemat_version,... 
                           'index_matrix_dimension', index_matrix_dimension,...
                           'groupname' ,groupname,...
                           'write_frequency', write_frequency,...
                           'nsteps', nsteps,...
                           'sysanr', sysanr,...
                           'fmdim', fmdim,...
                           'intsize', intsize,...
                           'realsize', realsize,...
                           'sizeof_int', sizeof_int,...
                           'sizeof_real', sizeof_real,...
                           'headsize', headsize,... 
                           'content', content,...
                           'precision', precision);
        
        % Check forcemat version
        if ((fm_header.forcemat_version) ~= 1.2 && (fm_header.forcemat_version ~= 1.5))
            error('Forcemat version %1.1f is not supported', fm_header.forcemat_version);
        end 
    end 
%%% END Subroutines %%%

% Check for necessary functions
if (exist('parseArgs.m', 'file') == 0) || (exist('progressbar.m', 'file') == 0)
    error('Your are missing either parseArgs.m or progressbar.m. You can download it on "http://www.mathworks.com"');
end 

fid     = fmOpen(filename);   
header  = fmHeader(fid);

% Check for integer overflow
maxIntegerValue = 4294967295;
if (header.fmdim > sqrt(maxIntegerValue))
    warning('Warning: your system size very large. Integer overflow will cause corrupt atom indexing!');
end

end 