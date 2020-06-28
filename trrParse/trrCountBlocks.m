function nframes = trrCountBlocks(fid)
%trrCountBlocks Count the number of blocks in a file. Returns the number of blocks. 
% After iterating over all data blocks the file pointer is set back the
% beginning of the data block. Runs this function before trrGetFrame, because 
% it will set the file pointer to the beginning of the file!
%   
%   Input parameter:
%       fid:    File identifier from *.trr file.
%
%   Output parameter:
%       trajectory: Struct containing timesteps, box vectors, coordinates, velocities and forces
%
%   Example for call: 
%       fid = trrOpen('$MATLAB-Repo$\tests\parseTRR\1MI5_peptide.pr.protein.trr')
%       nframes = trrCountBlocks(fid);
%       fclose(fid);
%   
%   See also:           
%       trrOpen.m, trrGetFrame.m 
%
%   Dependencies
%       none
%
%   Created:            Reiner J. Ribarics
%   Last modified:      $LastChangedDate: $
%   Version:            $Id: $

%-----------------------------START Subroutines-------------------------------%
    function int = trrInt()
        % Read an integer (size 4 Bytes)
        int = fread(fid, 1, 'int32');
    end 

    function int = trrIntFull(len)
        % Read an integer (size 4 Bytes)
        int = fread(fid, len, 'int32');
    end 

    function real = trrReal(prec)
        % Read a real single or double value depending on the given precision
        if prec == 4
            real = fread(fid, 1, 'single');
        elseif prec == 8
            real = fread(fid, 1, 'double');
        end
    end 

    function str = trrString()
        % String size is stored in the integer position right after the
        % version number. Reading string_size number of strings. 
        string_size = trrInt();
        str = fread(fid, string_size, 'char*1');
    end 

    function h_header = trrHeader()             
        % Read trrMagic number and version number
        header_start = trrIntFull(2);
        mn           = header_start(1);
        h_header.vn  = header_start(1);
      
        if (mn ~= trrMagic)
            h_header='Error: corrupted file';
            return;
            %error('File %s may be corrupted. Wrong magic number', filename);
        end        
        
        % Read 
        % h_header.vn = trrInt();
        
        % Reads string
        str = trrString();
        
        header_data = trrIntFull(13);
        
        % Read in some other header information
        h_header.ir_size   = header_data(1);
        h_header.e_size    = header_data(2);
        h_header.box_size  = header_data(3);
        h_header.vir_size  = header_data(4);
        h_header.pres_size = header_data(5);
        h_header.top_size  = header_data(6);
        h_header.sym_size  = header_data(7);
        h_header.x_size    = header_data(8);
        h_header.v_size    = header_data(9);
        h_header.f_size    = header_data(10);
        h_header.natoms    = header_data(11);
        h_header.step      = header_data(12);
        h_header.nre       = header_data(13);

        % Check if there are atoms
        if (h_header.natoms == 0)
            error('There are no atoms in %s.', filename);
        end 
        
        % Try to determine precision (float or double?)
        prec_x = h_header.x_size / (h_header.natoms * 3);
        prec_v = h_header.v_size / (h_header.natoms * 3);
        prec_f = h_header.f_size / (h_header.natoms * 3);
        
        if (prec_x == 0) && (prec_v ~= 0) && (prec_f ~= 0) 
            error('Failed to determine precision format.');
        else 
            h_header.prec_x = prec_x;
            h_header.prec_v = prec_v;
            h_header.prec_f = prec_f;
        end 
        
        % Read timestep and lambda value
        h_header.timestep = trrReal(h_header.prec_x);
        h_header.lambda   = trrReal(h_header.prec_x);
    end 

    function headsize = trrHeadsize()
        currentPosition = ftell(fid); 
        fseek(fid, 0, 'bof');
        trrHeader(); 
        headsize = ftell(fid); % and remember its end
        fseek(fid, currentPosition, 'bof');
    end 
%------------------------------END Subroutines--------------------------------%

trrMagic = 1993;
nframes  = 0;

% Jump to EOF to get file's bytesize   
fseek(fid, 0, 'eof');
EOF = ftell(fid);
fseek(fid, 0, 'bof'); % Set file position indicator back to beginning of file

% Filesize will get incremented by each jump and compared to EOF
filesize = 0;

% Assuming each header has the same size as the first one
cb_headsize = trrHeadsize();

while(1)
    % Determine size of header
    cb_header   = trrHeader(); 
    
    % Look out for corrupted header. If a file corruption is found return
    % error and corrupted frame
    if strcmp(cb_header, 'Error: corrupted file')
        nframes=sprintf('Error: corruption in frame %i', nframes);
        return
    end 

    % 3 + 3 + 3 Box = 9 additional entries
    jumpBytes   = 3 * cb_header.natoms * cb_header.prec_x * (numel(find([cb_header.prec_x cb_header.prec_v cb_header.prec_f]))) + (9 * cb_header.prec_x);
    filesize    = filesize + jumpBytes + cb_headsize;
    fseek(fid, jumpBytes, 'cof');

    nframes = nframes + 1;
    if (EOF - filesize == 0)
        break
    end
end 

% Set file position indicator back to beginning of file for trrGetFrame to
% work correctly
fseek(fid, 0, 'bof');

end