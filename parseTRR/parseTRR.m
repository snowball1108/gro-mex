function [trajectory] = parseTRR(filename, varargin)
%PARSETRR Parses GROMACS *.trr file and returns coordinates, velocities and
%forces. Additional timesteps in ps and box size are read. 
%
%   In-params           Description                                     Default value
%   start               First frame                                     1
%   end                 Last frame                                      nframes (i.e. all frames)
%   skip                Skip every nth frame                            1 (no frames will be skipped) 
%   coor                Read coordinates                                True
%   vel                 Read velocities                                 False
%   force               Read forces and                                 False
%   
%	Notes: Counting frames starts with 1 and ends with the total number of frames
%   (eg. 1 to 51 when frameCount = 51). Per default, velocities and forces
%   will not be read. Use the flag options 'vel' or 'force' to enable
%   parsing. 
%
%   Out-params:
%   trajectory:         Struct containing timesteps, box vectors, coordinates and if enabled velocities and forces 
%
%   Example for call:   trajectory = parseTRR('testInputParseTRR.trr', 'start', 5, 'end', 25, 'skip', 2);   
%   See also:           trrOpen, trrGetFrame, trrCountBlocks
%   Dependencies:       parseArgs.m, progressbar.m
%
%   Version:            11/2015
%   Author:             Reiner J. Ribarics, Medical University of Vienna

%-----------------------------START Subroutines-------------------------------%

    function fid = trrOpen(filename)
        fid = fopen(filename,'r', 'b'); % Open as big-endian (default)

        if (fid < 0) % Does the given filename exist?
            error('Could not open file "%s"', filename);
        end
        
        % Check magic number 
        mn = fread(fid, 1, 'int32');
        if (mn ~= trrMagic) % If magic number does not match, try little-endian format
            fclose(fid);
            fid = fopen(filename,'r', 'l');
            if (fid < 0) % Does the given filename exist?
                error('Could not open file "%s"', filename);
            end
            
            % Read magic number again
            mn = fread(fid, 1, 'int32');
            if (mn ~= trrMagic)
                error('File %s may be corrupted. Wrong magic number', filename);
            end 
        end 
        fseek(fid, 0, 'bof'); % Rewind file position indicator (set to beginning of file ('bof'))
    end 

    function int = trrInt()
        % Read an integer (size 4 Bytes) number
        int = fread(fid, 1, 'int32');
    end 

    function int = trrIntFull(len)
        % Read integer (size 4 Bytes) numbers of length 'len'
        int = fread(fid, len, 'int32');
    end 

    function h_header = trrHeader()             
        % Read magic number
        mn = trrInt();
        if (mn ~= trrMagic)
            error('File %s may be corrupted. Wrong magic number in frame %i', filename, cb_nframes+1);
        end        
        h_header.mn = mn;
        
        % Read version number
        h_header.vn = trrInt();
        
        % Read string
        trrString();
        
        header_data = trrIntFull(13);
        
        % Read more header information
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

        % Check if there are atoms saved in the file
        if (h_header.natoms == 0)
            error('There are no atoms in %s.', filename);
        end
        
        % Try to determine precision (float or double?)
        prec_x = h_header.x_size / (h_header.natoms * 3);
        prec_v = h_header.v_size / (h_header.natoms * 3);
        prec_f = h_header.f_size / (h_header.natoms * 3);
        
        if (prec_x == 0) && (prec_v ~= 0) && (prec_f ~= 0) 
            error('Failed to determine precision format. There are no coordinates save in your file.');
        else 
            h_header.prec_x = prec_x;
            h_header.prec_v = prec_v;
            h_header.prec_f = prec_f;
        end 
        
        % Read timestep and lambda value
        h_header.timestep = trrReal(h_header.prec_x);
        h_header.lambda   = trrReal(h_header.prec_x);
    end 

    function real = trrReal(prec)
        % Read a real single or double value depending on the given precision
        if prec == 4
            real = fread(fid, 1, 'single');
        elseif prec == 8
            real = fread(fid, 1, 'double');
		else
			error('Unknown precision.');
        end
    end 

    function str = trrString()
        % String size is stored in the integer position right after the
        % version number. 
		% Read string_size number
        string_size = trrInt();
		% Read a number of strings
        str = fread(fid, string_size, 'char*1');
    end 

    function vec = trrRvector(prec)
        x = trrReal(prec);
        y = trrReal(prec);
        z = trrReal(prec);
        vec = [x y z];
    end  

    function vec = trrRvectorFull(ts_natoms, prec)
        % Read a real single or double value depending on the given precision
        if prec == 4
            real  = fread(fid, ts_natoms*3, 'single');
        elseif prec == 8
            real  = fread(fid, ts_natoms*3, 'double');
        else
			error('Unknown precision.');
        end
        vec = rot90(reshape(real,3,numel(real)/3,1)',2); 
    end  
    
    function [nframes, cb_blocksize] = trrCountBlocks()
        nframes 	= 0;
        cb_nframes 	= 0;

        % Jump to end of file to get file bytesize   
        fseek(fid, 0, 'eof');
        EOF = ftell(fid);
        
        % Set file position indicator back to beginning of file
        fseek(fid, 0, 'bof');
        
        % Filesize will get incremented by each jump and compared to EOF
        filesize = 0;
        
		% Blocksize contains the size of each timestep data block
        % Blocksize will be initialised generously
        cb_blocksize = zeros(50000,1);
        
        while(1)
            % Determine size of header
            cb_header   = trrHeader(); 
            cb_headsize = trrHeadsize();
            cb_nframes 	= cb_nframes + 1;

            % DEBUG
            % fprintf('Frame: %i | Magic number: %i | Filesize start: %i | ', cb_nframes, cb_header.mn, filesize)
			
            filesize_start = filesize;
            %jumpBytes   = 3 * cb_header.natoms * cb_header.prec_x * (numel(find([cb_header.prec_x cb_header.prec_v cb_header.prec_f]))) + (9 * cb_header.prec_x);
            jumpBytes = (3 * cb_header.natoms * cb_header.prec_x) + ...
                        (3 * cb_header.natoms * cb_header.prec_v) + ...
                        (3 * cb_header.natoms * cb_header.prec_f) + ...
                        (9 * cb_header.prec_x); % 3 + 3 + 3 Box = 9 additional entries
            
			filesize    = filesize + jumpBytes + cb_headsize;
            
			% Jump to header of the next block
			fseek(fid, jumpBytes, 'cof');
            
            % DEBUG
            % fprintf('Filesize end: %i |Filesize diff: %i\n', filesize, (filesize-filesize_start))
                 
            nframes = nframes + 1;
			
			% Difference between current filesize and filesize before the jump equals to current block size
            cb_blocksize(nframes) = filesize-filesize_start;
			
			% Stop loop of filesize reaches end of file
            if (EOF - filesize == 0)
                break
            end 
        end 
        
        % Trim trailing zeros
        cb_blocksize(cb_blocksize==0)=[];
		
        % Set file position indicator back to beginning of file
        fseek(fid, 0, 'bof');
    end

    function headsize = trrHeadsize()
            currentPosition = ftell(fid); 
            fseek(fid, 0, 'bof');
            trrHeader();
            headsize = ftell(fid); % and remember its end
            fseek(fid, currentPosition, 'bof'); 
    end 

    function trrJumpToFrame(targetFrame)
        % Target frame has to be decreased by 1
        targetFrame = targetFrame - 1;
        if (targetFrame == 0)
            fseek(fid, 0, 'bof'); % Rewind to beginning of file
            return
        end
        
        % Jump to end of file to get file bytesize   
        fseek(fid, 0, 'eof');
        EOF = ftell(fid);
        
        % Set file position indicator back to beginning of file
        fseek(fid, 0, 'bof');
        
        % Filesize will get incremented by each jump and compared to EOF
        filesize    = 0;
        jtf_frames  = 0;
        
        for i = 1:targetFrame
            % Read header
            jtf_header = trrHeader(); 

            % Jumps one frame 
            % Calculate  the data block size between current and next header
            jumpBytes = (3 * jtf_header.natoms * jtf_header.prec_x) + ...
                        (3 * jtf_header.natoms * jtf_header.prec_v) + ...
                        (3 * jtf_header.natoms * jtf_header.prec_f) + ...
                        (9 * jtf_header.prec_x); % 3 + 3 + 3 Box = 9 additional entries
            
            % 
            filesize = filesize + jumpBytes + headsize;

            if filesize > EOF
                error('Function trrJumpToFrame exceeded filesize.');
            end 

            fseek(fid, jumpBytes, 'cof');       
            jtf_frames = jtf_frames + 1;
        end
    end

    function timestep = trrTimestep()
        % Read header
        ts_header       = trrHeader();        
        timestep.step   = ts_header.step;
        
        % Read box size
        ts_box(1,:) = trrRvector(ts_header.prec_x);
        ts_box(2,:) = trrRvector(ts_header.prec_x);
        ts_box(3,:) = trrRvector(ts_header.prec_x);
        timestep.box = ts_box;
        
        %%% COORDINATES %%%
        % Read coordinates
        if (ts_header.x_size > 0) && (argStruct.coor == true)
           ts_coordinates      = zeros(ts_header.natoms, 3);
           ts_coordinates(:,:) = trrRvectorFull(ts_header.natoms, ts_header.prec_x);
           
           if (ts_header.prec_x == 4)
                timestep.coordinates = single(ts_coordinates);
           elseif (ts_header.prec_x == 8)
                timestep.coordinates = ts_coordinates;
           end    

        % Skip coordinates   
        elseif (ts_header.x_size > 0) && (argStruct.coor == false)
                jumpBytes =  3 * ts_header.natoms * ts_header.prec_x; 
                fseek(fid, jumpBytes, 'cof');
        end    
        
        %%% VELOCITIES %%%
        % Read velocities
        if (ts_header.v_size > 0) && (argStruct.vel == true)
            ts_velocities       = zeros(ts_header.natoms, 3);
            ts_velocities(:,:)  = trrRvectorFull(ts_header.natoms, ts_header.prec_v);   
            
            if (ts_header.prec_v == 4)
                timestep.velocities = single(ts_velocities);
            elseif (ts_header.prec_v == 8)
                 timestep.velocities = ts_velocities;
            end 
        % Skip velocities if necessary   
        elseif (ts_header.v_size > 0) && (argStruct.vel == false)
            jumpBytes =  3 * ts_header.natoms * ts_header.prec_v;
            fseek(fid, jumpBytes, 'cof'); 
        end 
        
        %%% FORCES %%%
        % Read forces
        if (ts_header.f_size > 0) && (argStruct.force == true)
            ts_forces       = zeros(ts_header.natoms, 3);
            ts_forces(:,:)  = trrRvectorFull(ts_header.natoms, ts_header.prec_f);   
            
            if (ts_header.prec_f == 4)
                timestep.forces = single(ts_forces);
            elseif (ts_header.prec_f == 8)
                 timestep.forces = ts_forces;
            end 
        % Skip force if necessary
        elseif (ts_header.f_size > 0) && (argStruct.force == false)
            jumpBytes =  3 * ts_header.natoms * ts_header.prec_f;
            fseek(fid, jumpBytes, 'cof');             
        end
    end

    function trrJumpOneFrame()
        jof_header = trrHeader(); 
        jumpBytes = (3 * jof_header.natoms * jof_header.prec_x) + ...
                    (3 * jof_header.natoms * jof_header.prec_v) + ...
                    (3 * jof_header.natoms * jof_header.prec_f) + ...
                    (9 * jof_header.prec_x); % 3 + 3 + 3 Box = 9 additional entries      
        fseek(fid, jumpBytes, 'cof');
    end       
%------------------------------END Subroutines--------------------------------%

% Check for necessary functions
if (exist('parseArgs.m', 'file') == 0) || (exist('progressbar.m', 'file') == 0)
    error('Your are missing either parseArgs.m or progressbar.m. You can download it on "http://www.mathworks.com"');
end

% Initialize global variables
global            cb_nframes;
trrMagic        = 1993;
fid             = trrOpen(filename);
headsize        = trrHeadsize();
%[path file ext] = fileparts(filename);


% Input handling
argStruct = struct('start',     1,       ...
                   'end',       0,       ...
                   'skip',      1,       ...
                   'verbose',   false,   ...
                   'progress',  false,   ...
                   'vel',       false,   ...
                   'force',     false,   ...
                   'coor',      true); 
               
argStruct = parseArgs(varargin, argStruct, {'verbose'; 'progress'; 'vel'; 'force'; 'coor'});

% Get number of frames from user input or count it
[nframes, ~] = trrCountBlocks();

% If end frame = 0 then no user input was given and end frame will be set
% equal to number of frames
if (argStruct.end == 0)
    argStruct.end = nframes;
end

% FIXME: check if value of argStruct.trim does make sense!
% Checking sense of input

% Is your value for end frame greater or equal than 0?
if (sign(argStruct.end) < 0)
    error('Check your input of end frame!');
end

% Is your value for start frame greater or equal than 1?
if (sign(argStruct.start) < 1)
    error('Check your input of start frame!');
end

% Is your value for skip frame greater or equal than 1?
if (sign(argStruct.skip) < 0)
    error('Check your input of skip frame!');
end

% Is start frame greater than end frame?
if argStruct.start > argStruct.end
    error('Start frame is greater than end frame!'); 
end 

% Is end frame greater than the number of frames?
if argStruct.end > nframes
    error('The trajectory contains %i nframes. You entered %i!', nframes, argStruct.end);
end
    

% If end frame is 
readnFrames = (argStruct.end - argStruct.start + 1);

 
% Retrieving header information
currentPos  = ftell(fid); 
header      = trrHeader();
fseek(fid, currentPos, 'bof'); % Rewind to currentPos 


% Output data to command line if verbose mode is selected
if (argStruct.verbose)
    fprintf('; Trajectory contains %i frames\n', nframes);
    fprintf('; Trajectory contains %i atoms per frame\n', header.natoms);
end 


% Allocating memory
if (argStruct.trim > 0)
    coordinates = zeros(argStruct.trim, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    if (header.v_size > 0)
        velocities = zeros(argStruct.trim, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    end 
    if (header.f_size > 0)
        forces = zeros(argStruct.trim, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    end 
else
    coordinates = zeros(header.natoms, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    if (header.v_size > 0)
        velocities = zeros(header.natoms, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    end 
    if (header.f_size > 0)
        forces = zeros(header.natoms, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));
    end 
end


timesteps   = zeros(ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip), 1);
box         = zeros(3, 3, ceil(readnFrames/argStruct.skip) - mod(readnFrames, argStruct.skip));


% Initialize progressbar
if (argStruct.progress)
    set(0,'defaulttextinterpreter','none'); % Turn off Tex interpreter for better visualisation
    progressbar(['Reading trajectory ', file]);
end 


% Jump to start frame
trrJumpToFrame(argStruct.start);


% Parsing nframes
ii = 1;
for counter = 1:readnFrames
    if (mod(counter, argStruct.skip) == 0) % Skipping frames
        timestep            = trrTimestep();
        coordinates(:,:,ii) = timestep.coordinates;
        box(:,:,ii)         = timestep.box;
        timesteps(ii)       = timestep.step;
        
        if (isfield(timestep, 'velocities'))
            velocities(:,:,ii)  = timestep.velocities;
        end 
        if (isfield(timestep, 'forces'))
            forces(:,:,ii)      = timestep.forces;
        end
        
        % Calculate percentage of job and passing it to progessbar
        if (argStruct.progress) && (mod(counter, ceil(readnFrames/100))==1)
            progressbar(counter/readnFrames);
        end
        
        ii = ii + 1;
    else
        trrJumpOneFrame();
    end
end
progressbar(1);


% Group variables into struct object
trajectory.timesteps    = timesteps;
trajectory.box          = box;

if (isfield(timestep, 'coordinates'))
    trajectory.coordinates  = coordinates;
end 
if (header.v_size > 0)
    trajectory.velocities = velocities;
end 
if (header.f_size > 0)
    trajectory.forces = forces;
end 

fclose(fid);

end

%{
Header information
-----------------------------------------------------------------------------
Header
-----------------------------------------------------------------------------
datatype    count       description
int         1           mn          Magic number (1993)
int         1           vn          Version number
int         1           strLen      Length of following string
char        strLen      str         Hell knows what is stored in here
int         1           ir_size     -
int         1           e_size      -
int         1           box_size    Size of box
int         1           vir_size    -
int         1           pres_size   -
int         1           top_size    -
int         1           sym_size    -
int         1           x_size      Number of stored coordinates
int         1           v_size      Number of stored velocities
int         1           f_size      Number of stored forces
int         1           natoms      Number of atoms
int         1           step        Step number
int         1           nre         - 
int         1           timestep    Timestep  
int         1           lambda      Lamda
-----------------------------------------------------------------------------

Data Block organisation
-----------------------------------------------------------------------------
datatype    count       description
real         3x3        Box vectors
real         3xnatoms   Coordinates
-----------------------------------------------------------------------------
optional:
real         3xnatoms   Velocities
real         3xnatoms   Forces
-----------------------------------------------------------------------------
%}