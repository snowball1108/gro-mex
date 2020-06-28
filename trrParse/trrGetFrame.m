function trajectory = trrGetFrame(fid, varargin)
%trrGetFrame
% Function reads one frame of a *.trr file and gives back coordinates and,
% if desired, velocities and forces as well. 
%   
%   Input parameter:                                        Default value   
%       coor:   Will read coordinates and allocate memory.      True
%       vel:    Will read velocities and allocate memory.       False
%       force:  Will read forces and allocate memory.           False
%   
%   Output parameter:
%       trajectory: Struct containing timesteps, box vectors, coordinates, velocities and forces
%
%   Example for call: 
%       fid = trrOpen('$MATLAB-Repo$\tests\parseTRR\1MI5_peptide.pr.protein.trr')
%       trajectory = trrGetFrame(fid, 'trim', 10);   
%       fclose(fid);
%   
%   See also:           
%       trrOpen.m, trrCountbBlocks.m 
%
%   Dependencies
%       parseArgs.m
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

    function h_header = trrHeader()             
        % Read trrMagic number
        mn = trrInt();
        if (mn ~= trrMagic)
            error('File %s may be corrupted. Wrong magic number', filename);
        end        
        
        % Read version number/
        h_header.vn = trrInt();
        
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

    function vec = trrRvector(prec)
        x = trrReal(prec);
        y = trrReal(prec);
        z = trrReal(prec);
        vec = [x y z];
    end  

    function vec = trrRvectorFull(ts_natoms, prec)
        
        % Read a real single or double value depending on the given precision
        if prec == 4
            real = fread(fid, ts_natoms*3, 'single');
            real1 = reshape(real, 3, numel(real)/3, 1); 
            real2 = rot90(real1, 3);
            real3 = fliplr(real2);
        elseif prec == 8
            real = fread(fid, ts_natoms*3, 'double');
            real1 = reshape(real, 3, numel(real)/3, 1); 
            real2 = rot90(real1, 3);
            real3 = fliplr(real2);
        end
       vec = real3;
    end  

    function timestep = trrTimestep()
        ts_header       = trrHeader();        
        timestep.step   = ts_header.step;
        
        % Read the box size
        ts_box(1,:) = trrRvector(ts_header.prec_x);
        ts_box(2,:) = trrRvector(ts_header.prec_x);
        ts_box(3,:) = trrRvector(ts_header.prec_x);
        timestep.box = ts_box;
        
        % FIXME: implement trim function for velocities and forces
        % Read coordinates
        if (ts_header.x_size > 0) && (argStruct.coor == true)
           if (argStruct.trim > 0) 
               natoms   = argStruct.trim;
           else
               natoms   = ts_header.natoms;
           end 
           diff = ts_header.natoms - natoms;

           ts_coordinates      = zeros(natoms, 3);
           ts_coordinates(:,:) = trrRvectorFull(natoms, ts_header.prec_x);
           
           if (ts_header.prec_x == 4)
                timestep.coordinates = single(ts_coordinates);
           elseif (ts_header.prec_x == 8)
                timestep.coordinates = ts_coordinates;
           end    
           % Skipping rest of coordinates if necessary
           if (diff ~= 0)
               jumpBytes =  3 * diff * ts_header.prec_x;
               fseek(fid, jumpBytes, 'cof');
           end 
        elseif (ts_header.x_size > 0) && (argStruct.coor == false)
                % Skip coors
                jumpBytes =  3 * ts_header.natoms * ts_header.prec_x; 
                fseek(fid, jumpBytes, 'cof'); 
        end    
        
        % Read velocities
        if (ts_header.v_size > 0) && (argStruct.vel == true)
            ts_velocities       = zeros(ts_header.natoms, 3);
            ts_velocities(:,:)  = trrRvectorFull(ts_header.natoms, ts_header.prec_v);   
            
            if (ts_header.prec_v == 4)
                timestep.velocities = single(ts_velocities);
            elseif (ts_header.prec_v == 8)
                 timestep.velocities = ts_velocities;
            end 
        elseif (ts_header.v_size > 0) && (argStruct.vel == false)
            % Skip vel
            jumpBytes =  3 * ts_header.natoms * ts_header.prec_v;
            fseek(fid, jumpBytes, 'cof'); 
        end 
        
        % Read forces
        if (ts_header.f_size > 0) && (argStruct.force == true)
            ts_forces       = zeros(ts_header.natoms, 3);
            ts_forces(:,:)  = trrRvectorFull(ts_header.natoms, ts_header.prec_f);   
            
            if (ts_header.prec_f == 4)
                timestep.forces = single(ts_forces);
            elseif (ts_header.prec_f == 8)
                 timestep.forces = ts_forces;
            end 
        elseif (ts_header.f_size > 0) && (argStruct.force == false)
            % Skip force
            jumpBytes =  3 * ts_header.natoms * ts_header.prec_f;
            fseek(fid, jumpBytes, 'cof');             
        end
    end    
%------------------------------END Subroutines--------------------------------%

trrMagic = 1993;

% Input handling
argStruct = struct('trim',  0, ...
                   'vel',   false, ....
                   'force', false, ...
                   'coor',  true); 
argStruct = parseArgs(varargin, argStruct, {'vel'; 'force'; 'coor'});

% Parsing frame
timestep            = trrTimestep;
trajectory.timestep = timestep.step;
trajectory.box      = timestep.box;

if (isfield(timestep, 'coordinates'))
    trajectory.coordinates = timestep.coordinates;
end 

if (isfield(timestep, 'velocities'))
    trajectory.velocities = timestep.velocities;
end 

if (isfield(timestep, 'forces'))
    trajectory.forces = timestep.forces;
end 

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