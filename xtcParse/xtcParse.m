% Mex-version of the gromacs xtc-parser
% Created: Hung Dien, 2012
% Last modified: 4.9.2013
% Version 1.0

function [pdb] = xtcParse(filename, varargin)
    argStruct = struct('start',   1, ...
                       'end',     0, ...
                       'skip',    1) ;

    argStruct = parseArgs(varargin, argStruct, {'verbose'; 'progress'; 'raw'}); % The cell specifies optional arguments that have to be included in argStruct


%
% check arguments
%
    if(exist('filename') == 0  || isnumeric(filename) == 1)
       disp('ERROR: Filename argument is not a string!');
       return;
    end

    if (isnumeric(argStruct.start) == 0)
       disp('ERROR: Start frame argument is not a number!');
       return;
    end

    if (isnumeric(argStruct.end) == 0)
       disp('ERROR: End frame argument is not a number!');
       return;
    end

   if(argStruct.end < 0)
        disp('ERROR: end frame has to be higher than zero!');
        return;
   end

    if(argStruct.start < 1)
        disp('ERROR: start frame has to be higher than zero!');
         return;
    end

    if(argStruct.start > argStruct.end && argStruct.end ~= 0)
         disp('ERROR: end frame has to be higher than start frame!');
          return;
    end


%
% end check arguments
%
   [coords, steps, time, precision, box, boxCoords] = gromacsMex2(argStruct.start, argStruct.end, filename);

  if (argStruct.skip == 1)
        pdb.coordinates = coords;
        pdb.timesteps = steps;
        pdb.time = time;
        pdb.precision = precision;
        pdb.box = boxCoords;
   elseif (argStruct.skip > 1)
        selection=(1:argStruct.skip:numel(time));
        pdb.coordinates = coords(:,:,selection);
        pdb.timesteps = steps(selection);
        pdb.time = time(selection);
        pdb.precision = precision(selection);
        pdb.box = boxCoords(selection);
   end

end


