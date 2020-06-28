function stats = pfxCheck(filename, varargin)
% PFXCHECK Checks file integrity of PF2 Gromacs force files (*.pfa and *.pfr)
%
%
% Types of atom atom interactions
%
% #define PF_INTER_NONE		0       0
% #define PF_INTER_BOND		0x01    1
% #define PF_INTER_ANGLE	0x02    2
% #define PF_INTER_DIHEDRAL	0x04    4 
% #define PF_INTER_POLAR	0x08    8
% #define PF_INTER_COULOMB	0x10    16
% #define PF_INTER_LJ       0x20    32
% #define PF_INTER_NB14		0x40    64


argStruct = struct('vector', false) ; 
argStruct = parseArgs(varargin, argStruct, {'vector'}); % The cell specifies optional arguments that have to be included in argStruct


fid = pfaOpen(filename);

stats.nAtoms    = [];
stats.nFrames   = 0;
stats.intType   = []; 
stats.nInt      = []; 

while 1==1
    stats.nFrames = stats.nFrames + 1;
    fprintf('Reading frame %i\n', stats.nFrames)
    
    if argStruct.vector == true
        data = pfaGetFrame(fid, 'vector');
    else
        data = pfaGetFrame(fid);
    end 
    
    % When EOF is reached result is set to -1 by pfxGetFrame 
    if (data.result == -1)
        break
    end
    
    % Check if data are without NaN or Inf entries
    if ( sum(sum(isnan(data.result))) ~= 0 || sum(sum(isinf(data.result))) ~= 0 )
        error('File corruption detected at frame %i', stats.nFrames)
    end
    
    stats.nInt      = [stats.nInt, length(data.result)];
    stats.intType   = unique([stats.intType, unique(data.result(:,end))]);
    stats.nAtoms    = [stats.nAtoms, numel(unique([data.result(:,1); data.result(:,2)]))];
end

fclose(fid);

fprintf('\n')
fprintf('Frames:                        %i\n', stats.nFrames)
fprintf('Interation types:              %s\n', num2str(reshape(stats.intType, 1,numel(stats.intType))))
fprintf('Mean number of interactions:   %i\n', round(mean(stats.nInt)))
fprintf('Mean number of atoms:          %i\n', round(mean(stats.nAtoms)))
fprintf('\n')

end 