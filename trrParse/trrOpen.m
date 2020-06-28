function fid = trrOpen(filename)
%
%
trrMagic          = 1993;
fid = fopen(filename,'r', 'b'); % Open as big-endian (default)

if (fid < 0) % Does the given filename exist?
    error('Could not open file "%s"', filename);
end

% Check trrMagic number 
mn = fread(fid, 1, 'int32');
if (mn ~= trrMagic) % If trrMagic numer does not match try little-endian
    fclose(fid);
    fid = fopen(filename,'r', 'l');
    if (fid < 0) % Does the given filename exist?
        error('Could not open file "%s"', filename);
    end

    % Read trrMagic bit again
    mn = fread(fid, 1, 'int32');
    if (mn ~= trrMagic)
        error('File %s may be corrupted. Wrong magic number', filename);
    end 
end 
fseek(fid, 0, 'bof'); % Set file position indicator to beginning of file ('bof')

end 