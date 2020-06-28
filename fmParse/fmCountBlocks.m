function frames = fmCountBlocks(fid, header)
%fmCountBlocks Count the number of blocks in a file. Returns the number of blocks. 
% After iterating over all data blocks the file pointer is set back the
% beginning of the data block. Runs this function before fmGetFrame!

fmMagic = -280480;
sizeof_int = header.intsize;
sizeof_real = header.realsize;
sizeof_char = 1;
frames = 0;

% Jump to EOF to get file's bytesize   
fseek(fid, 0, 'eof');
EOF = ftell(fid);
% Set file position indicator back to beginning of data blocks
% right after the header
fseek(fid, header.headsize, 'bof');
% Filesize will get incremented by each jump and compared to EOF
filesize = header.headsize;

while (1)
    if (EOF - filesize == 3)
        break
    end 
    iter  = fread(fid, 1, header.sizeof_int); % read timestep
    entries = fread(fid, 1, header.sizeof_int); % read number of entries

    offset = sizeof_int * header.fmdim + entries * (sizeof_int + sizeof_real + sizeof_char);
    % Three times four bytes (12) equals the current writestep, 
    % the number of entries and the magic number (aka NEW_ENTRY)
    filesize = filesize + uint64(offset + 3*sizeof_int); 
    if (filesize >= EOF)
        error('Calculated file size exceeds EOF. File may be corrupted!');
    end 
    fseek(fid, offset, 'cof');

    % Check for possible file corruption
    magic_number = fread(fid, 1, header.sizeof_int);
    if magic_number ~= fmMagic
        warning('File corrupted, only counted the first %i frames', frames);
    end
    frames = frames + 1;
end 
fseek(fid, header.headsize, 'bof'); % Reset the file position indicator to the beginning of the data blocks
end 