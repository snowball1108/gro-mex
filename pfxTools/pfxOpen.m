function fid = pfaOpen(filename)

fid = fopen(filename,'r');

if (fid < 0) % Does the given filename exist?
    error('Could not open file "%s"', filename);
end

end 