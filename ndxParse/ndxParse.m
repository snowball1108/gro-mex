function topology = ndxParse(filename)
%NDXPARSE Reads in parameters from *.ndx Gromacs topology file
%
%   In-params
%   fileName            Path to text file containing parameters to read
%   
%   Out-params:
%   params              A struct object containing key/value pairs
%
%   Example for call:   params = ndxParse('C:\test.ndx')
%   See also:           
%   Dependencies:       none
%
%   Version:            $Id$

fid = fopen(filename);
d=fread(fid, 'char=>char');
d=d';
fclose(fid);

str = regexprep(d,'\s+',' ');

groupsi = [find(ismember(str, '['))', find(ismember(str, ']'))'];
len = length(groupsi);

for i = 1:len
    group = str((groupsi(i, 1)+2):(groupsi(i,2)-2));
    group(regexp(group,'[-,+]'))='_';
    
    if i == len
        block = str((groupsi(i, 2)+2):end);
    else
        block = str((groupsi(i, 2)+2):(groupsi(i+1, 1)-2));
    end 

    topology.(group) = str2num(block);
end