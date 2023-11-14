function outstring = clean_string(instring,newchar)
% Replace undesirable characters in a string
% Designed to make filenames avoid problematic characters
if (nargin() < 2), newchar = '_'; end
    
outstring = strrep(instring,'/', newchar);
outstring = strrep(outstring,'\',newchar);
outstring = strrep(outstring,' ',newchar);
outstring = strrep(outstring,'<',newchar);
outstring = strrep(outstring,'>',newchar);
outstring = strrep(outstring,':',newchar);
return