function timestamps=opentfile(tfile)
% pass tfile as string. Returns timestamps

if nargin~=1 %if file is not specified
[fileName,folderName] = uigetfile({'*.t','All t Files';...
          '*.*','All Files' },'Select file to open',...
          'C:\Data\');
tfile=[folderName fileName];
end

fid=fopen(tfile);
timestamps=fread(fid, 'uint32=>double','b'); 
fclose(fid);
% precision: source=>output. Remove if data is already in double 
% (may have been converted to uint32 when writing to file)

%then convert to whatever time scale 
% e.g. if data is at sampling frequency of 10kHz and want to convert to
% seconds: timestamps=timestamps/10000;