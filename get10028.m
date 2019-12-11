% This script is used to download 5 micrographs of the EMPIAR 10028
% dataset. See the Example section on the user manual for more details.

mkdir ./example/micrographs
Nfiles=5;

ftpServer = 'ftp.ebi.ac.uk';
fprintf('Establishing an FTP connection with the EMPIAR server ...\n');
ngdc = ftp(ftpServer);
fprintf('FTP connection was established.\n');

ftpAddress = '/empiar/world_availability/10028/data/Micrographs/Micrographs_part1';
fprintf('Change direcotry to %s\n',ftpAddress);
cd(ngdc,ftpAddress);

for k=1:Nfiles
    fname=sprintf('%03d.mrc',k);    
    fprintf('Downloading %s\n',fname);
    micrograph =mget(ngdc, fname,'./example/micrographs/');
end

fprintf('Closing the FTP connection...\n');
close(ngdc);
fprintf('The FTP connection was closed.\n');
