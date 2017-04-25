% array = readMDArray(file, type)
%   INPUTS:
%       file = filename, [file.dat, file.hdr] contain data/header
%       type =
%   OUTPUTS:
%       array = array data
function array = readMDArray(file, type)

    file_hdr = [file, '.hdr'];
    file_dat = [file, '.dat'];

    % header
    dims = load(file_hdr); %[];

    % file data
    fip = fopen(file_dat, 'rb');
    array = fread(fip, type);
    array = reshape(array, dims);
    fclose(fip);

end
