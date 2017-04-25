%  writeMDArray(w, wfile, D)
% INPUTS:
%   w     = weighting
%   wfile = file name
%   D     = number of dimensions of w, pad with ones
function writeMDArray(w, wfile, D)
    %fprintf('writing w...');
    if nargin < 3
        D = ndims(w);
    end
    assert(D >= ndims(w));

    writeMDArrayDat(w, wfile);
    writeMDArrayHdr(w, wfile, D);
    
    fprintf('done\n');
end


function writeMDArrayHdr(w, wfile, D)
    % Write header
    fip_hdr = fopen([wfile, '.hdr'], 'w');
    for i = 1:D
        fprintf(fip_hdr, '%d ', size(w,i));
    end
    %fwrite(fip_hdr, mdsize(w, 1:D), 'long');
    fclose(fip_hdr);
end

function writeMDArrayDat(w, wfile)
    % Write data
    fip_dat = fopen([wfile, '.dat'], 'wb');
    fwrite(fip_dat, w(:), class(w));
    fclose(fip_dat);
end
