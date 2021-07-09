% [pat, cost] = best_candidate(bcfast_path, w, totSamples, ... )
%
% Best candidate sampling, default uses exact method, but 
% options specify approximate method.
%
% INPUTS:
%   bcfast_path = path to bcfast executable
%   w           = weighting function
%   totSamples  = total samples for all frames
% Optional INPUTS:
%   'Max Per Frame' = set a max # samples per frame
%   'T', T          = threshold for w
%   'M', M          = number of neighbors for w
% OUTPUTS:
%   pat     = pattern
function [pat] = best_candidate(bcfast_path, w, totSamples, varargin)

nt = size(w,3);

% Options
opt.exact = true;
opt.totSamples = totSamples;
i = 1;
while i <= length(varargin)
    if strcmp(varargin{i}, 'Max Per Frame')
        opt.maxPerFrame = varargin{i+1};
        i = i + 2;
    elseif strcmp(varargin{i}, 'T')
        opt.T = varargin{i+1};
        opt.exact = false;
        i = i + 2;
    elseif strcmp(varargin{i}, 'K')
        opt.K = varargin{i+1};
        opt.exact = false;
        i = i + 2;
    else
        error('option not recognized');
    end
end
if ~isfield(opt, 'K')
    opt.K = 0;
end
if ~isfield(opt, 'T')
    opt.T = 0;
end
if ~isfield(opt, 'maxPerFrame')
    opt.maxPerFrame = totSamples;
end
opt.v2 = opt.T || opt.K;

% Checks
assert(all(isreal(w(:))));
assert(ndims(w) <= 5);

% Input file
wfile = [tempname, '_egl'];
writeMDArray(w, wfile, 5);

% Output file
patfile = [tempname, '_egl.txt'];

% Argument list
args = ['--w ', wfile, ' --pat ' patfile, ' --max-per-frame ', num2str(opt.maxPerFrame), ' --S ' num2str(opt.totSamples), ' '];
if opt.exact
    args = [args, '--e 1 '];
else
    args = [args, '--e 0 '];
end
if opt.T
    args = [args, '--T ', num2str(opt.T), ' '];
end
if opt.K
    args = [args, '--K ', num2str(opt.K), ' '];
end

% Run the sampling algorithm
run_sysstr([bcfast_path, ' ', args]);

% Read pattern
pat = read_pat(patfile);

% Clean up input and output
delete(patfile);
delete([wfile, '.dat']);
delete([wfile, '.hdr']);

end % end function

function run_sysstr(comm)
    fprintf('executing: %s\n', comm);
    system(comm);
end

function pat = read_pat(patfile)
    %fprintf('reading pattern file ');
        
    fip = fopen(patfile, 'r');

    line = fgets(fip);
    %fprintf('dims: %s ', line);
    dims = str2num(line);
    pat = zeros(dims);
    line = fgets(fip);
    while ischar(line)
        k = str2num(line) + 1;
        %fprintf('%s -> %d\n', line, k);
        pat(k) = pat(k) + 1;
        line = fgets(fip);
    end

    fclose(fip);
    %fprintf('done\n');
end

