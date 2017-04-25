% yz = zpad(y, dimsNew);
% 
% zero padding 
% dimsNew = new zero pad dimensions
function yNew = zpad(y, dimsNew)
    for d = 1:length(dimsNew)
        dims(d) = size(y,d);
    end
    if ndims(y) > length(dimsNew)
        dims
        dimsNew
        error('Cannot pad dims (first) and dimsNew (second)');
    end

    yNew = zeros(dimsNew(:)');

    for k = 1:length(dims)
        middle{k} = ((-dims(k)/2 + 1):dims(k)/2) + dimsNew(k)/2; 
    end

    comm = 'yNew(middle{1}';
    for i = 2:length(dimsNew)
        comm = [comm, sprintf(', middle{%d}', i)];
    end
    comm = [comm, ') = y;'];
    eval(comm);
end
