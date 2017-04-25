% res = repmat2(data, dims, repvec)
% 
% Repeat along dims repvec times
% e.g.   sns_maps = repmat2(sns_maps, [TE_DIM, TIME_DIM], [nechoes, nt])
function res = repmat2(data, dims, repvec)
    
    repvec2 = ones(1,max(max(dims(:)),ndims(data)));
    repvec2(dims) = repvec;
    res = repmat(data, repvec2);

end

