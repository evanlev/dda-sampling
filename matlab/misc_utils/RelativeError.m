% Evan Levine
% 
% err = RelativeError(xref, xn)
%     returns || xref - xn ||2^2/||xref||2^2 (dB)
function err = RelativeError(mref, xn)
    normmref = norm(mref(:),2);
    if normmref == 0
        error('Reference image has 0 norm');
    end
    err = 20*log10(norm(mref(:)-xn(:),2)/normmref);
end
