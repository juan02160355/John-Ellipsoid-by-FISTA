function out = phiprim(z)
    if z<0
        out = 0;
    elseif z>1
        out = 1;
    else
        out = z;
    end
end