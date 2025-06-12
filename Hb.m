function out = Hb(x)
    if x < 0
        out = 0;
    elseif x > 1
        out = x-0.5;
    else
        out = 0.5*x*x;
    end
end