function out = MIP(A,B)
    C = A.*B;
    D = sum(C);
    out = sum(D);
end