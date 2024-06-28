function outvector = fourpammapA(invector)
    l = length(invector);

    outvector = zeros(1,(l/2));

    for k = 1:2:l
        ind = (k+1)/2;
        if invector(k) == '1'
            outvector(ind) = str2double(invector(k+1)) + 2;
            continue;
        end
        outvector(ind) = str2double(invector(k+1));
    end
    
end