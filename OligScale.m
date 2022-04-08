function f = OligScale(scales,N)

%f = 1 when m+n = 2 
%f = 0 when m+n = 21
switch scales
    case 'yes'
f = -1/19 * (N) + 21/19;

    case 'no'
f=1; 
end


end 