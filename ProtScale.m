function f = ProtScale(scales, N)

%f = 1 when N = 1 
%f = 0 when N = 11
switch scales
    case 'yes'
        
f = -1/10 * (N) + 11/10; 

    case 'no'
        
f=1; 

end

end