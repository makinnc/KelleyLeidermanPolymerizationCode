function V = FiberVolume2(sol,R) 

L = 22.5;% * (sol.l(end) +1); % in nm

%R is in nm
% R=R/2;

V = pi .* (R).^2 .* L; % in nm^3



end 