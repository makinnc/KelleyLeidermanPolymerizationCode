function [Pgg, Pgp, Pgpgp]= ProbGamPfun(fracgp)
Pgg = (1-fracgp)^2;
Pgp = 2 * (1-fracgp)*fracgp; 
Pgpgp =  fracgp^2; 



end