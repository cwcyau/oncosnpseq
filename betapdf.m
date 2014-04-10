% OncoSNP
% 
% Created by Christopher Yau, Department of Statistics, University of Oxford
%
% (c) Copyright 2010 University of Oxford. All rights reserved.
%
% HE INSTALLATION AND USE OF THIS SOFTWARE AND SOURCE CODE CONSTITUTES THE ACCEPTANCE OF THE LICENCE AGREEMENT. A COPY OF THE LICENCE IS CONTAINED IN THE FILE licence.txt THATACCOMPANIES THIS PACKAGE
%
function y = betapdf(x,a,b)

tmp = (a - 1).*log(x) + (b - 1).*log((1 - x)) - betaln(a,b);
y = exp(tmp);


	

