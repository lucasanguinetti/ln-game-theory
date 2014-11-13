% findGamma      Compute the SINR x such that in=f(x)/x, where f(x) is computed 
%                using efficiencyFunction, for a given in
%
% inputs:        in (desired value for the utility )
%                L (number of information bits per packet)
%                gammaStar (optimal SINR that maximizes f(x)/x)
%
% outputs:       out (SINR x such that in=f(x)/x)
%
function out=findGamma(in, L, gammaStar)
threshold=1e-2;
x=linspace(0,gammaStar, 1001); % support of the utility function
f=efficiencyFunction(x, L)./x-in;
[fmin,imin]=min(abs(f));
out=x(imin);

while (fmin>threshold)
    x=linspace(out*0.8,out*1.2, 1001);
    f=efficiencyFunction(x, L)./x-in;
    [fmin,imin]=min(abs(f));
    out=x(imin);
end