% computeMu  Compute mu_k(s_{-k}) as defined in eq. (1)
%
% inputs:        Gamma (spreading factor)
%                h (2x2 matrix of channel power gains)
%                s (power used by player -k)
%                k (player index)
%
% outputs:       out (mu_k(s_{-k}))
%
function out=computeMu(Gamma, h, s, k)

if k==1
    out=Gamma*h(1,1)/(1+h(2,1)*s);
else
    out=Gamma*h(2,2)/(1+h(1,2)*s);
end