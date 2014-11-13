% findGammaStar  Compute the max and the argmax of the normalized utility
%                function, defined as u=f(x)/x, where f(x) is computed 
%                using efficiencyFunction
%
% inputs:        L (number of information bits per packet
%
% outputs:       uStar (maximum value of the normalized utility)
%                gammaStar (argmax of the normalized utility)
%
function [uStar, gammaStar]=findGammaStar(L)
x=linspace(0, L, 10000); % support of the utility function
u=efficiencyFunction(x, L)./x; % utility function
[uStar,imax]=max(u); % computing max and argmax
tmp=x(imax); % temporary gammaStar

isFinished=0; % loop control initialization
while (isFinished==0), % gammaStar neds to be refined
    x=linspace(tmp*0.8,tmp*1.2, 10000); % refining the support of the utility function
    u=efficiencyFunction(x, L)./x; % utility function
    [uStar,imax]=max(u); % computing max and argmax
    gammaStar=x(imax); % assigning gammaStar the current argmax value
    if (abs(tmp-gammaStar)/tmp<1e-3) % checking the accuracy of the argmax solution
        isFinished=1; % exiting the loop
    end
end