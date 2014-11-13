% pricingGame    Compute the social efficiency of the interference channel
%                game with pricing as a function of the (normalized)
%                pricing factor (as suggested at the end of Section
%                'Pricing the strategies'
% 
%                The system parameters can be changed to modify the nature
%                of the game
%
%
function pricingGame

close all; clc

fprintf('\n*** COMPUTING THE SOCIAL EFFICIENCY OF THE INTERFERENCE CHANNEL GAME WITH PRICING ***\n\n\n');


%% system parameters 
h=[0.75 0.25; 0.50 1.00]; %% channel power gains
Gamma=4; %% spreading gain
p=5*10^0; %% maximum power (all powers normalized to the AWGN power)
L=20; %% number of information data bits per packet

alpha=[0:0.01:0.10 0.1025:0.0025:0.13 0.131:0.001:0.15]; %% pricing factor (normalized to the AWGN power)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the NE with pricing (NEP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% since the NE is computed , using the best-response (BR) approach as a numerical 
%% search over the 2D set \mathcal{S}_1 x \mathcal{S}_2, where \mathcal{S}_k 
%% is continuous, for convenience we will discretize it over a finite number
%% of points powerGridPoints
powerGridPoints=10001; % the higher, the more accurate, but the slower
s1=linspace(0,p,powerGridPoints); % \mathcal{S}_1
s2=linspace(0,p,powerGridPoints); % \mathcal{S}_2

sumUtilities=zeros(1,length(alpha));

for index_alpha=1:length(alpha)
    
    fprintf('\n*** alpha = %.04f *** \n              ', alpha(index_alpha));

    %% computing player 2's BR to all powers s_1
    b2=zeros(1,length(s1)); %% b_2(s_1)
    for i=1:length(s1);
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(i/length(s1)/3)));
        mu2=computeMu(Gamma, h, s1(i), 2); %% m_2(s_1(i))
        u2=efficiencyFunction(mu2*s2,L)./s2 - alpha(index_alpha)*s2; %% \tilde{u}_2(s_1(i),s2) (vector of powerGridPoints points)
        [~,index_b2]=max(u2); %% index_b2=argmax(\tilde{u}_2)
        b2(i)=s2(index_b2); %% storing player 2's BR to s_1
    end

    %% computing player 1's BR to all powers s_2
    b1=zeros(1,length(s2)); %% b_1(s_2)
    for i=1:length(s2);
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(1/3+i/length(s2)/3)));
        mu1=computeMu(Gamma, h, s2(i), 1);  %% m_1(s_2(i))
        u1=efficiencyFunction(mu1*s1,L)./s1 - alpha(index_alpha)*s1; %% \tilde{u}_1(s_1,s2(i)) (vector of powerGridPoints points)
        [~,index_b1]=max(u1); %% index_b1=argmax(\tilde{u}_1)
        b1(i)=s1(index_b1);  %% storing player 1's BR to s_2
    end

    %% finding the fixed point of the BRs
    distance=zeros(1,length(s1)); %% distance across the two lines (BRs)
    for i=1:length(s1)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(2/3+i/length(s1)/3)));
        tmp=(s2-b2(i)).^2+(s1(i)-b1).^2; %% computing the distance of s2 wrt b2(i)
        distance(i)=min(tmp); %% selecting the minimum distance across the two lines (given b2(i))
    end
    [~,index]=min(distance); %% finding the point that minimizes the distance between BRs (i.e., the crossing pont)
    sNEP=[s1(index) b2(index)]; %% NEP profile \tilde{s}^*

    %% computing the performance of the NEP over the utility plan (considering the *original* utility)
    mu1=computeMu(Gamma, h, sNEP(2), 1); %% mu_1(s_2^*)
    uNEP(1)=efficiencyFunction(mu1*sNEP(1),L)/sNEP(1); % u_1(\tilde{s}^*)
    mu2=computeMu(Gamma, h, sNEP(1), 2); %% mu_2(s_1^*)
    uNEP(2)=efficiencyFunction(mu2*sNEP(2),L)/sNEP(2);  % u_2(\tilde{s}^*)
    
    %% storing the sum-utility as a function of alpha
    sumUtilities(index_alpha)=sum(uNEP);

end

%% plotting the results
plot(alpha, sumUtilities, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
hold on; grid on; box on;
title('Social efficiency of the continuous-power game with pricing');
xlabel('normalized pricing factor \alpha/\sigma^2');
ylabel('normalized sum-utility (u_1(s^*)+u_2(s^*))/(\sigma^2 t)');

%% plotting the best pricing factor alpha that maximizes the social efficiency of the NE
[maxSumUtilities, index]=max(sumUtilities);
plot(alpha(index)*[1 1], [min(sumUtilities) maxSumUtilities], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5);

fprintf('\n\n*** INTERFERENCE CHANNEL GAME WITH PRICING SOLVED! ***\n\n');
