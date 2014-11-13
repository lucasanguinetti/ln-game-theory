% plotUtilityPlan    Plot the utility plan for the interference channel
%                    game (Pareto frontier, performance @ Nash equilibrium,
%                    social optimum, Nash equilibrium of the priced game,
%                    Nash equilibrium of the repeated game, and Nash
%                    bargaining solution)
% 
%                    The system parameters can be changed to modify the
%                    nature of the game and all the associated performance
%
%
function plotUtilityPlan

close all; clc

fprintf('\n*** SOLVING THE INTERFERENCE CHANNEL GAME ***\n\n\n');


figure;
hold on; grid on; box on;

%% dummy points to correctly display the graph legend
plot(-1, -1, 'Marker', 'd', 'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 6.0, 'LineStyle', 'none');
plot(-1, -1, 'Marker', 's', 'MarkerEdgeColor', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 6.0, 'LineStyle', 'none');
plot(-1, -1, 'Marker', 'x', 'MarkerEdgeColor', [0 1 1], 'MarkerFaceColor', [0 1 1], 'MarkerSize', 6.0, 'LineStyle', 'none');
plot(-1, -1, 'Marker', '*', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 6.0, 'LineStyle', 'none');
plot(-1, -1, 'Marker', 'o', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 6.0, 'LineStyle', 'none');


%% system parameters 
h=[0.75 0.25; 0.50 1.00]; %% channel power gains
Gamma=4; %% spreading gain
p=5*10^0; %% maximum power (all powers normalized to the AWGN power)
L=20; %% number of information data bits per packet
[uStar, gammaStar]=findGammaStar(L); %% max and argmax of the normalized utility


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the Pareto frontier (PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncomputing the Pareto frontier (PF)...\n\n              ');

FP_points = 1001; % number of points of the Pareto frontier (the higher, the more accurate, but the slower)

% setting up the upper bound for player 1'2 utility (when s_2=0)
mu1=computeMu(Gamma, h, 0, 1); % mu_1(0) -- see eq. (1)
s1=gammaStar/mu1; % s_1^*(s_2=0)
u1max=efficiencyFunction(gammaStar,L)/s1; % player 1's maximum utility

% computing the vector for the x axis of the utility plan (linearly spaced)
u1PO=linspace(0, u1max, FP_points); 

% preparing the vector for the y axis of the utility plan (empty
% initialization)
u2PO=zeros(1, FP_points);

%% since the PF is computed as a numerical search over the 2D set \mathcal{S}_1 x \mathcal{S}_2,
%% where \mathcal{S}_k is continuous, for convenience we will discretize it over a finite number
%% of points powerGridPoints
powerGridPoints=1001; % as before, the higher, the more accurate, but the slower

maxIter=10; % maximum number of iterations (to prevent numerical problems due to finite arithmetics in the numerical search)
threshold=1e-2; % accuracy threshold for computing u2PO (the lower, the more accurate, but the slower)

for index_u1=FP_points-1:-1:1 % computing max(u2) as a function of u1PO(index_u1) ( u2PO(FP_points)=0 by definition (s_2=0), and hence starting from the (N-1)th point) )
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(1-index_u1/FP_points)));
    
    isFinished=0; % initialization for the loop control 
    iter=0; % initialization of the number of iterations
    oldU2PO=0; %% initialization for u_2
    
    %% initializing the subset of \mathcal{S}_2 that provides max(u2) as a function of u1PO(index_u1) (to be further refined)
    sLeft=eps; %% left margin
    sRight=p; %% right margin

    while (and(isFinished==0,iter<maxIter)), %% the iterative numerical search is on
        
        iter=iter+1; % updating the number of iterations
        s2=linspace(sLeft,sRight,powerGridPoints); % setting up the subset of \mathcal{S}_2 (linearly spaced)
        u2=zeros(1,length(s2)); % setting up the vector containing u_2([s_1,s_2]) (same length as s_2)

        for index_s2=1:length(s2) %% scanning the vector s_2
            
            mu1=computeMu(Gamma, h, s2(index_s2), 1); %% computing mu_1(s_2)
            u=u1PO(index_u1)/mu1; %% computing the normalized utility when considering u1PO(index_u1) and s_2(index_s2) (that modifies mu_1(s_2))
            
            if (u<=uStar) %% if u>uStar, there exists no gamma1 such that u1=u1PO(index_u1) -> break (since u2(index_s2) is already set to 0, it will be discarded when computing the maximum)
                
                gamma1=findGamma(u1PO(index_u1)/mu1, L, gammaStar); %% computing gamma_1 such that u_1([s_1,s_2])=u1PO(index_u1)
                s1=gamma1/mu1; %% computing s_1 such that u_1([s_1,s_2])=u1PO(index_u1)
                
                if (s1<=p) %% if s_1>p, there exists no s1 such that u1=u1PO(index_u1) -> break (since u2(index_s2) is already set to 0, it will be discarded when computing the maximum)
                    
                    mu2=computeMu(Gamma, h, s1, 2); %% computing mu_2(s_1), where s_1 is a function of both u1PO(index_u1) and s2
                    u2(index_s2)=efficiencyFunction(mu2*s2(index_s2),L)/s2(index_s2); %% computing u2([s_1,s_2]), where s_1 is a function of both u1PO(index_u1) and s2
                    
                else
                    
                    break;
                    
                end
                
            else
                
                break;
   
            end
        end
        
        [newU2PO, imax]=max(u2); % computing the maximum achievable utility (and its argamx) given u1PO(index_u1) and the current subset [sLeft, sRight]
        
        if (abs(newU2PO-oldU2PO)/oldU2PO>threshold) % accuracy check (if beyond the threshold, there is room for refining the solution by shrinking the subset [sLeft, sRight]

            sLeft=max(eps,s2(imax)*0.8); %% new left margin (across the current argmaximizer)
            sRight=min(p,s2(imax)*1.2); %% new right margin (across the current argmaximizer)
            oldU2PO=newU2PO; %% saving the current value for u2PO
        
        else
            
            isFinished=1; %% exit the loop, as the solution found is accurate enough
        
        end
        
    end
    
    if (iter<maxIter) %% the solution found is accurate enough -> store the result
        u2PO(index_u1)=newU2PO;
    end %% if the numerical solution is not converging, due to finite arithmetics -> assign 0 (already pre-assigned to u2PO(index_u1))

end

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b 100%% completed\n');

%% plotting the PF
area_plot=area(u1PO, u2PO);
set(area_plot(1),'FaceColor',0.9*[1 1 1]) %% filling the utility plan with shaded background

clear s1; clear s2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the Nash equilibrium (NE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncomputing the Nash equilibrium (NE)...\n\n              ');

%% since the NE is computed , using the best-response (BR) approach as a numerical 
%% search over the 2D set \mathcal{S}_1 x \mathcal{S}_2, where \mathcal{S}_k 
%% is continuous, for convenience we will discretize it over a finite number
%% of points powerGridPoints
powerGridPoints=10001; % as before, the higher, the more accurate, but the slower
s1=linspace(0,p,powerGridPoints); % \mathcal{S}_1
s2=linspace(0,p,powerGridPoints); % \mathcal{S}_2

%% computing player 2's BR to all powers s_1
b2=zeros(1,length(s1)); %% b_2(s_1)
for i=1:length(s1);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(i/length(s1)/3)));
    mu2=computeMu(Gamma, h, s1(i), 2); %% m_2(s_1(i))
    u2=efficiencyFunction(mu2*s2,L)./s2; %% u_2(s_1(i),s2) (vector of powerGridPoints points)
    [~,index_b2]=max(u2); %% index_b2=argmax(u_2)
    b2(i)=s2(index_b2); %% storing player 2's BR to s_1
end

%% computing player 1's BR to all powers s_2
b1=zeros(1,length(s2)); %% b_1(s_2)
for i=1:length(s2);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(1/3+i/length(s2)/3)));
    mu1=computeMu(Gamma, h, s2(i), 1);  %% m_1(s_2(i))
    u1=efficiencyFunction(mu1*s1,L)./s1; %% u_1(s_1,s2(i)) (vector of powerGridPoints points)
    [~,index_b1]=max(u1); %% index_b1=argmax(u_1)
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
sNE=[s1(index) b2(index)]; %% NE profile

%% computing the performance of the NE over the utility plan
mu1=computeMu(Gamma, h, sNE(2), 1); %% mu_1(s_2^*)
uNE(1)=efficiencyFunction(mu1*sNE(1),L)/sNE(1); % u_1(s^*)
mu2=computeMu(Gamma, h, sNE(1), 2); %% mu_2(s_1^*)
uNE(2)=efficiencyFunction(mu2*sNE(2),L)/sNE(2);  % u_2(s^*)
% plotting the result
plot(uNE(1), uNE(2), 'Marker', 'd', 'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 6.0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the social welfare (SW)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\ncomputing the social welfare (SW)...\n\n');

sumUtilities=u1PO+u2PO; %% since the SO point belongs to the PF, we can just consider the sum-utility over the PF
[SW, index_SO]=max(sumUtilities); %% the SO point is by definition the one that maximizes the sum-utility
% plotting the result
plot(u1PO(index_SO), u2PO(index_SO), 'Marker', 's', 'MarkerEdgeColor', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 10.0);
% plotting the tangent line with slope -1
plot([0 SW], [SW 0], 'Color', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 0.5);

% plotting the performance of the NE for the repeated game (which coincides
% with the SW)
plot(u1PO(index_SO), u2PO(index_SO), 'Marker', '*', 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 6.0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the NE with pricing (NEP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ncomputing the NE with pricing...\n\n              ');

%% since the NE is computed , using the best-response (BR) approach as a numerical 
%% search over the 2D set \mathcal{S}_1 x \mathcal{S}_2, where \mathcal{S}_k 
%% is continuous, for convenience we will discretize it over a finite number
%% of points powerGridPoints
powerGridPoints=10001; % as before, the higher, the more accurate, but the slower
s1=linspace(0,p,powerGridPoints); % \mathcal{S}_1
s2=linspace(0,p,powerGridPoints); % \mathcal{S}_2

alpha=0.12; %% pricing factor (normalized to the AWGN power)

%% computing player 2's BR to all powers s_1
b2=zeros(1,length(s1)); %% b_2(s_1)
for i=1:length(s1);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(i/length(s1)/3)));
    mu2=computeMu(Gamma, h, s1(i), 2); %% m_2(s_1(i))
    u2=efficiencyFunction(mu2*s2,L)./s2 - alpha*s2; %% \tilde{u}_2(s_1(i),s2) (vector of powerGridPoints points)
    [~,index_b2]=max(u2); %% index_b2=argmax(\tilde{u}_2)
    b2(i)=s2(index_b2); %% storing player 2's BR to s_1
end

%% computing player 1's BR to all powers s_2
b1=zeros(1,length(s2)); %% b_1(s_2)
for i=1:length(s2);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b %03d%% completed', ceil(100*(1/3+i/length(s2)/3)));
    mu1=computeMu(Gamma, h, s2(i), 1);  %% m_1(s_2(i))
    u1=efficiencyFunction(mu1*s1,L)./s1 - alpha*s1; %% \tilde{u}_1(s_1,s2(i)) (vector of powerGridPoints points)
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
% plotting the result
plot(uNEP(1), uNEP(2), 'Marker', 'x', 'MarkerEdgeColor', [0 1 1], 'MarkerSize', 8.0);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computing the Nash bargaining solution (NBS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\ncomputing the Nash bargaining solution (NBS)...\n\n\n');

%% identifying the PO points that lie in the region of Pareto improvement
indexLeft=FP_points+1-sum(u1PO>=uNE(1)); % left margin
indexRight=sum(u2PO>=uNE(2)); % right margin
nbsUtilities=(u1PO(indexLeft:indexRight)-uNE(1)).*(u2PO(indexLeft:indexRight)-uNE(2)); %% computing (u_1-u_1(s^*))*(u_2-u_2(s^*)) for the PO points lying in the region of Pareto improvement
[NBS, index_NBS]=max(nbsUtilities); % computing the one that corresponds to the NBS (by definition, it is the argmaximizer of the above product)
% plotting the result
plot(u1PO(indexLeft-1+index_NBS), u2PO(indexLeft-1+index_NBS), 'Marker', 'o', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 6.0);
% plotting the hyperbola
ux=linspace(uNE(1)+eps, 0.8, 1000); %% support for the hyperbola
uy=NBS./(ux-uNE(1))+uNE(2); %% equation of a hyperbola with x=uNE(1) and y=uNE(2) as the asymptotes
plot(ux, uy, 'Color', [0 0 1], 'LineStyle', '--', 'LineWidth', 0.5);

%% enriching the picture
axis([0.0 ceil(max(u1PO(end),u2PO(1))) 0.0 ceil(max(u1PO(end),u2PO(1)))]);
text(u1PO(end)/5, u2PO(1)/10, 'PF');
arrow([u1PO(end)/3.5 u2PO(1)/10], [u1PO(ceil(FP_points*3.5/5)) u2PO(ceil(FP_points*3.5/5))], 12);
legend('NE', 'SO', 'NE (w/ pricing)', 'NE (repeated game)', 'NBS', 'Location', 'SouthEast');
title('Normalized utility plan (continuous-power IC game)');
xlabel('normalized utility u_1(s) \cdot \sigma^2/t');
ylabel('normalized utility u_2(s) \cdot \sigma^2/t');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% displaying the inset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% computing the boundaries of the Pareto-improvement region
[~,index_left]=min(abs(u1PO-uNE(1)));
[~,index_right]=min(abs(u2PO-uNE(2)));

% specifying position and size of the rectangle
r_x=0.5*(u1PO(index_right)+u1PO(index_left));
r_w=1.25*(u1PO(index_right)-u1PO(index_left));
r_y=0.5*(u2PO(index_left)+u2PO(index_right));
r_h=1.25*(u2PO(index_left)-u2PO(index_right));
rectangle('Position', [r_x-r_w/2, r_y-r_h/2, r_w, r_h], 'EdgeColor', [0, 0, 0], 'LineWidth', 1, 'LineStyle', ':'); %% placing the box on the main figure

% specifying position and size of the inset
a_x = 0.55; a_y = 0.5; a_w = 0.3; a_h = a_w*r_h/r_w;
axes('Units', 'Normalized', 'Position', [a_x, a_y, a_w, a_h], 'XTick', [], 'YTick', [], 'Box', 'on', 'LineWidth', 1, 'Color', [255, 251, 237]/255);
hold on;

% plotting the Pareto-improvement region
X=u1PO(index_left:index_right);
Y=[u2PO(index_left:index_right), uNE(2)*ones(1,length(X))];
X=[X,fliplr(X)];
fill(X,Y,0.95*[1 1 1]);

% plotting the PF
plot(u1PO, u2PO, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);

% plotting the projection of the NE over the PF
plot([uNE(1) max(u1PO(end), u2PO(1))*1.2],[uNE(2) max(u1PO(end), u2PO(1))*1.2+(uNE(2)-uNE(1))], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5);

% plotting the five relevant points
plot(uNE(1), uNE(2), 'Marker', 'd', 'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'MarkerSize', 6.0);
plot(u1PO(index_SO), u2PO(index_SO), 'Marker', 's', 'MarkerEdgeColor', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 10.0);
plot(uNEP(1), uNEP(2), 'Marker', 'x', 'MarkerEdgeColor', [0 1 1], 'MarkerSize', 8.0);
plot(u1PO(index_SO), u2PO(index_SO), 'Marker', '*', 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 6.0);
plot(u1PO(indexLeft-1+index_NBS), u2PO(indexLeft-1+index_NBS), 'Marker', 'o', 'MarkerEdgeColor', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerSize', 6.0);

text(r_x, r_y-r_h*0.3, 'PI');
axis([r_x-r_w/2, r_x+r_w/2, r_y-r_h/2, r_y+r_h/2]);

fprintf('\n*** INTERFERENCE CHANNEL GAME SOLVED! ***\n\n');
