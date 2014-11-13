% utilityFunction    Plot the throughput and the utility function as
%                    functions of the SINR, using the models described in
%                    Section 'Introducing continuous powers'
% 
%                    The system parameters can be changed to modify the
%                    nature of the game and all the associated performance
%
%
function utilityFunction

L=20; % number of information bits per packet
sinrdB=-10:.1:20; % range of SINR values [dB]
sinr=10.^(0.1*sinrdB); % range of SINR values [linear]
throughput=(1-exp(-sinr)).^L; % throughput

[~,index]=max(throughput./sinr); % numerical search of the optimal SINR that maximizes the utility function
sinrStardB=10*log10(sinr(index)); % optimal SINR [dB]

hold on; grid on; box on;

figure(1);
[haxes,hline1,hline2] = plotyy(sinrdB, throughput, sinrdB, throughput./sinr,'plot','plot'); % plotting the two curves: throughput and utility function
ylabel(haxes(1),'normalized throughput') % label left y-axis
ylabel(haxes(2),'normalized utility') % label right y-axis
xlabel(haxes(2),'SINR [dB]') % label x-axis
%% setting the colors: throughput-> red, utility->blue
set(haxes(1),'YColor',[255 0 0]/255);
set(hline1, 'LineStyle', '-', 'Color', [255 0 0]/255, 'LineWidth', 2.0);
set(haxes(2),'YColor',[0 64 128]/255);
set(hline2, 'LineStyle', '-', 'Color', [0 64 128]/255, 'LineWidth', 2.0);
%% reporting the optimal SINR on the graph
plot(sinrStardB*[1 1], [0 0.89], 'Color', [0 64 128]/255, 'LineWidth', 2.0, 'LineStyle', '--')
text(sinrStardB*1.1,.05,'xS')
