close all % close all figures
clear % clear workspace

n = 1; % index value to be incremented for plot legend array
legendArray = cell(4,1); % array for plot legend

k_T = 5; % 2*k_T represents the number of bit periods over which to compute the pulse
R = 1000; % bit rate in bits per second
T_b = 1/R; % bit period in seconds
sampPerBit = 64; % number of samples per bit

for r = 0.25:0.25:1 % loop over possible roll off factors, r

    [h, t] = rrcpulse(k_T, T_b, sampPerBit, r); % call RRC pulse generator function with specified parameters
    
    % h is a vector of RRC function values
    % t is a vector of times at which the corresponding values in h occur
    
    % plot the RRC function that was produced
    figure(1)
    plot(t,h,'linewidth',2)
    hold on % plot subsequent pulses on the same pair of axes
    
    legendArray(n) = { strcat('r = ', num2str(r)) }; % append r value to legend
    n = n+1; % increment index n
end

grid on
title('Root Raised Cosine Rolloff Pulse', 'Interpreter', 'latex','FontSize',24)
xlabel('{\boldmath $t$}', 'Interpreter', 'latex','FontSize',18)
ylabel('{\boldmath $h(t)$}', 'Interpreter', 'latex','FontSize',18)
legend(legendArray) % add legend to the plot