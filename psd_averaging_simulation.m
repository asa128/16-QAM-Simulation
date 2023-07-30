close all % close all figures
clear % clear workspace

numTrials = 1000; % number of trials for which to average PSD over
N = 50*16; % number of bits per trial (must be a multiple of 16)

A = 0.67; % peak amplitude of I or Q channel
sampPerSym = 64; % number of samples per symbol
k_T = 20; % 2*k_T represents the number of bit periods over which to compute the pulse

% carrier frequencies in Hz
f_c_1 = 6250;
f_c_2 = 8750;
f_c_3 = 11250;
f_c_4 = 13750;

B = 1250; % baseband bandwidth of each carrier in Hz
r = 0.5; % roll-off factor

D = 2*B/(r+1); % symbol rate
T_sym = 1/D; % symbol period in seconds
sampPerSecond = D * sampPerSym; % sample rate
t_step = 1/sampPerSecond; % time step
f_c = [f_c_1, f_c_2, f_c_3, f_c_4]; % array of carrier frequencies


% CHECK IF NUMBER OF BITS ARE A MULTIPLE OF 16
if mod(N, 16) ~= 0
    msgbox("Number of bits must be a multiple of 16.")
    return
end


% CHECK NYQUIST CRITERION
f_max = max(f_c) + B; % determine maximum frequency of signal
f_nyquist = sampPerSecond/2; % determine nyquist frequency
message = sprintf(strcat("WARNING: Your signal contains frequencies greater than the Nyquist frequency. Try decreasing the carrier frequency or increasing the number of samples per bit.\n\n", "Current Maximum Frequency: ", num2str(f_max), " Hz\n", "Current Maximum Carrier Frequency: ", num2str(max(f_c)), " Hz\n", "Current Nyquist Frequency: ", num2str(f_nyquist), " Hz\n"));
if f_max > f_nyquist % display message if maximum frequency of signal is above nyquist frequency
    msgbox(message)
end


% GENERATE PULSE
[h, t] = rrcpulse(k_T, T_sym, sampPerSym, r); % generate RRC pulse

t_s = -k_T*T_sym : t_step : (N/16-1)*T_sym+k_T*T_sym - t_step; % time samples for signal
t_b = 0:T_sym:(N*T_sym - T_sym); % time samples for bits

% GENERATE MULTIPLE SIGNALS AND AVERAGE THEIR PSDs
P_s = 0;
P_alt = 0;
for n=1:1:numTrials
    bits_tx = randbin(N, 0, 1); % generate vector a of N random bits

    % SEPARATE BITS, BY GROUPS OF 4, INTO EACH CARRIER
    bits_carrier_tx = zeros(4, N/4); % each row represents a different carrier
    for counter = 1:1:N/16
        bits_carrier_tx(1, counter*4-3:counter*4) = bits_tx(counter*16-15:counter*16-12); % bits for carrier 1
        bits_carrier_tx(2, counter*4-3:counter*4) = bits_tx(counter*16-11:counter*16-8); % bits for carrier 2
        bits_carrier_tx(3, counter*4-3:counter*4) = bits_tx(counter*16-7:counter*16-4); % bits for carrier 3
        bits_carrier_tx(4, counter*4-3:counter*4) = bits_tx(counter*16-3:counter*16); % bits for carrier 4
    end
    
    
    % GENERATE I-Q VALUES FOR EACH CARRIER
    % each row of the following arrays represents a different carrier
    carrier_i_tx = zeros(4, N/16);
    carrier_q_tx = zeros(4, N/16);
    carrier_bitstring_tx = strings(4, N/16); % strings of bit values (for placing on plots)
    
    for carrier = [1 2 3 4]
        for counter = 1:1:N/16
            tempArray = bits_carrier_tx(carrier, (counter*4-3) : (counter*4) );
            
            % Gray-coded constellation
            if tempArray == [0, 0, 0, 0]
                carrier_i_tx(carrier, counter) = -A;
                carrier_q_tx(carrier, counter) = -A;
                carrier_bitstring_tx(carrier, counter) = "0000";
            elseif tempArray == [0, 0, 0, 1]
                carrier_i_tx(carrier, counter) = -A;
                carrier_q_tx(carrier, counter) = -A/3;
                carrier_bitstring_tx(carrier, counter) = "0001";
            elseif tempArray == [0, 0, 1, 0]
                carrier_i_tx(carrier, counter) = -A;
                carrier_q_tx(carrier, counter) = A;
                carrier_bitstring_tx(carrier, counter) = "0010";
            elseif tempArray == [0, 0, 1, 1]
                carrier_i_tx(carrier, counter) = -A;
                carrier_q_tx(carrier, counter) = A/3;
                carrier_bitstring_tx(carrier, counter) = "0011";
            elseif tempArray == [0, 1, 0, 0]
                carrier_i_tx(carrier, counter) = -A/3;
                carrier_q_tx(carrier, counter) = -A;
                carrier_bitstring_tx(carrier, counter) = "0100";
            elseif tempArray == [0, 1, 0, 1]
                carrier_i_tx(carrier, counter) = -A/3;
                carrier_q_tx(carrier, counter) = -A/3;
                carrier_bitstring_tx(carrier, counter) = "0101";
            elseif tempArray == [0, 1, 1, 0]
                carrier_i_tx(carrier, counter) = -A/3;
                carrier_q_tx(carrier, counter) = A;
                carrier_bitstring_tx(carrier, counter) = "0110";
            elseif tempArray == [0, 1, 1, 1]
                carrier_i_tx(carrier, counter) = -A/3;
                carrier_q_tx(carrier, counter) = A/3;
                carrier_bitstring_tx(carrier, counter) = "0111";
            elseif tempArray == [1, 0, 0, 0]
                carrier_i_tx(carrier, counter) = A;
                carrier_q_tx(carrier, counter) = -A;
                carrier_bitstring_tx(carrier, counter) = "1000";
            elseif tempArray == [1, 0, 0, 1]
                carrier_i_tx(carrier, counter) = A;
                carrier_q_tx(carrier, counter) = -A/3;
                carrier_bitstring_tx(carrier, counter) = "1001";
            elseif tempArray == [1, 0, 1, 0]
                carrier_i_tx(carrier, counter) = A;
                carrier_q_tx(carrier, counter) = A;
                carrier_bitstring_tx(carrier, counter) = "1010";
            elseif tempArray == [1, 0, 1, 1]
                carrier_i_tx(carrier, counter) = A;
                carrier_q_tx(carrier, counter) = A/3;
                carrier_bitstring_tx(carrier, counter) = "1011";
            elseif tempArray == [1, 1, 0, 0]
                carrier_i_tx(carrier, counter) = A/3;
                carrier_q_tx(carrier, counter) = -A;
                carrier_bitstring_tx(carrier, counter) = "1100";
            elseif tempArray == [1, 1, 0, 1]
                carrier_i_tx(carrier, counter) = A/3;
                carrier_q_tx(carrier, counter) = -A/3;
                carrier_bitstring_tx(carrier, counter) = "1101";
            elseif tempArray == [1, 1, 1, 0]
                carrier_i_tx(carrier, counter) = A/3;
                carrier_q_tx(carrier, counter) = A;
                carrier_bitstring_tx(carrier, counter) = "1110";
            elseif tempArray == [1, 1, 1, 1]
                carrier_i_tx(carrier, counter) = A/3;
                carrier_q_tx(carrier, counter) = A/3;
                carrier_bitstring_tx(carrier, counter) = "1111";
            end
        end
    end
    
    
    % SYNTHESIZE BASEBAND I-Q SIGNALS BY SUMMING TIME SHIFTED RRC PULSES
    % each row of the following arrays represents a different carrier
    carrier_i_tx_bb = zeros(4, length(t_s));
    carrier_q_tx_bb = zeros(4, length(t_s));
    
    for carrier = [1 2 3 4]
        for k = 0:1:N/16-1
            for i = 1:1:(2*k_T*sampPerSym)
                carrier_i_tx_bb(carrier, i+k*sampPerSym) = carrier_i_tx_bb(carrier, i+k*sampPerSym) + carrier_i_tx(carrier, k+1)*h(i);
                carrier_q_tx_bb(carrier, i+k*sampPerSym) = carrier_q_tx_bb(carrier, i+k*sampPerSym) + carrier_q_tx(carrier, k+1)*h(i);
            end
        end
    end
    
    
    % GENERATE MODULATED SIGNALS
    s_mod = zeros(4, length(t_s)); % each row represents a different carrier
    
    for carrier = [1 2 3 4] % modulate each carrier separately
        s_mod(carrier,:) = carrier_i_tx_bb(carrier,:).*cos(2.*pi.*f_c(carrier).*t_s) - carrier_q_tx_bb(carrier,:).*sin(2.*pi.*f_c(carrier).*t_s);
    end
    
    s = sum(s_mod, 1); % sum the modulated carriers to form the signal
    
    % DETERMINE ENERGY IN THE SIGNAL
    E = 0;
    for k = 1:1:length(s) % calcuate total energy in the signal using Riemann sum
        E = E + (s(k))^2*(t_step);
    end

    P_alt_i = E/ (N/16*T_sym); % alternative power calculation method, dividing energy by the data time

    % SIMULATION OF THE PSD OF THE MODULATED SIGNAL
    S = t_step .* fftshift( fft(s) ); % calculate DFT of s, center it around 0, scale it appropriately
    P_s_i = (abs(S)).^2 / (N/16 * T_sym); % compute approximate PSD of s
    f_step = 1 / (length(S) * t_step); % frequency step for approximating CTFT from DFT
    f = -f_step*(length(S)/2) : f_step : f_step*(length(S)/2-1); % calculate frequencies corresponding to the spectrum, S
    
    P_s = P_s + (1./numTrials).*P_s_i; % add PSD instance to sum for averaging
    
    P_alt = P_alt + (1./numTrials).*P_alt_i; % add alternative power calculation to sum for averaging

    disp(strcat("Progress = ", num2str(n/numTrials*100), "%")) % display progress in loop
end

P_s_dB = 10.*log10(P_s);

figure(3)

subplot(2,1,1)
plot(f, P_s,'linewidth',2) % plot PSD in absolute units
xlim([-1.5.*max(f_c) 1.5.*max(f_c)])
title('{\boldmath $P_{s}(f)$}', 'Interpreter', 'latex','FontSize',24)
xlabel('{\boldmath $f$} (Hz)', 'Interpreter', 'latex','FontSize',18)
ylabel('{\boldmath $P_{s}(f)$} (W/Hz)', 'Interpreter', 'latex','FontSize',18)
hold off

subplot(2,1,2)
plot(f, P_s_dB,'linewidth',2) % plot PSD in decibel units
xlim([-1.5.*max(f_c) 1.5.*max(f_c)])
title('{\boldmath $P_{s}(f)$}', 'Interpreter', 'latex','FontSize',24)
xlabel('{\boldmath $f$} (Hz)', 'Interpreter', 'latex','FontSize',18)
ylabel('{\boldmath $P_{s}(f)$} (dBW/Hz)', 'Interpreter', 'latex','FontSize',18)
hold off


% DETERMINE SIGNAL POWER
P = 0;
for k = 1:1:length(P_s) % calcuate total power using Riemann sum
    P = P + P_s(k) * f_step;
end

P_inband = 0;
for k = 1:1:length(P_s) % calcuate in-band power using Riemann sum
    if (f(k) >= 5000 && f(k) <= 15000)
        P_inband = P_inband + 2*P_s(k) * f_step;
    end
end

% display total and in-band power values
disp(strcat("Total signal power (W) = ", num2str(P)))
disp(strcat("Normalized in-band signal power (W) = ", num2str(P_inband)))
disp(strcat("Alternative power calculation method (W) = ", num2str(P_alt)))
