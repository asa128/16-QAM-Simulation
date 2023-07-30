close all % close all figures
clear % clear workspace

noiseVar = 0; % noise variance
bits_tx = [1 0 0 1 1 0 1 0 0 1 1 0 1 1 0 1 0 1 1 1 0 1 0 0 1 0 0 0 1 1 1 0]; % bits to transmit (number must be a multiple of 16)
N = length(bits_tx);
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


% PLOT BASEBAND SIGNALS AND BITS
figure(1)
for carrier = [1 2 3 4]
    subplot(2,2,carrier) % subplot for each carrier
    plot(t_s, carrier_i_tx_bb(carrier,:),'linewidth',2) % I-channel
    hold on
    plot(t_s, carrier_q_tx_bb(carrier,:),'linewidth',2) % Q-channel
    grid on
    xlim([-5*T_sym (N/16+5)*T_sym])
    ylim([-2*A 2*A])
    title(strcat("Baseband Input - Carrier ", num2str(carrier)), 'Interpreter', 'latex','FontSize',20)
    xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
    ystring = strcat("{\boldmath $b_{", num2str(carrier), "}(t)$}");
    ylabel(ystring, 'Interpreter', 'latex','FontSize',16)
    legend('I-channel', 'Q-channel')
    for cnt = 1:1:N/16
        text((cnt-1)*T_sym, 1.5*A, carrier_bitstring_tx(carrier,cnt), 'HorizontalAlignment', 'center')
    end
end


% GENERATE MODULATED SIGNALS
s_mod = zeros(4, length(t_s)); % each row represents a different carrier

for carrier = [1 2 3 4] % modulate each carrier separately
    s_mod(carrier,:) = carrier_i_tx_bb(carrier,:).*cos(2.*pi.*f_c(carrier).*t_s) - carrier_q_tx_bb(carrier,:).*sin(2.*pi.*f_c(carrier).*t_s);
end

% PLOT EACH MODULATED CARRIER
figure(2)
for carrier = [1 2 3 4]
    subplot(2,2,carrier) % subplot for each carrier
    plot(t_s, s_mod(carrier,:),'linewidth',2) % I-channel
    hold on
    grid on
    xlim([-5*T_sym (N/16+5)*T_sym])
    ylim([-2.5*A 2.5*A])
    title(strcat("Modulated Carrier ", num2str(carrier)), 'Interpreter', 'latex','FontSize',20)
    xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
    ystring = strcat("{\boldmath $s_{", num2str(carrier), "}(t)$}");
    ylabel(ystring, 'Interpreter', 'latex','FontSize',16)
    for cnt = 1:1:N/16
        text((cnt-1)*T_sym, 2*A, carrier_bitstring_tx(carrier,cnt), 'HorizontalAlignment', 'center')
    end
end

s = sum(s_mod, 1); % sum the modulated carriers to form the signal


% PLOT TRANSMITTED SIGNAL
figure(3)
subplot(2,1,1)
plot(t_s, s,'linewidth',2)
grid on
xlim([-5*T_sym (N/16+5)*T_sym])
title('Transmitted Signal', 'Interpreter', 'latex','FontSize',20)
xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
ylabel('{\boldmath $s(t)$}', 'Interpreter', 'latex','FontSize',16)


% DETERMINE AVERAGE ENERGY PER BIT IN SIGNAL
E = 0;
for k = 1:1:length(s) % calcuate total energy in the signal using Riemann sum
    E = E + (s(k))^2*(t_step);
end
Eb = E/N; % obtain energy per bit by dividing the total energy by the number of bits


% SIMULATION OF THE PSD OF THE MODULATED SIGNAL
S = t_step .* fftshift( fft(s) ); % calculate DFT of s, center it around 0, scale it appropriately
P_s = (abs(S)).^2 / (N/16 * T_sym); % compute approximate PSD of s
f_step = 1 / (length(S) * t_step); % frequency step for approximating CTFT from DFT
f = -f_step*(length(S)/2) : f_step : f_step*(length(S)/2-1); % calculate frequencies corresponding to the spectrum, S

P_s_dB = 10.*log10(P_s); % PSD in decibel units

figure(4)
subplot(2,1,1) % plot PSD in absolute units
plot(f, P_s,'linewidth',2)
xlim([-1.5.*max(f_c) 1.5.*max(f_c)])
title('{\boldmath $P_{s}(f)$}', 'Interpreter', 'latex','FontSize',24)
xlabel('{\boldmath $f$} (W/Hz)', 'Interpreter', 'latex','FontSize',18)
ylabel('{\boldmath $P_{s}(f)$}', 'Interpreter', 'latex','FontSize',18)
hold off
subplot(2,1,2) % plot PSD in decibel units
plot(f, P_s_dB,'linewidth',2)
xlim([-1.5.*max(f_c) 1.5.*max(f_c)])
title('{\boldmath $P_{s}(f)$} (dB)', 'Interpreter', 'latex','FontSize',24)
xlabel('{\boldmath $f$} (dBW/Hz)', 'Interpreter', 'latex','FontSize',18)
ylabel('{\boldmath $P_{s}(f)$} (dB)', 'Interpreter', 'latex','FontSize',18)
hold off


% ADD NOISE TO TRANSMITTED SIGNAL TO REPRESENT RECEIVED SIGNAL
r = s + sqrt(noiseVar)*randn(1, length(s)); % add noise to the signal


% PLOT RECEIVED SIGNAL (BEFORE FILTERING)
figure(3)
subplot(2,1,2)
plot(t_s, r,'linewidth',2)
grid on
xlim([-5*T_sym (N/16+5)*T_sym])
title('Received Signal (Transmitted Signal Plus Noise)', 'Interpreter', 'latex','FontSize',20)
xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
ylabel('{\boldmath $r(t)$}', 'Interpreter', 'latex','FontSize',16)


% I-Q DOWNCONVERT EACH CARRIER
% each row of the following arrays represents a different carrier
d_I = zeros(4, length(r));
d_Q = zeros(4, length(r));

for carrier = [1 2 3 4]
    d_I(carrier,:) = 2.*r.*cos(2.*pi.*f_c(carrier)*t_s);
    d_Q(carrier,:) = -2.*r.*sin(2.*pi.*f_c(carrier)*t_s);
end


% PLOT DOWNCOVERTED BASEBAND OUTPUTS
figure(5)
for carrier = [1 2 3 4]
    subplot(2,2,carrier) % subplot for each carrier
    plot(t_s, d_I(carrier,:),'linewidth',2) % I-channel
    hold on
    plot(t_s, d_Q(carrier,:),'linewidth',2) % Q-channel
    grid on
    xlim([-5*T_sym (N/16+5)*T_sym])
    title(strcat("Baseband Output (Before Filtering) - Carrier ", num2str(carrier)), 'Interpreter', 'latex','FontSize',20)
    xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
    ystring = strcat("{\boldmath $d_{", num2str(carrier), "}(t)$}");
    ylabel(ystring, 'Interpreter', 'latex','FontSize',16)
    legend('I-channel', 'Q-channel')
end


% FILTER EACH DOWNCONVERTED CARRIER
h_padded = [h zeros(1, length(r)-length(h))]; % pad impulse response h with zeros for convolution
t_filtered = -2*k_T*T_sym : t_step : 2*((N/16-1)*T_sym+k_T*T_sym - t_step); % corresponding time samples for filtered signal

% each row of the following arrays represents a different carrier
y_I = zeros(4, 2*length(h_padded)-1);
y_Q = zeros(4, 2*length(h_padded)-1);

for carrier = [1 2 3 4] % filter each downconverted carrier
    y_I(carrier,:) = (1./sampPerSym).*conv(d_I(carrier,:), h_padded);
    y_Q(carrier,:) = (1./sampPerSym).*conv(d_Q(carrier,:), h_padded);
end


% DETERMINE RECEIVED SYMBOL VALUES
% each row of the following arrays represents a different carrier
carrier_i_rx_raw = zeros(4, N/16);
carrier_q_rx_raw = zeros(4, N/16);

% capture values in y_I and y_Q that occurred at the symbol times
for carrier = [1 2 3 4]
    for count = 1:1:N/16
        carrier_i_rx_raw(carrier, count) = y_I(carrier, abs(t_filtered - (count-1)*T_sym) < 0.1*t_step );
        carrier_q_rx_raw(carrier, count) = y_Q(carrier, abs(t_filtered - (count-1)*T_sym) < 0.1*t_step );
    end
end


% PLOT RECEIVED CONSTELLATIONS
figure(6)
for carrier = [1 2 3 4]
    subplot(2,2,carrier)
    plot(carrier_i_rx_raw(carrier,:), carrier_q_rx_raw(carrier,:),'xr','linewidth',2)
    title(strcat("RX Constellation - Carrier ", num2str(carrier)), 'Interpreter', 'latex','FontSize',20)
    xlabel('{\boldmath $I$}', 'Interpreter', 'latex','FontSize',16)
    ylabel('{\boldmath $Q$}', 'Interpreter', 'latex','FontSize',16)
end


% DECISION RULE TO BEST GUESS TRANSMITTED SYMBOLS
roundTargets = [(-A) (-A/3) (A/3) (A)]; % rounding targets

% round all elements of the below arrays to the above targets
carrier_i_rx = interp1(roundTargets, roundTargets, carrier_i_rx_raw,'nearest','extrap');
carrier_q_rx = interp1(roundTargets, roundTargets, carrier_q_rx_raw,'nearest','extrap');


% CONVERT RECEIVED SYMBOL VALUES TO BITS
% each row of the following arrays represents a different carrier
bits_carrier_rx = zeros(4, N/4);
carrier_bitstring_rx = strings(4, N/16);

for carrier = [1 2 3 4]
    for counter = 1:1:N/16
        % Decode constellation
        if ( (carrier_i_rx(carrier, counter) == -A) && (carrier_q_rx(carrier, counter) == -A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 0, 0, 0];
            carrier_bitstring_rx(carrier, counter) = "0000";
        elseif ( (carrier_i_rx(carrier, counter) == -A) && (carrier_q_rx(carrier, counter) == -A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 0, 0, 1];
            carrier_bitstring_rx(carrier, counter) = "0001";
        elseif ( (carrier_i_rx(carrier, counter) == -A) && (carrier_q_rx(carrier, counter) == A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 0, 1, 0];
            carrier_bitstring_rx(carrier, counter) = "0010";
        elseif ( (carrier_i_rx(carrier, counter) == -A) && (carrier_q_rx(carrier, counter) == A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 0, 1, 1];
            carrier_bitstring_rx(carrier, counter) = "0011";
        elseif ( (carrier_i_rx(carrier, counter) == -A/3) && (carrier_q_rx(carrier, counter) == -A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 1, 0, 0];
            carrier_bitstring_rx(carrier, counter) = "0100";
        elseif ( (carrier_i_rx(carrier, counter) == -A/3) && (carrier_q_rx(carrier, counter) == -A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 1, 0, 1];
            carrier_bitstring_rx(carrier, counter) = "0101";
        elseif ( (carrier_i_rx(carrier, counter) == -A/3) && (carrier_q_rx(carrier, counter) == A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 1, 1, 0];
            carrier_bitstring_rx(carrier, counter) = "0110";
        elseif ( (carrier_i_rx(carrier, counter) == -A/3) && (carrier_q_rx(carrier, counter) == A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [0, 1, 1, 1];
            carrier_bitstring_rx(carrier, counter) = "0111";
        elseif ( (carrier_i_rx(carrier, counter) == A) && (carrier_q_rx(carrier, counter) == -A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 0, 0, 0];
            carrier_bitstring_rx(carrier, counter) = "1000";
        elseif ( (carrier_i_rx(carrier, counter) == A) && (carrier_q_rx(carrier, counter) == -A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 0, 0, 1];
            carrier_bitstring_rx(carrier, counter) = "1001";
        elseif ( (carrier_i_rx(carrier, counter) == A) && (carrier_q_rx(carrier, counter) == A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 0, 1, 0];
            carrier_bitstring_rx(carrier, counter) = "1010";
        elseif ( (carrier_i_rx(carrier, counter) == A) && (carrier_q_rx(carrier, counter) == A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 0, 1, 1];
            carrier_bitstring_rx(carrier, counter) = "1011";
        elseif ( (carrier_i_rx(carrier, counter) == A/3) && (carrier_q_rx(carrier, counter) == -A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 1, 0, 0];
            carrier_bitstring_rx(carrier, counter) = "1100";
        elseif ( (carrier_i_rx(carrier, counter) == A/3) && (carrier_q_rx(carrier, counter) == -A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 1, 0, 1];
            carrier_bitstring_rx(carrier, counter) = "1101";
        elseif ( (carrier_i_rx(carrier, counter) == A/3) && (carrier_q_rx(carrier, counter) == A) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 1, 1, 0];
            carrier_bitstring_rx(carrier, counter) = "1110";
        elseif ( (carrier_i_rx(carrier, counter) == A/3) && (carrier_q_rx(carrier, counter) == A/3) )
            bits_carrier_rx(carrier, (counter*4-3):(counter*4)) = [1, 1, 1, 1];
            carrier_bitstring_rx(carrier, counter) = "1111";
        end
    end
end


% PLOT FILTERED BASEBAND OUTPUTS AND BITS
figure(7)
for carrier = [1 2 3 4]
    subplot(2,2,carrier) % suplot for each carrier
    plot(t_filtered, y_I(carrier,:),'linewidth',2) % I-channel
    hold on
    plot(t_filtered, y_Q(carrier,:),'linewidth',2) % Q-channel
    grid on
    xlim([-5*T_sym (N/16+5)*T_sym])
    ylim([-2*A 2*A])
    title(strcat("Filtered Baseband Output - Carrier ", num2str(carrier)), 'Interpreter', 'latex','FontSize',20)
    xlabel('{\boldmath $t$} (s)', 'Interpreter', 'latex','FontSize',16)
    ystring = strcat("{\boldmath $y_{", num2str(carrier), "}(t)$}");
    ylabel(ystring, 'Interpreter', 'latex','FontSize',16)
    legend('I-channel', 'Q-channel')
    for cnt = 1:1:N/16
        if carrier_bitstring_rx(carrier, cnt) == carrier_bitstring_tx(carrier,cnt)
            text((cnt-1)*T_sym, 1.5*A, carrier_bitstring_rx(carrier,cnt), 'HorizontalAlignment', 'center')
        else
            text((cnt-1)*T_sym, 1.5*A, carrier_bitstring_rx(carrier,cnt), 'HorizontalAlignment', 'center', 'Color', 'r')
        end
    end
end


% COMBINE RECEIVED BITS FROM EACH CARRIER TO SINGLE BITSTREAM
bits_rx = zeros(1, N);

for counter = 1:1:N/16
    bits_rx(counter*16-15:counter*16-12) = bits_carrier_rx(1, counter*4-3:counter*4);
    bits_rx(counter*16-11:counter*16-8) = bits_carrier_rx(2, counter*4-3:counter*4);
    bits_rx(counter*16-7:counter*16-4) = bits_carrier_rx(3, counter*4-3:counter*4);
    bits_rx(counter*16-3:counter*16) = bits_carrier_rx(4, counter*4-3:counter*4);
end


%CALCULATE AND DISPLAY BIT ERROR RATE
BER = mean(bits_tx ~= bits_rx); % calculate BER by comparing input and output bit vectors
disp(strcat("BER = ", num2str(BER))) % display BER


%CALCULATE AND DISPLAY Eb/N0
N0 = 2*noiseVar/sampPerSecond;
disp(strcat("Eb/N0 = ", num2str(Eb/N0)))
