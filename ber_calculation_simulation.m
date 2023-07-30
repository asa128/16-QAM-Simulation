close all % close all figures
clear % clear workspace

Eb_N0_dB_min = -40; % minimum dB value for Eb/N0
Eb_N0_dB_max = 15; % maximum dB value for Eb/N0
numPoints = 100; % number of dB values to run calculation at

numTrials = 10; % number of trials for which to average each Eb/N0 case over
N = 20*16; % number of bits per trial (must be a multiple of 16)

A = 0.67; % peak amplitude of I or Q channel
sampPerBit = 64; % number of samples per bit
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
sampPerSecond = D * sampPerBit; % sample rate
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

Eb_N0_array_dB = linspace(Eb_N0_dB_min, Eb_N0_dB_max, numPoints); % Define Eb/N0 array in decibels

Eb_N0_array = 10.^(Eb_N0_array_dB./10); % Convert Eb/N0 array to abosulte units

BER_array = zeros(1, length(Eb_N0_array)); % create vector for bit error rates, same length as Eb/N0 vector

% GENERATE PULSE AND RANDOM BITS
[h, t] = rrcpulse(k_T, T_sym, sampPerBit, r); % generate RRC pulse

% LOOP OVER POSSIBLE Eb/N0 VALUES
for value = 1:1:length(Eb_N0_array)
    
    Eb_N0 = Eb_N0_array(value); % obtain the Eb/N0 value for the current loop index
    
    BER_trials = zeros(1,numTrials); % vector of length trials to store the BER of each trial
    
    % LOOP OVER NUMBER OF TRIALS FOR WHICH TO AVERAGE EACH Eb/N0 CASE OVER
    for num = 1:1:numTrials
        bits_tx = randbin(N, 0, 1); % generate vector a of N random bits

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
                for i = 1:1:(2*k_T*sampPerBit)
                    carrier_i_tx_bb(carrier, i+k*sampPerBit) = carrier_i_tx_bb(carrier, i+k*sampPerBit) + carrier_i_tx(carrier, k+1)*h(i);
                    carrier_q_tx_bb(carrier, i+k*sampPerBit) = carrier_q_tx_bb(carrier, i+k*sampPerBit) + carrier_q_tx(carrier, k+1)*h(i);
                end
            end
        end
        
        
        % GENERATE MODULATED SIGNALS
        s_mod = zeros(4, length(t_s)); % each row represents a different carrier
        
        for carrier = [1 2 3 4] % modulate each carrier separately
            s_mod(carrier,:) = carrier_i_tx_bb(carrier,:).*cos(2.*pi.*f_c(carrier).*t_s) - carrier_q_tx_bb(carrier,:).*sin(2.*pi.*f_c(carrier).*t_s);
        end
        
        s = sum(s_mod, 1); % sum the modulated carriers to form the signal
        
        
        % DETERMINE AVERAGE ENERGY PER BIT IN SIGNAL
        E = 0;
        for k = 1:1:length(s) % calcuate total energy in the signal using Riemann sum
            E = E + (s(k))^2*(t_step);
        end
        Eb = E/N; % obtain energy per bit by dividing the total energy by the number of bits
        
        % ADD NOISE TO TRANSMITTED SIGNAL TO REPRESENT RECEIVED SIGNAL
        N0 = 1/(Eb_N0)*Eb;
        r = s + sqrt(N0*sampPerSecond/2)*randn(1, length(s)); % add noise to the signal
        
        
        % I-Q DOWNCONVERT EACH CARRIER
        % each row of the following arrays represents a different carrier
        d_I = zeros(4, length(r));
        d_Q = zeros(4, length(r));
        
        for carrier = [1 2 3 4]
            d_I(carrier,:) = 2.*r.*cos(2.*pi.*f_c(carrier)*t_s);
            d_Q(carrier,:) = -2.*r.*sin(2.*pi.*f_c(carrier)*t_s);
        end
        
        
        % FILTER EACH DOWNCONVERTED CARRIER
        h_padded = [h zeros(1, length(r)-length(h))]; % pad impulse response h with zeros for convolution
        t_filtered = -2*k_T*T_sym : t_step : 2*((N/16-1)*T_sym+k_T*T_sym - t_step); % corresponding time samples for filtered signal
        
        % each row of the following arrays represents a different carrier
        y_I = zeros(4, 2*length(h_padded)-1);
        y_Q = zeros(4, 2*length(h_padded)-1);
        
        for carrier = [1 2 3 4] % filter each downconverted carrier
            y_I(carrier,:) = (1./sampPerBit).*conv(d_I(carrier,:), h_padded);
            y_Q(carrier,:) = (1./sampPerBit).*conv(d_Q(carrier,:), h_padded);
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
        

        % COMBINE RECEIVED BITS FROM EACH CARRIER TO SINGLE BITSTREAM
        bits_rx = zeros(1, N);
        
        for counter = 1:1:N/16
            bits_rx(counter*16-15:counter*16-12) = bits_carrier_rx(1, counter*4-3:counter*4);
            bits_rx(counter*16-11:counter*16-8) = bits_carrier_rx(2, counter*4-3:counter*4);
            bits_rx(counter*16-7:counter*16-4) = bits_carrier_rx(3, counter*4-3:counter*4);
            bits_rx(counter*16-3:counter*16) = bits_carrier_rx(4, counter*4-3:counter*4);
        end
        
        
        BER_trials(num) = mean(bits_tx ~= bits_rx); % calculate BER by comparing input and output bit vectors
        
    end

    %CALCULATE AVERAGE BER
    BER_array(value) = mean(BER_trials); % calculate average bit error rate over the trials for the current Eb/N0 value

    disp(strcat("Progress = ", num2str(100.0*value/length(Eb_N0_array)), "%")) % display progress in loop

end

figure(1)

% PLOT WITH A LINEAR SCALE FOR BER (P_e) AXIS
subplot(1,2,1)
%P_e = (3/8)*erfc(sqrt(0.4*Eb_N0_array))+(1/4)*erfc(3*sqrt(0.4*Eb_N0_array))-(1/8)*erfc(5*sqrt(0.4*Eb_N0_array));
P_e = (3/4)*qfunc(sqrt(0.8*Eb_N0_array))+(1/2)*qfunc(3*sqrt(0.8*Eb_N0_array))-(1/4)*qfunc(5*sqrt(0.8*Eb_N0_array));
plot(Eb_N0_array_dB, P_e, 'linewidth', 2)
hold on
plot(Eb_N0_array_dB, BER_array,'xr','linewidth',2)
title('Expected vs. Actual BER', 'Interpreter', 'latex','FontSize',20)
xlabel('{\boldmath $\frac{E_{b}}{N_{0}}$} (dB)', 'Interpreter', 'latex','FontSize',16)
xlim([Eb_N0_dB_min Eb_N0_dB_max])
ylabel('{\boldmath $BER$}', 'Interpreter', 'latex','FontSize',16)
legend('Expected BER', 'Actual BER')


% PLOT WITH A LOG SCALE FOR BER (P_e) AXIS
subplot(1,2,2)
semilogy(Eb_N0_array_dB, P_e, 'linewidth', 2)
hold on
semilogy(Eb_N0_array_dB, BER_array,'xr','linewidth',2)
title('Expected vs. Actual BER', 'Interpreter', 'latex','FontSize',20)
xlabel('{\boldmath $\frac{E_{b}}{N_{0}}$} (dB)', 'Interpreter', 'latex','FontSize',16)
xlim([Eb_N0_dB_min Eb_N0_dB_max])
ylabel('{\boldmath $BER$}', 'Interpreter', 'latex','FontSize',16)
legend('Expected BER', 'Actual BER')
