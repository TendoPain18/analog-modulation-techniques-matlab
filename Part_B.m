%%%%%%%%%%%%%%%%%%%%%%%%% Define Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
t = linspace(-0.001,2.001,1000); % Time Vector
Ac = 5;                 % Carrier Amplitude
fc = 10e3;              % Carrier Frequency


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Signals %%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the first message message m1(t)
x_m1 = zeros(size(t));
for i = 1:length(t)
    if t(i) >= 0 && t(i) < 0.5
        x_m1(i) = -2*t(i); 
    elseif t(i) >= 0.5 && t(i) < 1.5
        x_m1(i) = (-2*t(i)+2); 
    elseif t(i) >= 1.5 && t(i) < 2
        x_m1(i) = (-2*t(i) + 4); 
    end
end

% Generate the second message message m2(t)
x_m2 = zeros(size(t));
for i = 1:length(t)
    if t(i) == 0
        x_m2(i) = 1;
    elseif t(i) >= 0 && t(i) < 0.5
        x_m2(i) = 1;
    elseif t(i) >= 0.5 && t(i) < 1
        x_m2(i) = 0.5;
    elseif t(i) >= 1 && t(i) <= 1.5
        x_m2(i) = -0.5;
    elseif t(i) >= 1.5 && t(i) <= 2
        x_m2(i) = -1;
    elseif t(i) > 2
        x_m2(i) = 0;
    end
end

% Generate the carrier signal
carrier = Ac * cos(2*pi*fc*t);


%%%%%%%%%%%%%%%%%%%%%%%%%% Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generation of Modulated Signals
Modulated_signal_1 = x_m1 .* carrier;
Modulated_signal_2 = x_m2 .* carrier;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Modulated Signals Multiplied By Carrier
demodulated_signal_1 = Modulated_signal_1 .* carrier * (1/Ac^2); 
demodulated_signal_2 = Modulated_signal_2 .* carrier * (1/Ac^2);

S1_filtered = lowPassFilter(demodulated_signal_1, t, 60);
S2_filtered = lowPassFilter(demodulated_signal_2, t, 60);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% Plot the first message
subplot(2, 1, 1);
plot(t, x_m1, 'b', 'LineWidth', 1.5);
xlabel('Time (msec)');
ylabel('m1(t)');
title('Message 1 Plotting');
grid on;

% Plot the second message
subplot(2, 1, 2);
plot(t, x_m2, 'r', 'LineWidth', 1.5);
xlabel('Time (msec)');
ylabel('m2(t)');
title('Message 2 Plotting');
grid on;

% Plot m1(t) Signal and its Modulated Signal 1
figure;
subplot(2, 1, 1);
plot(t, x_m1, 'b', 'LineWidth', 1.5);
title('m1(t) Signal');
xlabel('Time (msec)');
ylabel('m1(t)');
grid on;

subplot(2, 1, 2);
plot(t, Modulated_signal_1, 'r', 'LineWidth', 1.5);
title('Modulated Signal 1');
xlabel('Time (msec)');
ylabel('m1(t) modulated');
grid on;

% Plot m2(t) Signal and its Modulated Signal 2
figure;
subplot(2, 1, 1);
plot(t, x_m2, 'b', 'LineWidth', 1.5);
title('m2(t) Signal');
xlabel('Time (msec)');
ylabel('m2(t)');
grid on;

subplot(2, 1, 2);
plot(t, Modulated_signal_2, 'r', 'LineWidth', 1.5);
title('Modulated Signal 2');
xlabel('Time (msec)');
ylabel('m2(t) modulated');
grid on;

% Plot Modulated Signals Multiplied By The Carrier
figure;
subplot(2, 1, 1);
plot(t, demodulated_signal_1, 'b', 'LineWidth', 1.5);
title('Demodulated Signal 1');
xlabel('Time (msec)');
ylabel('m1(t) demodulated');
grid on;

subplot(2, 1, 2);
plot(t, demodulated_signal_2, 'r', 'LineWidth', 1.5);
title('Demodulated Signal 2');
xlabel('Time (msec)');
ylabel('m2(t) demodulated');
grid on;


% Plot Modulated Signals Multiplied By The Carrier
fig = figure;
subplot(2, 1, 1);
plot(t, S1_filtered, 'b', 'LineWidth', 1.5);
title('Demodulated Signal 1 Filtered');
xlabel('Time (msec)');
ylabel('m1(t) filtered');
grid on;

subplot(2, 1, 2);
plot(t, S2_filtered, 'r', 'LineWidth', 1.5);
title('Demodulated Signal 2 Filtered');
xlabel('Time (msec)');
ylabel('m2(t) filtered');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = lowPassFilter(originalMessage, t, Wc)
    % Define parameters
    s = tf('s');

    % Calculate Laplace transform of the transfer function
    H1 = Wc / (s + Wc);
    H2 = Wc^2 / (s^2 + s*Wc*sqrt(2) + Wc^2);

    % Simulate the system response
    output = lsim(H1*H2, originalMessage, t);
end