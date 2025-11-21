%%%%%%%%%%%%%%%%%%%%%%%%% Define Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
t = linspace(-0.001,2.001,1000); % Time Vector
fc = 10000;                      % Carrier Frequency
Ac = 2;                          % Amplitude
Ka = 0.5;                        % Modulation Index


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Signals %%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the first message message m1(t)
m1_t = zeros(1,length(t));
for i = 1:length(t)
    if t(i) >= 0 && t(i) < 0.5
        m1_t(i) = -2*t(i); 
    elseif t(i) >= 0.5 && t(i) < 1.5
        m1_t(i) = (-2*t(i)+2); 
    elseif t(i) >= 1.5 && t(i) < 2
        m1_t(i) = (-2*t(i) + 4);
    end
end

% Generate the second message message m2(t)
m2_t = zeros(1,length(t));
for i = 1:length(t)
    if t(i) == 0
        m2_t(i) = 1;
    elseif t(i) >= 0 && t(i) < 0.5
        m2_t(i) = 1;
    elseif t(i) >= 0.5 && t(i) < 1
        m2_t(i) = 0.5;
    elseif t(i) >= 1 && t(i) <= 1.5
        m2_t(i) = -0.5;
    elseif t(i) >= 1.5 && t(i) <= 2
        m2_t(i) = -1;
    elseif t(i) > 2
        m2_t(i) = 0;
    end
end

% Generate the carrier signal
carrier = Ac * cos(2*pi*fc*t);


%%%%%%%%%%%%%%%%%%%%%%%%%% Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modulated signals for m1_t and m2_t
s_t1 = (1 + Ka * m1_t) .* carrier;
s_t2 = (1 + Ka * m2_t) .* carrier;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

% Plot the first message and its modulation
subplot(2,2,1);
plot(t,m1_t, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)'); 
ylabel('m1(t)');
title('First Message');
axis([0 2 -1 1]); 
yticks([-1 -0.5 0 0.5 1]);
xticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
grid on;

subplot(2,2,2);
plot(t, s_t1, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
xlim([0,1.99]);
ylabel('s(t)');
title('Modulated Signal for m_1(t)');
grid on;

% Plot the second message and its modulation
subplot(2,2,3);
plot(t,m2_t, 'r', 'LineWidth', 1.5);
xlabel('Time (ms)'); 
ylabel('m2(t)');
title('Second Message');
axis([0 2 -1 1]); 
xlim([-0.01, 2.01])
yticks([-1 -0.5 0 0.5 1]);
grid on;

subplot(2,2,4);
plot(t, s_t2, 'r', 'LineWidth', 1.5);
xlim([0,1.99]);
xlabel('Time (ms)');
ylabel('s(t)');
title('Modulated Signal for m_2(t)');
grid on;
