%%%%%%%%%%%%%%%%%%%%%%%%% Define Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
Am = 1;                                          % Message Amplitude
B = 1000;                                        % Message Bandwidth
Ac = 2;                                          % Carrier Amplitude
fc = 5000;                                       % Carrier Frequency
filterFrequency = 2*fc;                          % Filter Frequency
fs = (B+fc)*10;                                  % Sampling Rate
t = -0.005:1/fs:0.005;                           % Time vector
f = (-fs/2):(fs/length(t)):(fs/2-fs/length(t));  % Frequency vector
mfun = sinc(B*t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Signals %%%%%%%%%%%%%%%%%%%%%%%%%

% Generate message signal m3(t)
message = Am * mfun;
shiftedMessage = imag(hilbert(message));

% Generate the carrier signal
carrier = Ac * cos(2*pi*fc*t);
shiftedCarrier = imag(hilbert(carrier));


%%%%%%%%%%%%%%%%%%%%%%%%%% Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modulate the message signal m3(t) onto the carrier
[USB, LSB] = ssb_modulation(message, shiftedMessage, carrier, shiftedCarrier);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% spectrum of both signals
USB_spectrum = fft(USB);
LSB_spectrum = fft(LSB);

% Compute the two-sided spectrum
USB_spectrum_real = fftshift(real(USB_spectrum));
USB_spectrum_imag = fftshift(imag(USB_spectrum));

LSB_spectrum_real = fftshift(real(LSB_spectrum));
LSB_spectrum_imag = fftshift(imag(LSB_spectrum));


%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[demodUSB, demodLSB] = ssb_demodulation(USB, LSB, carrier, Ac);

outputUSB = lowPassFilter(demodUSB, t, filterFrequency);
outputLSB = lowPassFilter(demodLSB, t, filterFrequency);


%%%%%%%%%%%%%%%%%%%%%%%% Demodulation (f1 & f2) %%%%%%%%%%%%%%%%%%%%%%%%

c1 = Ac * cos(2*pi*(fc+0.1*B)*t);
c2 = Ac * cos(2*pi*(fc-0.2*B)*t);

[demodUSB_c1, demodLSB_c1] = ssb_demodulation(USB, LSB, c1, Ac);
[demodUSB_c2, demodLSB_c2] = ssb_demodulation(USB, LSB, c2, Ac);

outputUSB_c1 = lowPassFilter(demodUSB_c1, t, filterFrequency);
outputUSB_c2 = lowPassFilter(demodUSB_c2, t, filterFrequency);

outputLSB_c1 = lowPassFilter(demodLSB_c1, t, filterFrequency);
outputLSB_c2 = lowPassFilter(demodLSB_c2, t, filterFrequency);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotVectors({message, shiftedMessage}, ...
            t, 1, ...
            {'Time (s)'}, ...
            {'m(t)'}, ...
            {'Message Signal', 'Message Shifted Signal'}, ...
            {}, ...
            {});

plotVectors({carrier, shiftedCarrier}, ...
            t, 1, ...
            {'Time (s)'}, ...
            {'c(t)'}, ...
            {'Carrier Signal', 'Carrier Shifted Signal'}, ...
            {[-0.001, 0.001]}, ...
            {[-2.5, 2.5]});

plotVectors({USB, LSB}, ...
            t, 1, ...
            {'Time (s)'}, ...
            {'USB(t) & LSB(t)'}, ...
            {'Modulated Signals', 'LSP Signal'}, ...
            {}, ...
            {});

plotVectors({USB, LSB}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'USB(t)', 'LSB(t)'}, ...
            {'USP Signal', 'LSP Signal'}, ...
            {}, ...
            {});

plotVectors({USB_spectrum_real, USB_spectrum_imag}, ...
            f, 2, ...
            {'Freq (Hz)', 'Freq (Hz)'}, ...
            {'Real', 'Imaginary'}, ...
            {'Real part of the USB Spectrum', 'Imaginary part of the USB Spectrum'}, ...
            {}, ...
            {});

plotVectors({LSB_spectrum_real, LSB_spectrum_imag}, ...
            f, 2, ...
            {'Freq (Hz)', 'Freq (Hz)'}, ...
            {'Real', 'Imaginary'}, ...
            {'Real part of the LSB Spectrum', 'Imaginary part of the LSB Spectrum'}, ...
            {}, ...
            {});

figUSB = plotVectors({demodUSB, outputUSB}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From USB', 'Filtered Message From USB'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodUSB, t, figUSB));

figLSB = plotVectors({demodLSB, outputLSB}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From LSB', 'Filtered Message From LSB'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodLSB, t, figLSB));

figUSB_c1 = plotVectors({demodUSB_c1, outputUSB_c1}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From USB (c1)', 'Filtered Message From USB (c1)'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodUSB_c1, t, figUSB_c1));

figUSB_c2 = plotVectors({demodUSB_c2, outputUSB_c2}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From USB (c2)', 'Filtered Message From USB (c2)'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodUSB_c2, t, figUSB_c2));

figUSB_c1 = plotVectors({demodUSB_c1, outputLSB_c1}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From LSB (c1)', 'Filtered Message From LSB (c1)'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodUSB_c1, t, figUSB_c1));

figUSB_c2 = plotVectors({demodUSB_c2, outputLSB_c2}, ...
            t, 2, ...
            {'Time (s)', 'Time (s)'}, ...
            {'m(t)', 'm(t)'}, ...
            {'Non Filtered Message From LSB (c2)', 'Filtered Message From LSB (c2)'}, ...
            {}, ...
            {[], [-0.8, 1.1]});
createSlider(@(src, ~) updatePlotYData(src, demodUSB_c2, t, figUSB_c2));


%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [USB, LSB] = ssb_modulation(message, shiftedMessage, carrier, shiftedCarrier)
    USB = (message .* carrier) - (shiftedMessage .* shiftedCarrier);
    LSB = (message .* carrier) + (shiftedMessage .* shiftedCarrier);
end


function [demodUSB, demodLSB] = ssb_demodulation(USB, LSB, carrier, Ac)
    demodUSB = USB .* carrier .* (2/Ac^2);
    demodLSB = LSB .* carrier .* (2/Ac^2);
end


function output = lowPassFilter(originalMessage, t, Wc)
    % Define parameters
    s = tf('s');
    Rf = 33000;
    Ri = 33000;

    % Calculate Laplace transform of the transfer function
    H1 = -(Rf/Ri) * Wc / (s + Wc);
    H2 = -(Rf/Ri) * Wc^2 / (s^2 + s*Wc*sqrt(2) + Wc^2);

    % Simulate the system response
    output = lsim(H1*H2, originalMessage, t);
end


function updatePlotYData(slider, originalMessage, t, fig)
    new_Wc = slider.Value;
    output = lowPassFilter(originalMessage, t, new_Wc);
    updatePlot(fig, output);
    disp(new_Wc);
end


function updatePlot(fig, output)
    axes_handle = findobj(fig, 'Type', 'Axes');
    lines = findobj(axes_handle, 'Type', 'Line');
    line_handle = lines(1);
    set(line_handle, 'YData', output);
    drawnow;
end


function fig = plotVectors(vectors, t, mode, xlabels, ylabels, titles, xlimits, ylimits)
    
    validateMode(mode);

    if mode == 1
        fig = plotAllVectors(vectors, t, xlabels, ylabels, titles, xlimits, ylimits);
    elseif mode == 2
        fig = plotVectorsInSubplots(vectors, t, xlabels, ylabels, titles, xlimits, ylimits);
    end
end


function validateMode(mode)
    if ~(mode == 1 || mode == 2)
        error('Invalid mode. Mode must be 1 or 2.');
    end
end


function fig = plotAllVectors(vectors, t, xlabels, ylabels, titles, xlimits, ylimits)
    fig = figure;
    hold on;
    colors = {'r', 'b', 'g', 'y'};
    for i = 1:length(vectors)
        plot(t, vectors{i}, colors{i}, 'LineWidth', 1.5);
    end
    hold off;
    xlabel(xlabels{1});
    ylabel(ylabels{1});
    title(titles{1});
    legend(titles);
    grid on;
    
    % Set x-axis and y-axis limits for each vector
    if nargin >= 6 && ~isempty(xlimits)
        xlim(xlimits{1});
    end
    if nargin >= 7 && ~isempty(ylimits)
        ylim(ylimits{1});
    end
end


function fig = plotVectorsInSubplots(vectors, t, xlabels, ylabels, titles, xlimits, ylimits)
    numVectors = length(vectors);
    numRows = ceil(sqrt(numVectors));
    numCols = ceil(numVectors / numRows);
    colors = {'r', 'b', 'g', 'y'};
    fig = figure;
    for i = 1:numVectors
        subplot(numRows, numCols, i);
        plot(t, vectors{i}, colors{i}, 'LineWidth', 1.5);
        xlabel(xlabels{i});
        ylabel(ylabels{i});
        title(titles{i});
        grid on;
        
        % Set x-axis and y-axis limits for each vector
        if nargin >= 6 && ~isempty(xlimits) && ~isempty(xlimits{i})
            xlim(xlimits{i});
        end
        if nargin >= 7 && ~isempty(ylimits) && ~isempty(ylimits{i})
            ylim(ylimits{i});
        end
    end
end


function slider = createSlider(callback)
    fc = evalin('base', 'fc');
    slider = uicontrol('Style', 'slider', ...
        'Min', 0, 'Max', 50*fc, 'Value', 2*fc, ...
        'Units', 'normalized', 'Position', [0.05 0.01 0.9 0.05], ...
        'SliderStep', [1/1000, 1/1000], ...
        'BackgroundColor', [0.8 0.8 0.8], ...
        'ForegroundColor', [0.2 0.6 1], ...
        'Callback', callback);
end