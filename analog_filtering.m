%% Audio import

clc;
[a, Fs] = audioread("Tone.wav");

% Waveform plot

d = length(a)/Fs;
t = linspace(0, d, length(a));

figure();
plot(t, a);
title("Audio Waveform")
xlabel("time [s]");
ylabel("Amplitude");
grid on, grid minor;

%% Frecuency Spectrum 

A = fftshift(fft(a));

f = linspace(-Fs/2, Fs/2, length(A));

mag_A = abs(A);

mag_A_N = mag_A/max(mag_A);

figure();
plot(f, mag_A_N);
title("Audio Frecuency Spectrum");
xlabel("Frecuency [Hz]");
ylabel("Amplitude");
grid on, grid minor;

ax = gca;
ax.XAxis.Exponent = 3;

%% Filter Design

Wp = [400, 480].*2*pi;
Ws = [340, 540].*2*pi;
Rp = 0.5;
Rs = 40;

%% Butterworth
[n, Wn] = buttord(Wp, Ws, Rp, Rs, 's');
[num_b, den_b] = butter(n, Wn, 'bandpass', 's');

fn = Wn./(2*pi);
fprintf("Butterworth design: \n " ...
        + "-Order: %d\n "...
        + "-Cutoff Freq: %.2f, %.2f [Hz]\n\n", n, fn(1), fn(2));
    
%% Chebyshev Type I

[n, Wn] = cheb1ord(Wp, Ws, Rp, Rs, 's');
[num_c1, den_c1] = cheby1(n, Rp, Wn, 'bandpass', 's');

fn = Wn./(2*pi);
fprintf("Chebyshev Type I: \n " ...
        + "-Order: %d\n "...
        + "-Cutoff Freq: %.2f, %.2f [Hz]\n\n", n, fn(1), fn(2));

%% Chebyshev Type II

[n, Wn] = cheb2ord(Wp, Ws, Rp, Rs, 's');
[num_c2, den_c2] = cheby2(n, Rs, Wn, 'bandpass', 's');

fn = Wn./(2*pi);
fprintf("Chebyshev Type II: \n " ...
        + "-Order: %d\n "...
        + "-Cutoff Freq: %.2f, %.2f [Hz]\n\n", n, fn(1), fn(2));
    
%% Filter frequency spectrum     

[h_b,  wout] = freqs(num_b,  den_b,  f*2*pi);
[h_c1, wout] = freqs(num_c1, den_c1, f*2*pi);
[h_c2, wout] = freqs(num_c2, den_c2, f*2*pi);

fp = Wp./(2*pi);
fs = Ws./(2*pi);

% dB plot
figure();
subplot(2, 1, 1);
plot(f, mag2db(abs(h_b)), 'LineWidth', 1.2);
hold on;
plot(f, mag2db(abs(h_c1)), 'LineWidth', 1.2);
plot(f, mag2db(abs(h_c2)), 'LineWidth', 1.2);

xline(fp(1), '--k', "Fp");
xline(fp(2), '--k', "Fp", 'LabelHorizontalAlignment', 'left');

xline(fs(1), '-.k', "Fs");
xline(fs(2), '-.k', "Fs", 'LabelHorizontalAlignment', 'left');

yline(-Rp, '--k', "Rp", 'LabelHorizontalAlignment', 'left');
yline(-Rs, '--k', "Rs", 'LabelHorizontalAlignment', 'left');

title("Filter Frequency spectrum")
xlabel("Frequency [Hz]")
ylabel("Amplitude [dB]")
grid on, grid minor;

legend("Butterworth", "Chebyshev I", "Chebyshev II")

ax = gca;
ax.XAxis.Exponent = 3;

xlim([0, 1e3]);
ylim([-50, 10]); 

% Amplitude Plot
subplot(2, 1, 2);
plot(f, (abs(h_b)), 'LineWidth', 1.2);
hold on;
plot(f, (abs(h_c1)), 'LineWidth', 1.2);
plot(f, (abs(h_c2)), 'LineWidth', 1.2);

xline(fp(1), '--k', "Fp");
xline(fp(2), '--k', "Fp", 'LabelHorizontalAlignment', 'left');

xline(fs(1), '-.k', "Fs");
xline(fs(2), '-.k', "Fs", 'LabelHorizontalAlignment', 'left');

title("Filter Frequency spectrum")
xlabel("Frequency [Hz]")
ylabel("Amplitude ")
grid on, grid minor;

legend("Butterworth", "Chebyshev I", "Chebyshev II")

ax = gca;
ax.XAxis.Exponent = 3;

xlim([0, 1e3]);

%% Comparison 

figure();
plot(f, mag_A_N, 'LineWidth', 1.8);
hold on;
plot(f, (abs(h_b)), 'LineWidth', 1.1);
plot(f, (abs(h_c1)), 'LineWidth', 1.1);
plot(f, (abs(h_c2)), 'LineWidth', 1.1);

xlim([0, 2e3]);
legend("Audio", "Butterworth", "Chebyshev I", "Chebyshev II");
ax = gca;
ax.XAxis.Exponent = 3;
title("Frequency spectrum comparison");
xlabel("Frequency [Hz]")
ylabel("Amplitude")
grid on, grid minor;

%% Select Filter 
% Filter Chebyshev Type I selected

num = num_c1;
den = den_c1;

%% Filtering

lpf_sys = tf(num, den);
a_f = lsim(lpf_sys, a, t);

%% Filtered audio spectrum

A_f = fftshift(fft(a_f));

A_f_N = abs(A_f)./max(abs(A_f));

figure();
plot(f, A_f_N);
title("Filtered audio spectrum")
xlabel("Frequency [Hz]");
ylabel("Amplitude");
grid on, grid minor;
ax = gca;
ax.XAxis.Exponent = 3;

%% Audio playing

sound(a, Fs);
pause(d+0.5);

sound(a_f, Fs);
pause(d+0.5);

