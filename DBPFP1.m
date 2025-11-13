% ========================================================================
% IMU Power Line Denoising + Bandpass Filter
% Metode: FFT Domain Power Line Removal + Butterworth Bandpass Filter
% ========================================================================

clear; clc; close all;

%% 1. KONFIGURASI PARAMETER
filename = 'S1P1_G.csv';   % Ganti dengan file CSV Anda

fs = 100;                  % Sampling frequency (Hz)
powerline_freq = 50;       % Power line frequency (Hz)
notch_bandwidth = 2;       % Bandwidth frekuensi yang dihapus (Hz)

% Parameter Bandpass Filter
bpf_low = 0.5;             % Cutoff bawah (Hz)
bpf_high = 20;             % Cutoff atas (Hz)
bpf_order = 4;             % Order filter

fprintf('=== FILTER AKTIF ===\n');
fprintf('1. FFT Denoising Power Line (%d Hz)\n', powerline_freq);
fprintf('2. Butterworth Bandpass Filter (%.1f–%.1f Hz)\n\n', bpf_low, bpf_high);

%% 2. BACA DATA CSV
fprintf('Membaca data dari %s...\n', filename);
try
    data = readtable(filename);
    fprintf('Data berhasil dibaca (%d sampel, %d kolom)\n', height(data), width(data));
catch
    error('Gagal membaca file CSV. Pastikan file ada dan formatnya benar.');
end

%% 3. DETEKSI KOLOM IMU
col_names = data.Properties.VariableNames;
accel_cols = find(contains(lower(col_names), {'acc', 'accel'}));
gyro_cols = find(contains(lower(col_names), {'gyro', 'gyr'}));
all_sensor_cols = [accel_cols, gyro_cols];

fprintf('\nKolom Accelerometer: '); disp(col_names(accel_cols));
fprintf('Kolom Gyroscope: '); disp(col_names(gyro_cols));

if isempty(all_sensor_cols)
    error('Tidak ada kolom sensor IMU yang terdeteksi.');
end

%% 4. DESAIN BANDPASS FILTER
nyquist_freq = fs/2;
fprintf('\n=== Desain Bandpass Filter ===\n');
[b_bpf, a_bpf] = butter(bpf_order, [bpf_low bpf_high]/nyquist_freq, 'bandpass');
[h_bpf, w_bpf] = freqz(b_bpf, a_bpf, 1024, fs);
fprintf('Bandpass Filter: %.2f–%.2f Hz (order %d)\n', bpf_low, bpf_high, bpf_order);

%% 5. PROSES FILTERING (FFT + BPF)
data_denoised = data;
data_filtered = data;
fprintf('\nMemulai proses filtering...\n');

for i = all_sensor_cols
    sig = table2array(data(:, i));
    
    % === (1) FFT Domain Power Line Denoising ===
    N = length(sig);
    X = fft(sig);
    freq = (0:N-1) * fs / N;
    mask = ones(size(X));

    % Hilangkan 50 Hz dan harmonik (100 Hz, 150 Hz)
    for k = 1:3
        target = k * powerline_freq;
        mask(abs(freq - target) < notch_bandwidth) = 0;
        mask(abs(freq - (fs - target)) < notch_bandwidth) = 0;
    end
    
    X_filtered = X .* mask;
    sig_denoised = real(ifft(X_filtered));

    % === (2) Bandpass Filter ===
    sig_bpf = filtfilt(b_bpf, a_bpf, sig_denoised);
    
    % Simpan hasilnya ke tabel
    data_denoised.(col_names{i}) = sig_denoised;
    data_filtered.(col_names{i}) = sig_bpf;
end

%% 6. SIMPAN DATA KE CSV
fprintf('\nMenyimpan hasil...\n');
output_main = 'data_imu_filtered_BPF_Powerline.csv';
output_denoised = 'data_imu_denoised_FFT.csv';
output_combined = 'data_imu_all_results.csv';

writetable(data_filtered, output_main);
writetable(data_denoised, output_denoised);

% Gabungkan data: Original, Denoised FFT, dan Filtered (BPF)
combined = table;
for i = all_sensor_cols
    colname = col_names{i};
    combined.([colname '_Original']) = table2array(data(:, i));
    combined.([colname '_FFT_Denoised']) = table2array(data_denoised(:, i));
    combined.([colname '_BPF_Filtered']) = table2array(data_filtered(:, i));
end
writetable(combined, output_combined);

fprintf('✅ Data akhir tersimpan ke:\n');
fprintf('  1. %s (hasil akhir)\n', output_main);
fprintf('  2. %s (setelah FFT Denoising)\n', output_denoised);
fprintf('  3. %s (gabungan Original, FFT, dan BPF)\n', output_combined);

%% 7. VISUALISASI
fprintf('\nMenampilkan hasil...\n');
t = (0:height(data)-1) / fs;
sample_col = all_sensor_cols(1);
original = table2array(data(:, sample_col));
denoised = table2array(data_denoised(:, sample_col));
filtered = table2array(data_filtered(:, sample_col));

% Time Domain Comparison
figure('Name', 'Time Domain Comparison', 'Position', [100 100 1200 600]);
subplot(3,1,1); plot(t, original, 'b'); grid on;
title(sprintf('Original - %s', col_names{sample_col}));
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2); plot(t, denoised, 'y'); grid on;
title('After FFT Denoising'); xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,3); plot(t, filtered, 'r'); grid on;
title('After Bandpass Filter'); xlabel('Time (s)'); ylabel('Amplitude');

% Frequency Domain
figure('Name', 'Frequency Domain', 'Position', [200 200 1200 500]);
[pxx_o, f] = pwelch(original, hamming(512), 256, 2048, fs);
[pxx_f, ~] = pwelch(filtered, hamming(512), 256, 2048, fs);
plot(f, 10*log10(pxx_o), 'b', f, 10*log10(pxx_f), 'r', 'LineWidth', 1.2);
grid on; legend('Original', 'Filtered');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density Comparison');
xline(powerline_freq, '--k', sprintf('%d Hz', powerline_freq));
xline(2*powerline_freq, '--k', sprintf('%d Hz', 2*powerline_freq));
xline(3*powerline_freq, '--k', sprintf('%d Hz', 3*powerline_freq));

% Bandpass Filter Response
figure('Name', 'Bandpass Filter Response', 'Position', [300 300 800 500]);
subplot(2,1,1);
plot(w_bpf, 20*log10(abs(h_bpf)), 'b', 'LineWidth', 1.5); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Bandpass Magnitude Response');
xline(bpf_low, '--r', 'Low Cutoff');
xline(bpf_high, '--r', 'High Cutoff');

subplot(2,1,2);
plot(w_bpf, unwrap(angle(h_bpf))*180/pi, 'r', 'LineWidth', 1.5); grid on;
xlabel('Frequency (Hz)'); ylabel('Phase (deg)');
title('Bandpass Phase Response');

%% 8. ANALISIS NOISE REDUCTION
fprintf('\n=== Analisis Noise Reduction ===\n');
for i = all_sensor_cols
    orig = table2array(data(:, i));
    filt = table2array(data_filtered(:, i));
    noise = orig - filt;

    SNR_improve = 10*log10(var(filt)/var(noise));
    fprintf('%s: SNR Improvement = %.2f dB | RMS noise = %.4f\n', ...
        col_names{i}, SNR_improve, rms(noise));
end

fprintf('\n=== Proses selesai! ===\n');
