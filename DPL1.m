% ========================================================================
% IMU Power Line Noise Removal (Tanpa iirnotch)
% Menggunakan metode: Bandstop Filter, FFT, atau Adaptive Filter
% ========================================================================

clear; clc; close all;

%% 1. KONFIGURASI PARAMETER
% Ganti dengan nama file CSV Anda
filename = '.csv';

% Parameter filter
fs = 100;              % Sampling frequency (Hz) - sesuaikan dengan data Anda
powerline_freq = 50;   % Frekuensi power line (50 Hz untuk Indonesia/Eropa, 60 Hz untuk US)
notch_bandwidth = 2;   % Bandwidth notch filter (Hz)

% Pilih metode filtering (1, 2, atau 3)
filter_method = 1;     % 1: Butterworth Bandstop, 2: FFT Domain, 3: Moving Average

fprintf('=== METODE FILTERING ===\n');
fprintf('1: Butterworth Bandstop Filter\n');
fprintf('2: FFT Domain Filtering\n');
fprintf('3: Moving Average Filter\n');
fprintf('Metode dipilih: %d\n\n', filter_method);

%% 2. BACA DATA CSV
fprintf('Membaca data dari %s...\n', filename);
try
    data = readtable(filename);
    fprintf('Data berhasil dibaca!\n');
    fprintf('Jumlah kolom: %d\n', width(data));
    fprintf('Jumlah baris: %d\n', height(data));
    fprintf('Nama kolom: ');
    disp(data.Properties.VariableNames);
catch
    error('Gagal membaca file CSV. Pastikan file ada dan formatnya benar.');
end

%% 3. IDENTIFIKASI KOLOM DATA IMU
col_names = data.Properties.VariableNames;
accel_cols = [];
gyro_cols = [];

% Cari kolom accelerometer dan gyroscope
for i = 1:length(col_names)
    if contains(lower(col_names{i}), {'acc', 'accel'})
        accel_cols = [accel_cols, i];
    elseif contains(lower(col_names{i}), {'gyro', 'gyr'})
        gyro_cols = [gyro_cols, i];
    end
end

fprintf('\nKolom Accelerometer terdeteksi: ');
disp(col_names(accel_cols));
fprintf('Kolom Gyroscope terdeteksi: ');
disp(col_names(gyro_cols));

%% 4. DESAIN FILTER BERDASARKAN METODE
switch filter_method
    case 1  % === BUTTERWORTH BANDSTOP FILTER ===
        fprintf('\n=== Menggunakan Butterworth Bandstop Filter ===\n');
        
        nyquist_freq = fs/2;
        
        % Design bandstop filter untuk frekuensi fundamental
        f_low = powerline_freq - notch_bandwidth/2;
        f_high = powerline_freq + notch_bandwidth/2;
        
        % Validasi frekuensi
        if f_low > 0 && f_high < nyquist_freq
            [b1, a1] = butter(4, [f_low f_high]/nyquist_freq, 'stop');
            use_filter1 = true;
        else
            use_filter1 = false;
            fprintf('  WARNING: Fundamental frequency out of range, skipped\n');
        end
        
        % Design untuk harmonik ke-2
        f_low2 = 2*powerline_freq - notch_bandwidth/2;
        f_high2 = 2*powerline_freq + notch_bandwidth/2;
        
        if f_low2 > 0 && f_high2 < nyquist_freq
            [b2, a2] = butter(4, [f_low2 f_high2]/nyquist_freq, 'stop');
            use_filter2 = true;
        else
            use_filter2 = false;
            fprintf('  WARNING: 2nd harmonic out of range, skipped\n');
        end
        
        % Design untuk harmonik ke-3
        f_low3 = 3*powerline_freq - notch_bandwidth/2;
        f_high3 = 3*powerline_freq + notch_bandwidth/2;
        
        if f_low3 > 0 && f_high3 < nyquist_freq
            [b3, a3] = butter(4, [f_low3 f_high3]/nyquist_freq, 'stop');
            use_filter3 = true;
        else
            use_filter3 = false;
            fprintf('  WARNING: 3rd harmonic out of range, skipped\n');
        end
        
        fprintf('Bandstop ranges:\n');
        if use_filter1
            fprintf('  Fundamental: %.1f - %.1f Hz\n', f_low, f_high);
        end
        if use_filter2
            fprintf('  Harmonik 2: %.1f - %.1f Hz\n', f_low2, f_high2);
        end
        if use_filter3
            fprintf('  Harmonik 3: %.1f - %.1f Hz\n', f_low3, f_high3);
        end
        
    case 2  % === FFT DOMAIN FILTERING ===
        fprintf('\n=== Menggunakan FFT Domain Filtering ===\n');
        fprintf('Frekuensi yang akan dihilangkan:\n');
        fprintf('  %d Hz (fundamental)\n', powerline_freq);
        fprintf('  %d Hz (harmonik 2)\n', 2*powerline_freq);
        fprintf('  %d Hz (harmonik 3)\n', 3*powerline_freq);
        
    case 3  % === MOVING AVERAGE FILTER ===
        fprintf('\n=== Menggunakan Moving Average Filter ===\n');
        window_size = round(fs / powerline_freq);
        fprintf('Window size: %d samples\n', window_size);
end

%% 5. APLIKASIKAN FILTER KE DATA
data_filtered = data;
all_sensor_cols = [accel_cols, gyro_cols];

fprintf('\nMemfilter data sensor...\n');

for i = all_sensor_cols
    original_signal = table2array(data(:, i));
    filtered_signal = original_signal;
    
    switch filter_method
        case 1  % Butterworth Bandstop
            % Aplikasikan cascaded bandstop filter
            if use_filter1
                filtered_signal = filtfilt(b1, a1, filtered_signal);  % Fundamental
            end
            if use_filter2
                filtered_signal = filtfilt(b2, a2, filtered_signal);  % Harmonik 2
            end
            if use_filter3
                filtered_signal = filtfilt(b3, a3, filtered_signal);  % Harmonik 3
            end
            
        case 2  % FFT Domain
            % Transform ke frequency domain
            N = length(original_signal);
            X = fft(original_signal);
            freq = (0:N-1) * fs / N;
            
            % Buat notch mask untuk menghilangkan frekuensi power line
            notch_mask = ones(size(X));
            
            % Tentukan range frekuensi yang akan dihilangkan
            for harmonic = 1:3
                target_freq = harmonic * powerline_freq;
                freq_idx = abs(freq - target_freq) < notch_bandwidth;
                notch_mask(freq_idx) = 0;
                
                % Mirror untuk frekuensi negatif
                freq_idx_neg = abs(freq - (fs - target_freq)) < notch_bandwidth;
                notch_mask(freq_idx_neg) = 0;
            end
            
            % Aplikasikan mask dan transform kembali
            X_filtered = X .* notch_mask;
            filtered_signal = real(ifft(X_filtered));
            
        case 3  % Moving Average
            % Gunakan moving average untuk menghilangkan komponen periodik
            window_size = round(fs / powerline_freq);
            
            % Subtract moving average (high-pass effect untuk powerline)
            ma = movmean(original_signal, window_size);
            filtered_signal = original_signal - (ma - mean(original_signal));
            
            % Atau gunakan median filter sebagai alternatif
            % filtered_signal = original_signal - medfilt1(original_signal, window_size) + mean(original_signal);
    end
    
    data_filtered.(col_names{i}) = filtered_signal;
end

%% 6. SIMPAN DATA HASIL FILTERING
output_filename = 'data_imu_filtered.csv';
writetable(data_filtered, output_filename);
fprintf('\nData terfilter disimpan ke: %s\n', output_filename);

%% 7. VISUALISASI HASIL
fprintf('\nMembuat visualisasi...\n');

% Buat time vector
t = (0:height(data)-1) / fs;

% Pilih satu kolom untuk visualisasi detail
if ~isempty(accel_cols)
    sample_col = accel_cols(1);
elseif ~isempty(gyro_cols)
    sample_col = gyro_cols(1);
else
    error('Tidak ada kolom sensor yang terdeteksi');
end

% Figure 1: Time Domain Comparison
figure('Name', 'Power Line Denoising - Time Domain', 'Position', [100 100 1200 600]);

subplot(2,1,1);
plot(t, table2array(data(:, sample_col)), 'b', 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('Amplitude');
title(sprintf('Original Signal - %s', col_names{sample_col}));
grid on;

subplot(2,1,2);
plot(t, table2array(data_filtered(:, sample_col)), 'r', 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('Amplitude');
title(sprintf('Filtered Signal - %s', col_names{sample_col}));
grid on;

% Figure 2: Frequency Domain Analysis
figure('Name', 'Power Line Denoising - Frequency Domain', 'Position', [150 150 1200 600]);

original = table2array(data(:, sample_col));
filtered = table2array(data_filtered(:, sample_col));

% Compute power spectral density
[pxx_orig, f] = pwelch(original, hamming(512), 256, 2048, fs);
[pxx_filt, ~] = pwelch(filtered, hamming(512), 256, 2048, fs);

subplot(2,1,1);
plot(f, 10*log10(pxx_orig), 'b', 'LineWidth', 1.5);
hold on;
xline(powerline_freq, '--r', sprintf('%d Hz', powerline_freq), 'LineWidth', 2);
xline(2*powerline_freq, '--r', sprintf('%d Hz', 2*powerline_freq), 'LineWidth', 2);
xline(3*powerline_freq, '--r', sprintf('%d Hz', 3*powerline_freq), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Original Signal - PSD');
xlim([0 min(200, fs/2)]);
grid on;
legend('PSD', 'Power Line Frequencies', 'Location', 'best');

subplot(2,1,2);
plot(f, 10*log10(pxx_filt), 'r', 'LineWidth', 1.5);
hold on;
xline(powerline_freq, '--k', sprintf('%d Hz', powerline_freq), 'LineWidth', 2);
xline(2*powerline_freq, '--k', sprintf('%d Hz', 2*powerline_freq), 'LineWidth', 2);
xline(3*powerline_freq, '--k', sprintf('%d Hz', 3*powerline_freq), 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Filtered Signal - PSD');
xlim([0 min(200, fs/2)]);
grid on;
legend('PSD', 'Power Line Frequencies', 'Location', 'best');

% Figure 3: Comparison and Noise
figure('Name', 'Comparison - Original vs Filtered', 'Position', [200 200 1200 600]);

noise = original - filtered;

subplot(3,1,1);
plot(t, original, 'b', 'LineWidth', 0.8);
title('Original Signal');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(t, filtered, 'r', 'LineWidth', 0.8);
title('Filtered Signal');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, noise, 'g', 'LineWidth', 0.8);
title('Removed Noise (Original - Filtered)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% 8. ANALISIS NOISE REDUCTION
fprintf('\n=== Analisis Noise Reduction ===\n');

for i = all_sensor_cols
    original = table2array(data(:, i));
    filtered = table2array(data_filtered(:, i));
    noise = original - filtered;
    
    % Hitung SNR improvement
    signal_power = var(filtered);
    noise_power = var(noise);
    snr_improvement = 10*log10(signal_power / noise_power);
    
    % Hitung persentase noise reduction
    noise_reduction_pct = (rms(noise) / rms(original)) * 100;
    
    fprintf('\n%s:\n', col_names{i});
    fprintf('  SNR Improvement: %.2f dB\n', snr_improvement);
    fprintf('  RMS Original: %.6f\n', rms(original));
    fprintf('  RMS Filtered: %.6f\n', rms(filtered));
    fprintf('  RMS Noise Removed: %.6f\n', rms(noise));
    fprintf('  Noise Reduction: %.2f%%\n', noise_reduction_pct);
    
    % Analisis frekuensi power line
    [pxx_noise, f_noise] = pwelch(noise, [], [], [], fs);
    [~, idx_powerline] = min(abs(f_noise - powerline_freq));
    power_at_powerline = 10*log10(pxx_noise(idx_powerline));
    
    fprintf('  Power at %d Hz in noise: %.2f dB/Hz\n', powerline_freq, power_at_powerline);
end

fprintf('\n=== Proses Selesai! ===\n');
fprintf('Metode yang digunakan: ');
switch filter_method
    case 1
        fprintf('Butterworth Bandstop Filter\n');
    case 2
        fprintf('FFT Domain Filtering\n');
    case 3
        fprintf('Moving Average Filter\n');
end