% view_single_gyro.m
% Script untuk menampilkan isi dan grafik dari satu file .mat
% dan menyimpan hasilnya ke CSV
clc; clear; close all;

% --- 1️⃣ Tentukan lokasi file yang ingin kamu buka
filePath = '/Users/mbprom2512/merge_EMG_Gyro/S1P1_G.mat';  % ubah sesuai lokasi file kamu

% --- 2️⃣ Muat isi file
S = load(filePath);
disp('✅ File berhasil dimuat.');

% --- 3️⃣ Tampilkan daftar variabel di dalam file
disp('📦 Variabel yang ada di dalam file:');
disp(fieldnames(S));

% --- 4️⃣ Cek apakah variabel gyro ada
if isfield(S, 'gyroX') && isfield(S, 'gyroY') && isfield(S, 'gyroZ') && isfield(S, 'timestamps')
    % Ambil data
    gyroX = S.gyroX;
    gyroY = S.gyroY;
    gyroZ = S.gyroZ;
    timestamps = S.timestamps;
    
    % --- 5️⃣ Gabungkan ke dalam tabel
    data = table(timestamps, gyroX, gyroY, gyroZ);
    openvar('data');   % buka di Variable Editor (seperti Excel)
    
    % --- 6️⃣ Simpan ke file CSV
    % Buat nama file CSV otomatis berdasarkan nama file .mat
    [filePath_dir, fileName, ~] = fileparts(filePath);
    csvFileName = fullfile(filePath_dir, [fileName, '.csv']);
    
    % Tulis tabel ke CSV
    writetable(data, csvFileName);
    disp(['💾 Data berhasil disimpan ke: ', csvFileName]);
    
    % --- 7️⃣ Plot ketiga sumbu gyro terhadap waktu
    figure('Name', 'Gyroscope Data (Time Series)', 'NumberTitle', 'off');
    plot(timestamps, gyroX, 'r', 'DisplayName','gyroX'); hold on;
    plot(timestamps, gyroY, 'g', 'DisplayName','gyroY');
    plot(timestamps, gyroZ, 'b', 'DisplayName','gyroZ');
    xlabel('Time (s)');
    ylabel('Gyro Value');
    title('Gyroscope Data (X, Y, Z)');
    legend;
    grid on;
else
    disp('⚠️ File ini tidak memiliki variabel gyroX, gyroY, gyroZ, dan timestamps.');
end