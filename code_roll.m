% ANALISIS GYRO Z (Rotasi Roll) DENGAN RMS
try
    % Load file .mat
    load('S1P1_G.mat');
    disp('File S1P1_G.mat berhasil dimuat.');
    
    % --- Ambil data ---
    gyro_z = gyroZ;           % Rotasi roll (sumbu Z)
    t = timestamps;           % Waktu
    
    % Pastikan waktu dalam detik dan mulai dari nol
    t = t - t(1);
    if max(t) > 1e5
        t = t / 1000;  % jika timestamp dalam ms, ubah ke detik
    end
    
    % Hitung frekuensi sampling
    fs = 1 / mean(diff(t));
    
    % --- Hitung RMS dengan jendela 1 detik ---
    window_time = 1;                       % jendela waktu (detik)
    window_size = round(window_time * fs); % ubah ke sampel
    rms_values = sqrt(movmean(gyro_z.^2, window_size)); % RMS moving window
    
    % --- Hitung nilai rata-rata RMS ---
    mean_rms = mean(rms_values);
    
    % --- Plot hasil ---
    figure;
    subplot(2,1,1);
    plot(t, gyro_z, 'b');
    xlabel('Waktu (s)');
    ylabel('Gyro Z (°/s)');
    title('Sinyal Gyro Z (Rotasi Roll)');
    grid on;
    
    subplot(2,1,2);
    plot(t, rms_values, 'r', 'LineWidth', 1.5); hold on;
    yline(mean_rms, '--w', sprintf('Mean RMS = %.3f', mean_rms), ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
    xlabel('Waktu (s)');
    ylabel('Nilai RMS');
    title('Nilai RMS dari Gyro Z terhadap Waktu');
    legend('RMS(t)', 'Mean RMS', 'Location', 'best');
    grid on;
    
    disp(['Frekuensi sampling terdeteksi: ', num2str(fs), ' Hz']);
    disp('Analisis selesai.');
    
catch ME
    disp('Terjadi kesalahan:');
    disp(ME.message);
end