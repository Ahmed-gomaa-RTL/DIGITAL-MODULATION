clear
clc

% Number of bits to generate
Num_bits = 100000;
Bits = [0, 1];

% Generate random bit stream Tb = 1sec
time_NRZ = 0:1:Num_bits;
time_RZ = 0:0.5:Num_bits;
time_Ary = 0:1:Num_bits;
Stream_bits = datasample(Bits, Num_bits);

%Values of Eb/No
EbN0_dB = -10 : 2 : 6;

% Polar NRZ PCM Encoding
Encode_NRZ = Stream_bits * 2 - 1;
Encode_NRZ_plot = Encode_NRZ;
Encode_NRZ_plot(length(Stream_bits) + 1) = -1;

% Polar RZ PCM Encoding
Encode_RZ_plot = zeros(1, Num_bits * 2 + 1);
for i = 1:1:Num_bits
    Encode_RZ_plot(2 * i - 1) = Encode_NRZ(i);
end
Encode_RZ = Encode_RZ_plot(1:end - 1);

% 4-Ary PCM Encoding
Four_Ary_plot = zeros(1, Num_bits);
for j = 2:2:Num_bits
    if (Stream_bits(j - 1) == 0 && Stream_bits(j) == 0)
        Four_Ary_plot(j) = -3;
         Four_Ary_plot(j - 1) = -3;
    elseif (Stream_bits(j - 1) == 0 && Stream_bits(j) == 1)
        Four_Ary_plot(j) = -1;
        Four_Ary_plot(j - 1) = -1;
    elseif (Stream_bits(j - 1) == 1 && Stream_bits(j) == 1)
        Four_Ary_plot(j) = 1;
        Four_Ary_plot(j - 1) = 1;
    elseif (Stream_bits(j - 1) == 1 && Stream_bits(j) == 0)
        Four_Ary_plot(j) = 3;
        Four_Ary_plot(j - 1) = 3;
    end
end
Encode_4_Ary = Four_Ary_plot;
Four_Ary_plot(end + 1) = -3;

% Display generated and encoded information
fprintf('Random bit stream of %d bits generated and encoded using:\n', Num_bits);
fprintf('1. Polar NRZ PCM\n2. Polar RZ PCM\n3. 4-Ary PCM with Gray Coding\n');

% Plotting Section
% Plot Polar NRZ PCM
figure;
subplot ( 3 , 1 , 1 ) ;
stairs(time_NRZ, Encode_NRZ_plot, 'b', 'LineWidth', 2.5);
grid on;
box on;
title('Polar NRZ PCM', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([-1.5 1.5]);
xlim([0 10]);

% Plot Polar RZ PCM
subplot ( 3 , 1 , 2 ) ;
stairs(time_RZ, Encode_RZ_plot, 'g', 'LineWidth', 2.5);
grid on;
box on;
title('Polar RZ PCM', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([-1.5 1.5]);
xlim([0 10]);

% Plot 4-Ary PCM
subplot ( 3 , 1 , 3 ) ;
stairs(time_Ary, Four_Ary_plot, 'r', 'LineWidth', 2.5);
grid on;
box on;
title('4-Ary PCM', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([-3.5 3.5]);
xlim([0 10]);

% Plot Constellation
% Plot Polar NRZ PCM
figure;
subplot( 3 , 1 , 1 );
plot(Encode_NRZ, zeros(size(Encode_NRZ)), '.', 'LineWidth', 100);
title('Polar NRZ PCM Constallation', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
box on;

% Plot Polar RZ PCM
subplot( 3 , 1 , 2 );
plot(Encode_RZ(1:2:end) / sqrt(2), zeros(size(Encode_RZ(1:2:end))), '.', 'LineWidth', 100);
title('Polar RZ PCM Constallation', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
box on;

% Plot Polar 4-Ary PCM
subplot( 3 , 1 , 3 );
plot(Encode_4_Ary(1:2:end) * sqrt(2), zeros(size(Encode_4_Ary(1:2:end))), '.', 'LineWidth', 100);
title('Polar 4 Ary PCM Constallation', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
grid on;
box on;

% Probabilty of 0 and 1
Count0 = sum(Stream_bits == 0);
Probability0 = Count0 / length(Stream_bits);
Probability1 = 1 - Probability0;

% Transmit through AWGN Channel
for EbN0_dB_f = -10 : 2 : 6
    % Calculate noise power NRZ 
    EbN0_Power = 10^(EbN0_dB_f/10);
    noise_power_NRZ = 1 / (EbN0_Power);
    
    % Calculate noise power RZ  
    noise_power_RZ = 1 / (EbN0_Power);
    
    % Calculate noise power 4-Ary 
    noise_power_4Ary = 1 / (EbN0_Power);

    % Add AWGN to NRZ PCM
    noise_NRZ = sqrt(noise_power_NRZ / 2) * randn(size(Encode_NRZ));
    received_NRZ = Encode_NRZ + noise_NRZ;
    received_NRZ_plot = received_NRZ;
    received_NRZ_plot(end + 1) = -1;

    % Add AWGN to RZ PCM
    noise_RZ = sqrt(noise_power_RZ) * randn(size(Encode_RZ));
    received_RZ = Encode_RZ + noise_RZ;
    received_RZ_plot = received_RZ;
    received_RZ_plot(end + 1) = 0;

    % Add AWGN to 4-Ary PCM
    noise_Ary = sqrt(noise_power_4Ary / 4) * randn(size(Encode_4_Ary));
    received_Ary = Encode_4_Ary + noise_Ary;
    received_Ary_plot = received_Ary;
    received_Ary_plot(end + 1) = -3;
    
    % Optimum Threshold NRZ
    Threshold_NRZ = (noise_power_NRZ / 4) * log10(Probability0 / Probability1);

    % Plot diagrams for received signals
    figure;
    subplot(3, 1, 1);
    stairs(time_NRZ, received_NRZ_plot, 'b', 'LineWidth', 2.5);
    hold on;
    plot(time_NRZ, Threshold_NRZ * ones(size(time_NRZ)));
    title(['Received NRZ PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlim([0 10]);
    grid on;
    box on;

    subplot(3, 1, 2);
    stairs(time_RZ, received_RZ_plot, 'g', 'LineWidth', 2.5);
    hold on;
    plot(time_RZ, zeros(size(time_RZ)));
    title(['Received RZ PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlim([0 10]);
    grid on;
    box on;

    subplot(3, 1, 3);
    stairs(time_Ary, received_Ary_plot, 'r', 'LineWidth', 1.5);
    hold on;
    plot(time_Ary, zeros(size(time_Ary)));
    plot(time_Ary, 2 * ones(size(time_Ary)));
    plot(time_Ary, -2 * ones(size(time_Ary)), 'Color', [1, 0.41, 0.16]);
    title(['Received 4-Ary PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Amplitude', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Time', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlim([0 10]);
    grid on;
    box on;
    
    % Plot Constellation
    % Plot Polar NRZ PCM
    figure;
    subplot( 3 , 1 , 1 );
    plot(received_NRZ, zeros(size(received_NRZ)), '.', 'LineWidth', 100);
    title(['Received Constellation NRZ PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    box on;

    % Plot Polar RZ PCM
    subplot( 3 , 1 , 2 );
    plot(received_RZ(1:2:end) / sqrt(2), zeros(size(received_RZ(1:2:end))), '.', 'LineWidth', 100);
    title(['Received Constellation RZ PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    box on;

    % Plot Polar 4-Ary PCM
    subplot( 3 , 1 , 3 );
    plot(received_Ary(1:2:end) * sqrt(2), zeros(size(received_Ary(1:2:end))), '.', 'LineWidth', 100);
    title(['Received Constellation 4 Ary PCM (Eb/N0 = ', num2str(EbN0_dB_f), ' dB)'], 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    ylabel('Imagine', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel('Real', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    grid on;
    box on;
    
    % Create the matched filter
    % NRZ PCM
    matched_filter_NRZ = fliplr(Encode_NRZ);
    % RZ PCM
    matched_filter_RZ = fliplr(Encode_RZ);
    % 4-Ary PCM
    matched_filter_ARY = fliplr(Encode_4_Ary);

    % Apply the matched filter
    % NRZ PCM
    output_filter_NRZ = conv(received_NRZ, matched_filter_NRZ, 'same');
    % RZ PCM
    output_filter_RZ = conv(received_RZ, matched_filter_RZ, 'same');
    % 4-Ary PCM
    output_filter_ARY = conv(received_Ary, matched_filter_ARY, 'same');
    
    % Signal after threshold NRZ
    received_NRZ_threshold = (received_NRZ > Threshold_NRZ) * 1 + (received_NRZ < Threshold_NRZ) * -1;
    
    % Signal after threshold RZ threshold = 0
    received_RZ_threshold = (received_RZ(1:2:end) > 0) * 1 + (received_RZ(1:2:end) < 0) * -1;
    
    % Signal after threshold 4-Ary threshold = 0, 2 , 2
    received_Ary_threshold = (received_Ary > 0) .* ((received_Ary > 2) * 3 + (received_Ary < 2) * 1) + (received_Ary < 0) .* ((received_Ary > -2) * -1 + (received_Ary < -2) * -3);
    
    %Bit Error Rate NRZ
    BER_NRZ = sum (received_NRZ_threshold ~= Encode_NRZ) / length(received_NRZ_threshold);
    BER_VS_No_NRZ(6 - (EbN0_dB_f / -2)) = BER_NRZ;
    
    %Bit Error Rate RZ
    BER_RZ = sum (received_RZ_threshold ~= Encode_RZ(1:2:end)) / length(received_RZ_threshold);
    BER_VS_No_RZ(6 - (EbN0_dB_f / -2)) = BER_RZ;
    
    %Bit Error Rate 4-Ary
    BER_Ary = sum (received_Ary_threshold ~= Encode_4_Ary) / length(received_Ary_threshold);
    BER_VS_No_Ary(6 - (EbN0_dB_f / -2)) = BER_Ary;
    
    %Bit Error Rate Theoritical NRZ
    BER_Theoritical = (Probability0 * 0.5 * erfc((1 + Threshold_NRZ) / ((EbN0_Power)^-.5))) + (Probability1 * 0.5 * erfc((1 - Threshold_NRZ) / ((EbN0_Power)^-0.5)));
    BER_Theoritical_VS_No(6 - (EbN0_dB_f / -2)) = BER_Theoritical;
      
end    

% Plot BER VS Eb/no and BER theoritical
figure;
semilogy(EbN0_dB, BER_VS_No_NRZ, 'LineWidth', 1);
hold on 
semilogy(EbN0_dB, BER_Theoritical_VS_No, 'LineWidth', 1);
hold on
semilogy(EbN0_dB, BER_VS_No_RZ, 'LineWidth', 1);
hold on
semilogy(EbN0_dB, BER_VS_No_Ary, 'LineWidth', 1);
title('BER VS Eb/No ', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylabel('BER', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
xlabel('Eb/No', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
legend('BER VS Eb/No NRZ', 'BER Theoritical VS Eb/No NRZ', 'BER VS Eb/No RZ', 'BER VS Eb/No 4-Ary', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
ylim([10^-5  1]); 
xlim([-10 6]); 
grid on;
box on;
