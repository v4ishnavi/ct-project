clc;clear;close all;
%A/D Conversion
[audio1, fs] = audioread('project.wav');

% This is the A/D converter, completely copied honestly

audio2 = audio1(:,2);
audio2 = audio2(1:5000);
% disp(length(audio));

bl = 1000000; %block length
n = ceil(length(audio2)/bl);
% n = 5;
%number of fullsize block lengths
% bl = 100;
% n = 1;
% op_rect_ml = zeros(1,16*n*bl); %16 for 16 binary bits?
% op_rcos_ml = zeros(1, 16*n*bl);
% op_rect_m_ = zeros(1, 16*n*bl);
% op_rcos_m_ = zeros(1, 16*n*bl);

% op_rect_ml = zeros(1, 16*length(audio2)); %16 for 16 binary bits?
% op_rcos_ml = zeros(1, 16*length(audio2));
% op_rect_m_ = zeros(1, 16*length(audio2));
% op_rcos_m_ = zeros(1, 16*length(audio2));

op_rect_ml = '';
op_rcos_ml = '';
op_rect_m_ = '';
op_rcos_m_ = '';

for i = 1:n
    
    if i == n
        infocus = 1+(i-1)*bl:length(audio2);
        outfocus = 1+(i-1)*16*bl:16*length(audio2);
    else
        infocus = 1+(i-1)*bl:bl+(i-1)*bl;
        outfocus = 1+(i-1)*16*bl:16*bl+(i-1)*16*bl;         
    end
    disp(['Length of a segment: ',num2str(length(outfocus))]);
    audio = audio2(infocus);
    
    audio_normalized = int16(audio * 32767); 
    audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16); 
    binary_vector = audio_binary(:)';
    % disp(binary_vector);
    
    %------------------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             encoding                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    l = length(binary_vector);
    k_rcos = 8; 
    k_rect = 6;
    m = 4;
    
    %This section of the code encodes the binary vector using 4-PAM 
    symbols_vector = fourpammapA(binary_vector);
    
    %-------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             line encoding                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp(l);
    symbols_vector_up = zeros(1,m*l/2);
    for j = 1:l/2
        symbols_vector_up(j + (j-1)*(m-1)) = symbols_vector(j);
    end 
    % Rectangular pulse
    p1 = [ones(1,k_rect/2)];
    %raised cosine pulse 
    a_ = 0.5;
    [p2, ~] = raised_cosine(a_, 4, k_rcos/2);
    
    line_vector_rect = conv(symbols_vector_up,p1,'same');
    line_vector_rcos = conv(symbols_vector_up, p2, 'same');

    % break;

    %-------------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             modulation                              %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %rectangular encoded: 
     fc = 1e6;
     Tb = 2/fs;
    len_ip = length(symbols_vector_up); 
    %Tb = Ts*logM
    t = (0:len_ip-1)/(0.5*fs);
    carrier = cos(fc*2*pi*t);

    modulated_rect_vec = line_vector_rect .* carrier;
    modulated_rcos_vec = line_vector_rcos .* carrier;


    
    %_---------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                                 channel                             %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    snr_val = 20;
    % Channel memoryless: 
    channel_op_rcos = awgn(modulated_rcos_vec, snr_val);
    channel_op_rect = awgn(modulated_rect_vec, snr_val); 

    %channel memory: 
    % Parameters
    a = 0.9; 
   
    b = 10;
    t_h = 0:Tb:(10*Tb);  % Covering 0 to 2Tb for two delta functions
    h = a * (t_h == 0) + (1 - a) * (t_h == b*Tb);  % Impulse response
    
    % ch_m_o_rcos = conv(modulated_rcos_vec, h, 'same'); 
    % chmo_rcos = awgn(ch_m_o_rcos, snr_val);
    
    ch_m_o_rect = conv(modulated_rect_vec, h, 'same');
    chmo_rect = awgn(ch_m_o_rect, snr_val);
    chmo_rect = linear_equalizer(chmo_rect, modulated_rect_vec, 0.0001, 64);
    
    ch_m_o_rcos = conv(modulated_rcos_vec, h, 'same');
    chmo_rcos = awgn(ch_m_o_rcos, snr_val);
    chmo_rcos = linear_equalizer(chmo_rcos, modulated_rcos_vec, 0.0001, 64);


    % disp(chmo_rcos(1:100));
    % figure;
    % plot(chmo_rcos);
    % hold on
    % plot(modulated_rcos_vec);
    % hold off

    %---------------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             demodulation                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % demod_rect_vec = lowpass(channel_op_rect.*carrier, fc, 2*fc);
    demod_rect_vec = lowpass(channel_op_rect.*carrier,fc,2*fc+1);
    demod_rcos_vec = lowpass(channel_op_rcos.*carrier, fc,2*fc+1); 

  
    % figure;
    % plot(t, demod_rect_vec);
    % plot(t, demod_rcos_vec);
    
    demod_rect_vec_mem = lowpass(chmo_rect.*carrier, fc, 2*fc+1);
    demod_rcos_vec_mem = lowpass(chmo_rcos.*carrier, fc, 2*fc+1);
    %  figure;
    % plot((200:400), demod_rect_vec_mem(200:400), 'LineWidth', 1.5);
    % title("Demodulated output of rectangular pulse")
    % xlabel("time")
    % ylabel("amplitudes")
    % % ylim([-4 4])

    %     figure;
    %   plot((200:400), demod_rect_vec_mem(200:400), 'LineWidth', 1.5);
    % title("Demodulated output of rectangular pulse")
    % xlabel("time")
    % ylabel("amplitudes")
    % break;
    %----------------------------------------------------------------------------------
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                             line decoding                           %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    line_decoded_rect_vec = conv(demod_rect_vec, p1, 'same');
    % line_decoded_rect_vec = 3*line_decoded_rect_vec/max(line_decoded_rect_vec);
    
    line_decoded_rcos_vec = conv(demod_rcos_vec, p2, 'same');
    % line_decoded_rcos_vec = 3*line_decoded_rcos_vec/max(line_decoded_rcos_vec);
    % figure;
    % plot(line_decoded_rect_vec)
    % plot(line_decoded_rcos_vec);
    
    line_decoded_rect_vec_mem = conv(demod_rect_vec_mem,p1,'same');
    line_decoded_rcos_vec_mem = conv(demod_rcos_vec_mem,p2,'same');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                            4PAM decoding                            %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    decoded_rect = (line_decoded_rect_vec);
    decoded_downsample_rect = downsample(decoded_rect, m);
    
    % figure;
    % stem(real(decoded_downsample_rect),imag(decoded_downsample_rect));
    
    decoded_rcos = (line_decoded_rcos_vec);
    decoded_downsample_rcos = downsample(decoded_rcos, m);
    
    
    % figure;
    % stem(real(decoded_downsample_rcos),imag(decoded_downsample_rcos));
    
    
    decoded_rect_mem = (line_decoded_rect_vec_mem);
    decoded_downsample_rect_mem = downsample(decoded_rect_mem,m);
    
    % 
    % figure;
    % stem(real(decoded_downsample_rect_mem),imag(decoded_downsample_rect_mem));
    
    decoded_rcos_mem = (line_decoded_rcos_vec_mem);
    decoded_downsample_rcos_mem = downsample(decoded_rcos_mem,m);
    
    % figure;
    % stem(real(decoded_downsample_rcos_mem),imag(decoded_downsample_rcos_mem));
    
    final_output_rect = fourpamunmapA(decoded_downsample_rect,1);
    final_output_rcos = fourpamunmapA(decoded_downsample_rcos,0);
    
    final_output_rect_mem = fourpamunmapA(decoded_downsample_rect_mem,0);
    final_output_rcos_mem = fourpamunmapA(decoded_downsample_rcos_mem,0);
    
    count_rect = 0; 
    count_rcos = 0; 
    count_rect_mem = 0;
    count_rcos_mem = 0;
    % disp(length(binary_vector)); disp(length(final_output_rcos))
    
    % final_output_rect = binary_vector;
    
    for p = 1: length(binary_vector)
 
        if binary_vector(p) ~= final_output_rect(p)

            count_rect = count_rect + 1;
        end
        if binary_vector(p) ~= final_output_rcos(p)
            count_rcos = count_rcos + 1;
        end
        if binary_vector(p) ~= final_output_rect_mem(p)
            count_rect_mem = count_rect_mem + 1;
        end
        if binary_vector(p) ~= final_output_rcos_mem(p)
            count_rcos_mem = count_rcos_mem + 1;
        end
    end
    disp('-------------------------------------------------------')
    disp("BER OF THE FOLLOWING:")
    disp(['RCOS WITH MEMORY: ',num2str(count_rcos_mem/length(binary_vector))])
    disp(['RECT WITH MEMORY: ', num2str(count_rect_mem/length(binary_vector))])
    disp(['RCOS MEMORYLESS: ', num2str(count_rcos/length(binary_vector))]);
    disp(['RECT MEMORYLESS: ', num2str(count_rect/length(binary_vector))]);
    disp('--------------------------------------------------------')

    % op_rect_ml(outfocus) = final_output_rect;
    % op_rcos_ml(outfocus) = final_output_rcos;
    % op_rect_m_(outfocus) = final_output_rect_mem;
    % op_rcos_m_(outfocus) = final_output_rcos_mem;

    op_rect_ml(i,1:length(binary_vector)) = final_output_rect;
    op_rcos_ml(i,1:length(binary_vector)) = final_output_rcos;
    op_rect_m_(i,1:length(binary_vector)) = final_output_rect_mem;
    op_rcos_m_(i,1:length(binary_vector)) = final_output_rcos_mem;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                            D/A conversion                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sound(double(typecast(uint16(bin2dec(num2str(reshape(op_rcos_m_, [], 16)))), 'int16'))/32767, fs)

