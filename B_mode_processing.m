% DR/EP & DSC
clc; clear all;close all;

% SPEC

    f_sample    = 40e6;
    f_center    = 2.5e6;
    scanline    = 96;
    data_length = 10240;     % depth 100mm : *5
                                %       120mm : *6
                                %       160mm : *8
    
   
    CSPEED = 1540000;                  % [mm]
    lambda = CSPEED/f_center/1000;     % Wavelength    
    ele_pitch   = 2*0.26e-3;	   % Element Pitch [m]   

    depth       = 20;        % [cm]
    VIEW_ANGLE  = 80; 
    
    DynamicRange = 70;      % [db]
    view_width = scanline*ele_pitch*1000;    
    
% data loading
%     fid = fopen('Beamformed_data\sum_out_JH.bin','rb');
    fid = fopen('Beamformed_data\sum_out.bin','rb');
    data    = fread(fid,'double');
    fclose all;    
    
    data = data./max(abs(data))*2^15;
    
    data = reshape(data, length(data)/scanline,scanline);    
 
% DC cancel
    disp(sprintf('\n\tDc cancel'));
    h = fir1(128,[0.1],'high');
    dcc_out = (convn(data', h, 'same'))';
    
% tgc
    gain = ones(1,data_length);
    gain(1:1150) = (4:-3/1149:1);
    gain(1600:data_length) = (1:6/(data_length-1600):7);


    gain = (ones(scanline,1) * gain)';    
    tgc_out = dcc_out .* gain;    
%     tgc_out = dcc_out;        

%     figure(1)
%     subplot(2,1,1)
%     plot(dcc_out)
%     subplot(2,1,2)
%     plot(tgc_out)

% filter for harmonic
%     h1 = fir1(64,0.9,'low');
%     h2 = fir1(64,0.7,'high');
% %     freqz(h2)
%     tgc_out = (convn(tgc_out', h1, 'same'))';
%     tgc_out = (convn(tgc_out', h2, 'same'))';    
% 

% magnitude calculation and decimation 
    disp(sprintf('\n\tMagnitude calculation'));
    
%     f_start = 2.95e6;
%     f_end = 3.05e6;
%     f_start = 12.4e6;
%     f_end = 12.3e6;
% 
%     cos_t = cos((1:data_length)*2*pi/f_sample .* (f_start:-(f_start-f_end)/(data_length-1):f_end));
%     sin_t = sin((1:data_length)*2*pi/f_sample .* (f_start:-(f_start-f_end)/(data_length-1):f_end));    
%     
    cos_t = cos((1:data_length)*2*pi*f_center/f_sample);
    sin_t = sin((1:data_length)*2*pi*f_center/f_sample);   


for sn=1:scanline
            d_i(1:data_length,sn) = tgc_out(:,sn) .* sin_t';
            d_q(1:data_length,sn) = tgc_out(:,sn) .* cos_t';
    end        
    
    h = fir1(64,f_center/f_sample);
    
    H = abs(fft(h));
    H = 20*log10(H./max(H));    
    
    DQ = abs(fft(d_i(:,floor(scanline/2))));
    DQ = 20*log10(DQ./max(DQ));
    
    d_i_f = (convn(d_i', h, 'same'))';
    d_q_f = (convn(d_q', h, 'same'))';
    
    mag = sqrt(d_i_f.^2+d_q_f.^2);

%     subplot(2,1,1)    
%     plot(mag(:,95),'r');
%     
% log compression

    logDB = DynamicRange;

    Ymax = 2^8;
    
%    Xmax = 2^10;
    Xmax = max(max(mag));
    Xmin = 10^(-logDB/20)*Xmax;

    % mag = (mag > Xmax).*Xmax + (mag <= Xmax).*mag;
    logout = (mag >= Xmin).* (Ymax/(log10(Xmax/Xmin))*log10(mag/Xmin));
    logout(1,1)=255;  
    
%% Scan converter
    dsc_in = logout;
    disp(sprintf('\n\tScan converter\n'));
    size_out=size(dsc_in);
    
    fid=fopen('sc_interp_coef.dat','r');
    SC_Coef=fscanf(fid,'%f',[4 32]);
    fclose all;
    
    % function [sc_data]=scan_converter(input_buffer,register_file_name,...
    %                             array_type,interp_mode,scanline_num,sample_number,...
    %                             image_x,image_y,...
    %                             depth_st,depth_end,convex_radius,angle,c,...
    %                             fs,linear_array_length,hi_bound);
 
    [sc_data]=scan_converter_rev1(dsc_in,  SC_Coef', 1, 5,  size_out(2),    size_out(1),...
                            640 , 480, 0, depth*10,    0,  VIEW_ANGLE*2*pi/360,    CSPEED,...
                             f_sample,  view_width,  2^8);
                        
                         
    figure(1004)
    
    map = [0:255; 0:255; 0:255]'/255;
    sc_data(1,1)=255;
    x_label = (0:scanline-1)*ele_pitch*1000;    
    imagesc(x_label,0:depth*10,sc_data),colormap(map),truesize;


        F = getframe(gca);
%         imwrite(F.cdata,'image_temp.bmp','bmp');