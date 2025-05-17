clc; clear; close all;
format long g;

%% -------------------------------
% Parameter settings
%% -------------------------------
Beam_Num       = 96;
Elem_Pitch     = 0.26;           % [mm]
Elem_Num       = 64;
Ch_Num         = 16;
Sound_Speed    = 1540000;        % [mm/s]
Fs             = 40e6;           % [Hz]
data_length    = 10240;          % Original data length
Unit_Dis       = Sound_Speed/(2*Fs);  % Physical depth interval [mm]

% Interpolation (upsampling) parameters
inter_factor   = 4;              % Interpolation factor (enhance time resolution)
% Effective sampling frequency after interpolation
Fs_new = Fs * inter_factor;

View_Angle     = 80;             % Transmit angle range (-80 deg ~ 80 deg)
output_dir     = 'Beamformed_data';

% Element x-coordinates (centered)
x_element = ((0:Elem_Num-1) - (Elem_Num-1)/2) * Elem_Pitch;

%% -------------------------------
% Initialize output array (data_length x Beam_Num)
%% -------------------------------
Sum_out_2D = zeros(data_length, Beam_Num, 'single');
fprintf('Sum_out_2D size: %d x %d\n', size(Sum_out_2D,1), size(Sum_out_2D,2));

%% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% -------------------------------
% Start beamforming
%% -------------------------------
for beam_idx = 1:Beam_Num
    
    %% 1. Select active elements (maintaining original logic)
    active_center = round((beam_idx-1)*(Elem_Num-1)/(Beam_Num-1)) + 1;
    active_start  = max(1, active_center - floor(Ch_Num/2));
    active_end    = min(Elem_Num, active_start + Ch_Num - 1);
    active_start  = max(1, active_end - Ch_Num + 1);
    active_idx    = active_start:active_end;
    
    %% 2. Load RF data
    fid = fopen(sprintf('RF_data/RxScanline%d.bin', beam_idx-1), 'rb');
    rf_raw = fread(fid, [data_length, Ch_Num], 'int16');
    fclose(fid);
    
    %% 3. Interpolate (upsample) RF data
    % The physical signal length remains unchanged,
    % but internal calculations use the interpolated data (fine resolution)
    fine_length = data_length * inter_factor;  % Number of samples for internal calculations
    rf_upsampled = zeros(fine_length, Ch_Num, 'single');
    for ch = 1:Ch_Num
        temp = resample(double(rf_raw(:,ch)), inter_factor, 1);
        rf_upsampled(:, ch) = single(temp) / 32768.0;  % Convert int16 to float and adjust scaling
    end
    
    %% 4. Apodization (window function)
    apod = kaiser(Ch_Num, 3)';  % Use Hann window (as in the original example)
    
    %% 5. Prepare dynamic focusing and steering angle calculation
    % Effective view angle considering both sides (e.g., 160 degrees)
    effective_view_angle = View_Angle * 2;
    effective_delta_angle = effective_view_angle/(Beam_Num-1);
    theta_deg = -effective_view_angle/2 + (beam_idx-1)*effective_delta_angle;
    theta_rad = deg2rad(theta_deg);
    
    active_x = x_element(active_idx);
    center_x = mean(active_x);
    
    %% 6. Dynamic focusing: output remains at the original data_length
    beam_sum = zeros(data_length,1,'single');
    
    for sample_idx = 1:data_length
        
        % Current depth (physical depth) [mm]: using the original depth interval (Unit_Dis)
        current_depth = (sample_idx-1) * Unit_Dis;
        
        % Calculate focal coordinates considering beam steering
        x_focus = center_x + current_depth * sin(theta_rad);
        z_focus = current_depth * cos(theta_rad);
        
        % Reference distance (from array center to the focal point)
        R_ref = sqrt( (center_x - x_focus)^2 + z_focus^2 );
        
        % The corresponding base index in the interpolated data for the current sample in the original signal
        base_index_fine = (sample_idx-1)*inter_factor + 1;
        
        sum_value = 0;
        % For each channel, apply delay compensation and sum the contributions
        for ch = 1:Ch_Num
            x_ch = active_x(ch);
            
            % Calculate the actual distance from the element to the focal point
            R_ch = sqrt( (x_ch - x_focus)^2 + z_focus^2 );
            delta_R = R_ch - R_ref;
            
            % Time delay [s]
            delay_time = delta_R / Sound_Speed;
            
            % Delay in the interpolated data (in fine sample units)
            delay_samples_fine = delay_time * Fs_new;
            
            % Fractional index in the interpolated data to reference
            sample_in_fine = base_index_fine - delay_samples_fine;
            
            % Linear interpolation: only process if within valid range
            if sample_in_fine >= 1 && sample_in_fine < fine_length
                idx_floor = floor(sample_in_fine);
                idx_ceil  = idx_floor + 1;
                if idx_ceil > fine_length
                    idx_ceil = fine_length;
                end
                frac = sample_in_fine - idx_floor;
                value = (1-frac) * rf_upsampled(idx_floor, ch) + frac * rf_upsampled(idx_ceil, ch);
                
                sum_value = sum_value + value * apod(ch);
            end
        end
        
        beam_sum(sample_idx) = sum_value;
    end
    
    Sum_out_2D(:, beam_idx) = beam_sum;
    
    fprintf('Beam %d/%d completed (active elements: %d-%d)\n', beam_idx, Beam_Num, active_start, active_end);
end

%% -------------------------------
% Save final data (data_length x Beam_Num)
%% -------------------------------
fid = fopen(fullfile(output_dir, 'sum_out.bin'), 'wb');
fwrite(fid, Sum_out_2D, 'double');
fclose(fid);
