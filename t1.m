clear; clc; close all;

% Parameters
f_T = 2.4e9; % Frequency band
Pt_values = [2, 5, 10, 15]; % Transmit power values
G_I = 0;          % Sensor antenna gain
G_s = 15;         % Satellite antenna gain
G_g = 30;         % Ground station antenna gain (assumed)
T_sys = 1000;     % Receiver system noise temperature in K
Re = 6371e3;      % Earth radius in meters
h_s = 600e3;      % Altitude of satellite
h_u = 0;          % Altitude of IoT sensors
c = 3e8;          % Speed of light in m/s
sigma_shadow = 3; % Shadowing standard deviation (dB)
capture_th = 1;   % Capture threshold in dB

N = 100; % Number of sensors

% sensor positions 
sensor_lat = 180 * rand(1, N) - 90;   % Latitude between -90° and 90°
sensor_lon = 360 * rand(1, N) - 180;  % Longitude between -180° and 180°

MC = 1000; % Number of Monte Carlo simulations
delta = 1;  % AoI threshold
AoI_users = ones(1, N); % AoI for each user is 1 (fresh info)
AverageAoI = zeros(1, MC); % average AoI 
peak_AoI = 0;             % Peak AoI
peak_AoI_num = 0;         % Count of peak AoI occurrences

% Initialize AoI storage and elevation angle
AoI_phase2 = zeros(N, length(Pt_values));  
elev_angles = linspace(0, 90, N);

% Monte Carlo simulation
for mc = 1:MC
    % IoT sensors decide to transmit based on AoI threshold
    transmitting_users = (rand(1, N) < 0.5) .* (AoI_users >= delta);
    transmitting_indices = find(transmitting_users);

    % 1. Initialize received powers
    rx_powers_dBm = -inf(1, N); 
    for j = 1:length(Pt_values)
        Pt = Pt_values(j);
        P_T = Pt;  

        % 2. Calculate received power 
        for user = transmitting_indices
            theta = rand * (pi/2);  % Elevation angle 

            % Calculate distance
            d = calculate_distance(Re, h_u, h_s, theta);

            % Calculate K-factor
            P_r = calculate_received_power(P_T, G_I, G_s, G_g, 0, sigma_shadow);  
            K = calculate_K(f_T, Re, h_u, h_s, theta, c, P_r);

            % Log-normal shadowing
            shadowing_dB = sigma_shadow * randn;  % shadowing noise

            % Calculate received power in dBm
            rx_powers_dBm(user) = calculate_received_power(P_T, G_I, G_s, G_g, K, shadowing_dB);

            AoI_phase2(user, j) = AoI_users(user);
        end
    end

    % 3 Capture effect
    captured_user = 0;
        if ~isempty(transmitting_indices)
        % To find the strongest signal
            [max_power_dBm, max_idx] = max(rx_powers_dBm(transmitting_indices));
            max_user = transmitting_indices(max_idx);
    
            % Calculate interference from other transmitters
            other_powers = rx_powers_dBm(transmitting_indices);
            other_powers(max_idx) = [];  
    
        if isempty(other_powers)
            captured_user = max_user;
        else
            % Convert to linear for SIR calculation
            max_power_linear = 10.^(max_power_dBm / 10);
            interference_linear = sum(10.^(other_powers / 10));

        if interference_linear == 0
            SIR_dB = inf;  
        else
            % Signal-to-Interference Ratio
            SIR_dB = 10 * log10(max_power_linear / interference_linear);
        end
        
        % Check if SIR is above the capture threshold
        if SIR_dB >= capture_th
            captured_user = max_user;
            end
        end
    end

    % 4. Update AoI based on captured power
    for user = 1:N
        AoI_users(user) = AoI_users(user) + 1;  

        % Reset AoI for captured user
        if user == captured_user
            peak_AoI = peak_AoI + AoI_users(user);
            peak_AoI_num = peak_AoI_num + 1;
            AoI_users(user) = 0;
        end
    end

    % 5. Calculate average AoI
    AverageAoI(mc) = mean(AoI_users);

    if mod(mc, 100) == 0
        fprintf('Monte Carlo Iteration: %d\n', mc);
        fprintf('  Current Average AoI: %.2f\n', AverageAoI(mc));
        if captured_user > 0
            fprintf('  Captured sensor: %d at (%.2f°, %.2f°)\n', ...
                    captured_user, sensor_lat(captured_user), sensor_lon(captured_user));
        else
            fprintf('  No capture in this iteration\n');
        end
    end
end

% 6. Calculate distance
function d = calculate_distance(Re, h_u, h_s, theta)
    term1 = (Re + h_s)^2 - (Re + h_u)^2 * cos(theta)^2;
    term2 = (Re + h_u) * sin(theta);
    d = sqrt(term1) - term2;
end

% 7. Calculate K-factor
function K_final = calculate_K(f_T, Re, h_u, h_s, theta, c, P_r)
    numerator = c;
    denominator = 4 * pi * f_T * ...
                 ( - (Re + h_u) * sin(theta) + ...
                   sqrt( (Re + h_s)^2 - (Re + h_u)^2 * cos(theta)^2 ) );
               
    K_geom = (numerator / denominator)^2; 

    K_final = exp(-K_geom * (P_r + 1)) * besseli(0, 2 * K_geom * sqrt(P_r)); 
end

% 8. Calculate received power
function P_r_dBm = calculate_received_power(P_T, G_I, G_s, G_g, K, shadowing_dB)
    P_T_dBm = 10 * log10(P_T * 1000);  % Convert to dBm

    P_r_dBm = P_T_dBm + G_I + G_s + G_g + 10 * log10(K) + shadowing_dB;  % Including K and shadowing
end

% Plot AoI vs Elevation Angle
figure;
hold on;
plot(elev_angles, AoI_phase2(:, 1), '-o', 'LineWidth', 2);
plot(elev_angles, AoI_phase2(:, 2), '-s', 'LineWidth', 2);
plot(elev_angles, AoI_phase2(:, 3), '-^', 'LineWidth', 2);
plot(elev_angles, AoI_phase2(:, 4), '-d', 'LineWidth', 2);

xlabel('Elevation Angle (degrees)');
ylabel('Age of Information (seconds)');
title('AoI vs Elevation Angle');
grid on;
legend(arrayfun(@(pt) sprintf('Pt = %d', pt), Pt_values, 'UniformOutput', false));
hold off;