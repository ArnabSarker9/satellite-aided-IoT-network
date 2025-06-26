clear; clc; close all;

% parameters
f_T = 2.4e9; % Frequency band 
Pt_values = [2, 5, 10, 15];
G_I = 0;          % Sensor antenna gain [dBi] (from Table II)
G_s = 15;            % Satellite antenna gain [dBi] (from Table II)
G_g = 30;             % Ground station antenna gain [dBi] (assumed)
T_sys = 1000;         % Receiver system noise temperature in K
Re = 6371e3;          % Earth radius in meters (6,371 km)
h_s = 600e3;          % Altitude of satellite (600 km)
h_u = 0;              % Altitude of IoT sensors (fixed at ground level)
c = 3e8;              % Speed of light in m/s (300,000,000 m/s)
sigma_shadow = 3;     % Shadowing standard deviation (dB) (reduced for debugging)
capture_th = 1;       % Capture threshold in dB

N = 100; % number of sensors

% Generate random sensor positions (latitude/longitude)
sensor_lat = 180 * rand(1, N) - 90;   % Latitude between -90° and 90°
sensor_lon = 360 * rand(1, N) - 180;  % Longitude between -180° and 180°

MC = 1000; % number of monte carlo simulation
delta = 1; % AoI threshold
MC = 1000;              % Number of Monte Carlo simulations
delta = 1;              % AoI threshold
AoI_users = ones(1, N); % AoI for each user is 1 (fresh info)
AverageAoI = zeros(1, MC); % Store average AoI for each Monte Carlo iteration
peak_AoI = 0;           % peak AoI
peak_AoI_num = 0;       % Count of peak AoI occurrences

% Initialize AoI storage for each Pt and elevation angle
AoI_phase2 = zeros(N, length(Pt_values));  % Store AoI values for different Pt values
elev_angles = linspace(0, 90, N);  % Array for elevation angles (0° to 90°)

% monte carlo simulation
for mc = 1:MC

    % Phase 1: IoT sensors decide to transmit based on AoI threshold
    transmitting_users = (rand(1, N) < 0.5) .* (AoI_users >= delta);
    transmitting_indices = find(transmitting_users);

    % Initialize received powers
    rx_powers_dBm = -inf(1, N); % Initialize with very low power
    for j = 1: length(Pt_values)
        Pt = Pt_values(j);
        P_T = 10 * log10(Pt);  % transmit power in dBW

        % Phase 2: Calculate received power for transmitting users
        for user = transmitting_indices
        % Generate random elevation angle 
        theta = rand * (pi/2);
        
        % Calculate distance
        d = calculate_distance(Re, h_u, h_s, theta);
        
        % Calculate K-factor using the provided formula
        K = calculate_K(f_T, Re, h_u, h_s, theta, c);
        
        % log-normal shadowing
        shadowing_dB = sigma_shadow * randn;
        
        % Calculate received power in dBm
        rx_powers_dBm(user) = calculate_received_power(P_T, G_I, G_s, G_g, K, shadowing_dB);

        AoI_phase2(user, j) = AoI_users(user);
        end
    end
        % Phase 3: Capture effect - find strongest signal
    captured_user = 0;
    if ~isempty(transmitting_indices)
        % Find the strongest signal
        [max_power_dBm, max_idx] = max(rx_powers_dBm(transmitting_indices));
        max_user = transmitting_indices(max_idx);
        
        % Calculate interference from other transmitters
        other_powers = rx_powers_dBm(transmitting_indices);
        other_powers(max_idx) = []; % Remove the strongest signal
        
        if isempty(other_powers)
            % Only one transmitter - automatic capture
            captured_user = max_user;
        else
            % Convert to linear for SIR calculation
            max_power_linear = 10.^(max_power_dBm/10);
            interference_linear = sum(10.^(other_powers/10));
            
            % Calculate Signal-to-Interference Ratio (SIR)
            SIR_dB = 10*log10(max_power_linear / interference_linear);
            disp(['SIR for captured user ', num2str(max_user), ': ', num2str(SIR_dB), ' dB']);
            
            if SIR_dB >= capture_th
                captured_user = max_user;
            end
        end
    end

    % Phase 4: Update AoI based on captured power
    for user = 1:N
        AoI_users(user) = AoI_users(user) + 1;  % Increment AoI for all users
        
        % Reset AoI for captured user
        if user == captured_user
            peak_AoI = peak_AoI + AoI_users(user);
            peak_AoI_num = peak_AoI_num + 1;
            AoI_users(user) = 0;  % Reset AoI for successful transmission
        end
    end

        % Phase 5: Calculate average AoI for this iteration
    AverageAoI(mc) = mean(AoI_users);
    
    % Display progress
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


function d = calculate_distance(Re, h_u, h_s, theta)
    % Calculate distance using geometric formula
    term1 = (Re + h_s)^2 - (Re + h_u)^2 * cos(theta)^2;
    term2 = (Re + h_u) * sin(theta);
    d = sqrt(term1) - term2;
end

function K = calculate_K(f_T, Re, h_u, h_s, theta, c)
    % Calculate K-factor using the provided formula
    numerator = c;
    denominator = 4 * pi * f_T * ...
                 ( - (Re + h_u) * sin(theta) + ...
                   sqrt( (Re + h_s)^2 - (Re + h_u)^2 * cos(theta)^2 ) );
    K = (numerator / denominator)^2;
end

function P_r_dBm = calculate_received_power(P_T, G_I, G_s, G_g, K, shadowing_dB)
    % Convert transmit power to dBm
    P_T_dBm = 10*log10(P_T*1000);
    
    % Calculate received power (dBm)
    P_r_dBm = P_T_dBm + G_I + G_s + G_g + 10*log10(K) + shadowing_dB;
end

% Plot AoI vs Elevation Angle for different Pt values
figure;
hold on;

% Loop through each Pt value and plot AoI vs Elevation Angle
for j = 1:length(Pt_values)
    plot(elev_angles, AoI_phase2(:, j), '-o', 'LineWidth', 2);
end

% Customize the plot
xlabel('Elevation Angle (degrees)');
ylabel('Age of Information (seconds)');
title('AoI vs Elevation Angle');
grid on;
legend(arrayfun(@(pt) sprintf('Pt = %d', pt), Pt_values, 'UniformOutput', false));
hold off;
