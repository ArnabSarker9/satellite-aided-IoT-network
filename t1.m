clear; clc; close all;

% Parameters
f_T = 2400e6; % frequency 2.4 GHz
Re = 6378e3; % earth radius
h_s = 400e3; % satellite altitude
h_u = 0; % ground sensor position
G_t = 0; % sensor antenna gain
G_s = 10^(15/10); % satellite antenna gain
Lta = 10^(-1/10); % transmitter antenna losses
Lra = 10^(-0.1/10); % receiver antenna losses
La = 10^(-0.2/10); % atmospheric losses
c = 3e8; % speed of light
sigma_shadow = 3; % shadowing standard deviation
snr_th = 1; % SNR threshold
capture_th = 1; % capture threshold
aoiThreshold = 1; % AoI threshold
k_boltz = 1.380649e-23;% Boltzmann constant
T_sys = 1000;  % System temperature (K)
B = 128e3; % Bandwidth (Hz)

% number of sensors
N = 10;
activationProbability = 1/N;

% noise
noise_power = k_boltz * T_sys * B;

% monte carlo
MC = 1e5;

% transmit power
Pt = 10;

% elevation angel
ElevAngle = 10:10:90;
avgAoIresult = zeros(size(ElevAngle));

for angleElevation = 1:length(ElevAngle)
    ElevDeg = ElevAngle(angleElevation);
    ElevRad = deg2rad(ElevDeg);

    aoiUsers = ones(1, N);
    K = getKfromAngles(angleElevation);
    d = sqrt((Re + h_s)^2 - (Re + h_u)^2 * cos(angleElevation).^2) - ((Re + h_u) * sin (angleElevation)) + rand(1,N)*10; % distance
    pathLoss = (c ./ (4 * pi * f_T * d)).^2;
    rxPower = Pt * G_s * G_t * pathLoss;
    % receiverPDF = K * exp(-K .* (rxPower +1)) .* besseli (0, 2 * K .* sqrt(rxPower));

    aoiCalculation = ones(1, N);
    avgAoI = zeros(1, MC);

    % monte carlo simulation
    for mc = 1: MC
        FindActiveUsers = (aoiUsers >= delta) & (randn(1,N) <activationProbability);
        activeUsers = find(FindActiveUsers);
        receiverPowerUsers = zeros(1,N);
        captureSensor = 0;

        if ~isempty(activeUsers)
            K_linear = 10^(K/10);

            for j = 1:length(activeUsers)
                users = activeUsers (j);

                s_dB = rand(0, 1, [1,N]);
                shadowing = 10^(s_dB/10);

                receiverPowerUsers(users) = Pt * G_t * G_s * shadowing * Lta * Lra * La;
            end

            % capture effect
            activePower = receiverPowerUsers(activeUsers);
            [maxPower, highestReceivedPower] = max(activePower);

            % check SNR
            if maxPower < noise_power * snr_th
                captureSensor = 0;
            else
                % check SINR
                interference = sum(activePower) - maxPower;
                SINR = maxPower / (interference + noise_power);

                if SINR >= capture_th
                    captureSensor = activeUsers(highestReceivedPower);
                end
            end
        end
        %  AoI
        aoiCalculation = aoiCalculation + 1;

        if captureSensor > 0
            aoiCalculation(captureSensor) = 0;
        end
        % storing results
        avgAoI(mc) = mean(aoiCalculation);
    end
    % stored AoI
    avgAoIresult(angleElevation) = mean(avgAoI);
    fprintf('Elevation: %d° completed | Avg AoI: %.2f\n', ElevDeg, avgAoIresult);
end

figure;
plot(ElevAngle, avgAoIresult, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Elevation Angle (degrees)');
ylabel('Average Age of Information (slots)');
title('AoI Performance vs. Satellite Elevation Angle');
xticks(ElevAngle);

function K = getKfromAngles(angleElevation)
K_values = [0.9, 1.5, 2.2, 4.1, 8.9, 11.4, 13.5, 15.2, 18.6];
K_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
K = K_values(K_angles==angleElevation);
end
