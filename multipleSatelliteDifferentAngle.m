clear; clc; close all;

% Parameters
f_T = 2400e6; % frequency 2.4 GHz
Re = 6378e3; % earth radius
h_s = 400e3; % satellite altitude
h_u = 0; % ground sensor position
G_t = 1; % sensor antenna gain
G_s = 10^(15/10); % satellite antenna gain
Lta = 10^(-1/10); % transmitter antenna gain
Lra = 10^(-0.1/10); % receiver antenna gain
La = 10^(-0.2/10); % atmospheric losses
c = 3e8; % speed of light
sigma_shadow = 3; % shadowing standard deviation
snr_th = 1; % SNR threshold
capture_th = 1; % capture threshold
aoiThreshold = 1; % AoI threshold
k_boltz = 1.380649e-23; % Boltzmann constant
T_sys = 1000; % system temperature
B = 128e3; % bandwidth
SNR = 10; % dB
N = 50; % number of sensors
numSatellites = 50; % number of satellites
activationProbabilityBase = 1/N; % Base activation probability

noisePower = 1/10^(SNR/10); % Noise power based on SNR

% Monte Carlo iteration number
MC = 1e4;

% Transmit power
Pt = 10;
delta = 1;

% Different elevation angles for each satellite
minElevation = 10;
maxElevation = 90;
satelliteElevations = linspace(minElevation, maxElevation, numSatellites);

% Pre-calculation of satellites K-factor, distance, PDF, and CDF
K_values = zeros(1, numSatellites);
d = zeros(1, numSatellites);
satellitePDFs = cell(1, numSatellites);
satelliteCDFs = cell(1, numSatellites);

P = 0.001:0.001:10;

for sat = 1:numSatellites
    ElevDeg = satelliteElevations(sat);
    ElevRad = deg2rad(ElevDeg);

    K_values(sat) = getKfromAngles(ElevDeg);
    K_linear = 10.^(K_values(sat)/10);

    % PDF & CDF
    pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
    cdfRx = cumtrapz(P, pdf_RX);

    areaUnderPDF = trapz(P, pdf_RX);
    pdfNormalized = pdf_RX / areaUnderPDF;

    satellitePDFs{sat} = pdfNormalized;
    satelliteCDFs{sat} = cdfRx;
end

% For AoI tracking
aoiUsers = ones(1, N);
avgAoI = zeros(1, MC);

% Storage for successful captures
captureData = [];
captureCounter = 0;

% Monte Carlo simulation
for mc = 1:MC
    % Vary activation probability
    activationProbability = activationProbabilityBase * rand();

    FindActiveUsers = (aoiUsers >= delta) & (rand(1,N) < activationProbability);
    activeUsers = find(FindActiveUsers);
    receivePowerUsers = zeros(numSatellites, N);
    captureSensor = zeros(1, numSatellites); % to store captured user per satellite

    if ~isempty(activeUsers)
        % Generate received power for each active user at each satellite
        for j = 1:length(activeUsers)
            currentActiveUser = activeUsers(j);
            for sat = 1:numSatellites
                u = rand();
                cdfRx = satelliteCDFs{sat}; % Precomputed CDF for satellite
                possibleValues = find(cdfRx < u);
                if ~isempty(possibleValues)
                    samplePower = P(possibleValues(end));
                else
                    samplePower = P(1);
                end
                receivePowerUsers(sat, currentActiveUser) = samplePower;
            end
        end

        % Capture effect per satellite
        for sat = 1:numSatellites
            userForThisSat = activeUsers; % Consider all active users
            if ~isempty(userForThisSat)
                activePower = receivePowerUsers(sat, userForThisSat);
                [maxPower, positionOfMaxPower] = max(activePower);
                strongestUser = userForThisSat(positionOfMaxPower);

                % Check for deep fading
                if maxPower < capture_th
                    continue; % Skip to next satellite
                end

                % SINR calculation
                noisePowerIt = randn() * sqrt(noisePower);
                interference = sum(activePower) - maxPower;
                SINR = maxPower / (interference + noisePowerIt);

                % Check if SINR meets threshold (e.g., 3 dB)
                if 10*log10(SINR) >= 3
                    captureSensor(sat) = strongestUser;
                end
            end
        end
    end

    % AoI update
    aoiUsers = aoiUsers + 1;

    % Reset AoI for captured sensors and record capture data
    successfulCaptures = find(captureSensor > 0);
    for satelliteCapturingPacket = 1:length(successfulCaptures)
        sat = successfulCaptures(satelliteCapturingPacket);
        capturedUser = captureSensor(sat);
        currentAoI = aoiUsers(capturedUser);
        aoiUsers(capturedUser) = 0; % Reset AoI for captured user (only once per user)

        % Record capture data
        captureCounter = captureCounter + 1;
        captureData(captureCounter, :) = [satelliteElevations(sat), currentAoI, mc, sat, capturedUser];

        % Print capture information
        if mod(mc, 1000) == 0
            fprintf('MC: %d | Satellite %d (Elev: %.1f°) captured from Sensor %d | AoI: %.3f\n', ...
                mc, sat, satelliteElevations(sat), capturedUser, currentAoI);
        end
    end

    % Store average AoI for this iteration
    avgAoI(mc) = mean(aoiUsers);
end

% Display results
fprintf('Total successful captures: %d\n', captureCounter);
fprintf('Average AoI across all iterations: %.3f\n', mean(avgAoI));

% Satellite configuration
fprintf('\nSatellite Configuration:\n');
for sat = 1:numSatellites
    fprintf('Satellite %d: Elevation = %.1f°, K-factor = %.1f dB\n', sat, satelliteElevations(sat), K_values(sat));
end

% Plot AoI vs Elevation angle
figure;
avgAoIPerSatellite = zeros(1, numSatellites);
captureCountPerSatellite = zeros(1, numSatellites);

for sat = 1:numSatellites
    satCaptures = captureData(captureData(:,4) == sat, :);
    if ~isempty(satCaptures)
        avgAoIPerSatellite(sat) = mean(satCaptures(:,2));
        captureCountPerSatellite(sat) = size(satCaptures, 1);
    end
end

plot(satelliteElevations, avgAoIPerSatellite);
xlabel('Satellite Elevation Angle (degrees)');
ylabel('Average AoI at Capture');
title('Average AoI vs Satellite Elevation Angle');
grid on;

% Plot AoI vs Number of Satellites
figure;
satelliteRange = 1:1:numSatellites;
avgAoIPerSatelliteCount = zeros(1, length(satelliteRange));

for i = 1:length(satelliteRange)
    satCount = satelliteRange(i);
    satCaptures = captureData(captureData(:,4) <= satCount, :);
    if ~isempty(satCaptures)
        avgAoIPerSatelliteCount(i) = mean(satCaptures(:,2));
    end
end

plot(satelliteRange, avgAoIPerSatelliteCount);
xlabel('Number of Satellites');
ylabel('Average AoI at Capture');
title('Average AoI vs Number of Satellites');
grid on;

% Plot AoI vs Number of Devices
figure;
deviceRange = 1:1:N;
avgAoIPerDeviceCount = zeros(1, length(deviceRange));

for i = 1:length(deviceRange)
    deviceCount = deviceRange(i);
    deviceCaptures = captureData(captureData(:,5) <= deviceCount, :);
    if ~isempty(deviceCaptures)
        avgAoIPerDeviceCount(i) = mean(deviceCaptures(:,2));
    end
end

plot(deviceRange, avgAoIPerDeviceCount);
xlabel('Number of Devices');
ylabel('Average AoI at Capture');
title('Average AoI vs Number of Devices');
grid on;

% Function to get K-factor from elevation angle with validation
function K = getKfromAngles(angleElevation)
    K_values = [0.9, 1.5, 2.2, 4.1, 8.9, 11.4, 13.5, 15.2, 18.6];
    K_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];

    if angleElevation < 10
        K = K_values(1);
    elseif angleElevation > 90
        K = K_values(end);
    else
        K = interp1(K_angles, K_values, angleElevation, 'linear');
    end
end