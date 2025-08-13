clear; clc; close all;

%% Parameters
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
N = 50; % number of sensors (base value)
numSatellites = 50; % number of satellites (base value)
activationProbability = 1/500; % Base activation probability

noisePower = 1/10^(SNR/10); % Noise power based on SNR

% Monte Carlo iteration number
MC = 1e4;

% Transmit power
Pt = 10;
delta = 1;

% Satellite orbital parameters (rotational kinematics)
orbitPeriod = 90*60; % 90 minutes in seconds
omega = 2*pi/orbitPeriod; % angular velocity (rad/s)
timeThisMC = 575.5e-3; % time step for Monte Carlo iterations

% Pre-calculation of K-factor values for all possible angles
angleRange = 0:90;
K_values = getKfromAngles(angleRange);

%% Monte Carlo Loop 1: Varying Number of Satellites
fprintf('Running Monte Carlo Loop 1: Varying Number of Satellites\n');
satelliteRange = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
aoiVsSats = zeros(size(satelliteRange));

for satIdx = 1:length(satelliteRange)
    currentNumSats = satelliteRange(satIdx);
    fprintf('Processing %d satellites...\n', currentNumSats);
    
    % Initialize satellite positions for current number
    initialAngles = unifrnd(0, 2*pi, [1 currentNumSats]);
    
    % For AoI tracking
    aoiUsers = ones(1, N);
    totalAvgAoI = 0;
    
    % Monte Carlo simulation for current number of satellites
    for mc = 1:MC
        currentTime = mc * timeThisMC;
        
        % Update satellite positions and calculate elevation angles
        satelliteElevations = zeros(1, currentNumSats);
        for sat = 1:currentNumSats
            orbitAngle = initialAngles(sat) + omega * currentTime;
            orbitAngle = mod(orbitAngle, 2*pi);
            elevationAngle = rad2deg(orbitAngle);
            
            if elevationAngle > 90
                elevationAngle = 180 - elevationAngle;
            end
            
            if elevationAngle < 10
                elevationAngle = 0;
            end
            
            if elevationAngle > 10 && elevationAngle <= 90
                satelliteElevations(sat) = elevationAngle;
            end
        end
        
        % Remove satellites that are below horizon
        validSats = (satelliteElevations >= 10) & (satelliteElevations <= 90);
        satelliteElevations = satelliteElevations(validSats);
        numValidSats = sum(validSats);
        
        if numValidSats == 0
            aoiUsers = aoiUsers + 1;
            continue;
        end
        
        % Pre-calculate PDFs and CDFs for valid satellites
        satellitePDFs = cell(1, numValidSats);
        satelliteCDFs = cell(1, numValidSats);
        P = 0.001:0.001:10;
        
        for sat = 1:numValidSats
            ElevDeg = satelliteElevations(sat);
            K_dB = getKfromAngles(ElevDeg);
            K_linear = 10^(K_dB/10);
            
            % Calculate Rician PDF
            pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
            areaUnderPDF = trapz(P, pdf_RX);
            pdfNormalized = pdf_RX / areaUnderPDF;
            cdfRx = cumtrapz(P, pdfNormalized);
            
            satellitePDFs{sat} = pdfNormalized;
            satelliteCDFs{sat} = cdfRx;
        end
        
        % Find active users
        FindActiveUsers = (aoiUsers >= delta) & (rand(1,N) < activationProbability);
        activeUsers = find(FindActiveUsers);
        receivePowerUsers = zeros(numValidSats, N);
        captureSensor = zeros(1, numValidSats);

        if ~isempty(activeUsers)
            % Generate received power for each active user at each satellite
            for j = 1:length(activeUsers)
                currentActiveUser = activeUsers(j);
                for sat = 1:numValidSats
                    u = rand();
                    cdfRx = satelliteCDFs{sat};
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
            for sat = 1:numValidSats
                userForThisSat = activeUsers;
                if ~isempty(userForThisSat)
                    activePower = receivePowerUsers(sat, userForThisSat);
                    [maxPower, positionOfMaxPower] = max(activePower);
                    strongestUser = userForThisSat(positionOfMaxPower);

                    % Check for deep fading
                    if maxPower < capture_th
                        continue;
                    end

                    % SINR calculation
                    noisePowerIt = randn() * sqrt(noisePower);
                    interference = sum(activePower) - maxPower;
                    SINR = maxPower / (interference + noisePowerIt);

                    % Check if SINR meets threshold
                    if 10*log10(SINR) >= 3
                        captureSensor(sat) = strongestUser;
                    end
                end
            end
        end

        % AoI update
        aoiUsers = aoiUsers + 1;

        % Reset AoI for captured sensors
        successfulCaptures = find(captureSensor > 0);
        for idx = 1:length(successfulCaptures)
            sat = successfulCaptures(idx);
            capturedUser = captureSensor(sat);
            aoiUsers(capturedUser) = 0;
        end
        
        totalAvgAoI = totalAvgAoI + mean(aoiUsers);
    end
    
    aoiVsSats(satIdx) = totalAvgAoI / MC;
end

%% Monte Carlo Loop 2: Varying Activation Probability
fprintf('Running Monte Carlo Loop 2: Varying Activation Probability\n');
probRange = [1/2, 1/3, 1/5, 1/10, 1/20, 1/50, 1/100, 1/500];
aoiVsProb = zeros(size(probRange));

% Initialize satellite positions (using base number)
initialAngles = unifrnd(0, 2*pi, [1 numSatellites]);

for probIdx = 1:length(probRange)
    currentProb = probRange(probIdx);
    fprintf('Processing activation probability 1/%.0f...\n', 1/currentProb);
    
    % For AoI tracking
    aoiUsers = ones(1, N);
    totalAvgAoI = 0;
    
    % Monte Carlo simulation for current activation probability
    for mc = 1:MC
        currentTime = mc * timeThisMC;
        
        % Update satellite positions and calculate elevation angles
        satelliteElevations = zeros(1, numSatellites);
        for sat = 1:numSatellites
            orbitAngle = initialAngles(sat) + omega * currentTime;
            orbitAngle = mod(orbitAngle, 2*pi);
            elevationAngle = rad2deg(orbitAngle);
            
            if elevationAngle > 90
                elevationAngle = 180 - elevationAngle;
            end
            
            if elevationAngle < 10
                elevationAngle = 0;
            end
            
            if elevationAngle > 10 && elevationAngle <= 90
                satelliteElevations(sat) = elevationAngle;
            end
        end
        
        % Remove satellites that are below horizon
        validSats = (satelliteElevations >= 10) & (satelliteElevations <= 90);
        satelliteElevations = satelliteElevations(validSats);
        numValidSats = sum(validSats);
        
        if numValidSats == 0
            aoiUsers = aoiUsers + 1;
            continue;
        end
        
        % Pre-calculate PDFs and CDFs for valid satellites
        satellitePDFs = cell(1, numValidSats);
        satelliteCDFs = cell(1, numValidSats);
        P = 0.001:0.001:10;
        
        for sat = 1:numValidSats
            ElevDeg = satelliteElevations(sat);
            K_dB = getKfromAngles(ElevDeg);
            K_linear = 10^(K_dB/10);
            
            % Calculate Rician PDF
            pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
            areaUnderPDF = trapz(P, pdf_RX);
            pdfNormalized = pdf_RX / areaUnderPDF;
            cdfRx = cumtrapz(P, pdfNormalized);
            
            satellitePDFs{sat} = pdfNormalized;
            satelliteCDFs{sat} = cdfRx;
        end
        
        % Find active users
        FindActiveUsers = (aoiUsers >= delta) & (rand(1,N) < currentProb);
        activeUsers = find(FindActiveUsers);
        receivePowerUsers = zeros(numValidSats, N);
        captureSensor = zeros(1, numValidSats);

        if ~isempty(activeUsers)
            % Generate received power for each active user at each satellite
            for j = 1:length(activeUsers)
                currentActiveUser = activeUsers(j);
                for sat = 1:numValidSats
                    u = rand();
                    cdfRx = satelliteCDFs{sat};
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
            for sat = 1:numValidSats
                userForThisSat = activeUsers;
                if ~isempty(userForThisSat)
                    activePower = receivePowerUsers(sat, userForThisSat);
                    [maxPower, positionOfMaxPower] = max(activePower);
                    strongestUser = userForThisSat(positionOfMaxPower);

                    % Check for deep fading
                    if maxPower < capture_th
                        continue;
                    end

                    % SINR calculation
                    noisePowerIt = randn() * sqrt(noisePower);
                    interference = sum(activePower) - maxPower;
                    SINR = maxPower / (interference + noisePowerIt);

                    % Check if SINR meets threshold
                    if 10*log10(SINR) >= 3
                        captureSensor(sat) = strongestUser;
                    end
                end
            end
        end

        % AoI update
        aoiUsers = aoiUsers + 1;

        % Reset AoI for captured sensors
        successfulCaptures = find(captureSensor > 0);
        for idx = 1:length(successfulCaptures)
            sat = successfulCaptures(idx);
            capturedUser = captureSensor(sat);
            aoiUsers(capturedUser) = 0;
        end
        
        totalAvgAoI = totalAvgAoI + mean(aoiUsers);
    end
    
    aoiVsProb(probIdx) = totalAvgAoI / MC;
end

%% Monte Carlo Loop 3: Varying Number of Devices
fprintf('Running Monte Carlo Loop 3: Varying Number of Devices\n');
deviceRange = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
aoiVsDevices = zeros(size(deviceRange));

% Initialize satellite positions (using base number)
initialAngles = unifrnd(0, 2*pi, [1 numSatellites]);

for devIdx = 1:length(deviceRange)
    currentN = deviceRange(devIdx);
    fprintf('Processing %d devices...\n', currentN);
    
    % For AoI tracking
    aoiUsers = ones(1, currentN);
    totalAvgAoI = 0;
    
    % Monte Carlo simulation for current number of devices
    for mc = 1:MC
        currentTime = mc * timeThisMC;
        
        % Update satellite positions and calculate elevation angles
        satelliteElevations = zeros(1, numSatellites);
        for sat = 1:numSatellites
            orbitAngle = initialAngles(sat) + omega * currentTime;
            orbitAngle = mod(orbitAngle, 2*pi);
            elevationAngle = rad2deg(orbitAngle);
            
            if elevationAngle > 90
                elevationAngle = 180 - elevationAngle;
            end
            
            if elevationAngle < 10
                elevationAngle = 0;
            end
            
            if elevationAngle > 10 && elevationAngle <= 90
                satelliteElevations(sat) = elevationAngle;
            end
        end
        
        % Remove satellites that are below horizon
        validSats = (satelliteElevations >= 10) & (satelliteElevations <= 90);
        satelliteElevations = satelliteElevations(validSats);
        numValidSats = sum(validSats);
        
        if numValidSats == 0
            aoiUsers = aoiUsers + 1;
            continue;
        end
        
        % Pre-calculate PDFs and CDFs for valid satellites
        satellitePDFs = cell(1, numValidSats);
        satelliteCDFs = cell(1, numValidSats);
        P = 0.001:0.001:10;
        
        for sat = 1:numValidSats
            ElevDeg = satelliteElevations(sat);
            K_dB = getKfromAngles(ElevDeg);
            K_linear = 10^(K_dB/10);
            
            % Calculate Rician PDF
            pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
            areaUnderPDF = trapz(P, pdf_RX);
            pdfNormalized = pdf_RX / areaUnderPDF;
            cdfRx = cumtrapz(P, pdfNormalized);
            
            satellitePDFs{sat} = pdfNormalized;
            satelliteCDFs{sat} = cdfRx;
        end
        
        % Find active users
        FindActiveUsers = (aoiUsers >= delta) & (rand(1,currentN) < activationProbability);
        activeUsers = find(FindActiveUsers);
        receivePowerUsers = zeros(numValidSats, currentN);
        captureSensor = zeros(1, numValidSats);

        if ~isempty(activeUsers)
            % Generate received power for each active user at each satellite
            for j = 1:length(activeUsers)
                currentActiveUser = activeUsers(j);
                for sat = 1:numValidSats
                    u = rand();
                    cdfRx = satelliteCDFs{sat};
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
            for sat = 1:numValidSats
                userForThisSat = activeUsers;
                if ~isempty(userForThisSat)
                    activePower = receivePowerUsers(sat, userForThisSat);
                    [maxPower, positionOfMaxPower] = max(activePower);
                    strongestUser = userForThisSat(positionOfMaxPower);

                    % Check for deep fading
                    if maxPower < capture_th
                        continue;
                    end

                    % SINR calculation
                    noisePowerIt = randn() * sqrt(noisePower);
                    interference = sum(activePower) - maxPower;
                    SINR = maxPower / (interference + noisePowerIt);

                    % Check if SINR meets threshold
                    if 10*log10(SINR) >= 3
                        captureSensor(sat) = strongestUser;
                    end
                end
            end
        end

        % AoI update
        aoiUsers = aoiUsers + 1;

        % Reset AoI for captured sensors
        successfulCaptures = find(captureSensor > 0);
        for idx = 1:length(successfulCaptures)
            sat = successfulCaptures(idx);
            capturedUser = captureSensor(sat);
            aoiUsers(capturedUser) = 0;
        end
        
        totalAvgAoI = totalAvgAoI + mean(aoiUsers);
    end
    
    aoiVsDevices(devIdx) = totalAvgAoI / MC;
end

%% Results Plotting
% 1. Plot AoI vs Number of Satellites
figure;
plot(satelliteRange, aoiVsSats, '-o', 'LineWidth', 2);
xlabel('Number of Satellites');
ylabel('Average AoI');
title('Average AoI vs Number of Satellites');
grid on;

% 2. Plot AoI vs Activation Probability
figure;
probLabels = {'1/2', '1/3', '1/5', '1/10', '1/20', '1/50', '1/100', '1/500'};
plot(1:length(probRange), aoiVsProb, '-o', 'LineWidth', 2);
xlabel('Activation Probability');
ylabel('Average AoI');
title('Average AoI vs Activation Probability');
set(gca, 'XTick', 1:length(probRange), 'XTickLabel', probLabels);
grid on;

% 3. Plot AoI vs Number of Devices
figure;
plot(deviceRange, aoiVsDevices, '-o', 'LineWidth', 2);
xlabel('Number of Devices');
ylabel('Average AoI');
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