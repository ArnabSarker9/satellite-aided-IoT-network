clear; clc; close all;

% parameters
f_T = 2400e6; % frequency 2.4 GHz
Re = 6378e3; % earth radius
h_s = 400e3; % satellite altitude
h_u = 0; % ground sensor position
G_t = 1; % sensor antenna gain
G_s = 10^(15/10); % satellite antenna gain
Lta = 10^(-1/10); % transmitter antenna gain
Lra = 10^(-0.1/10); % receiver antenna gain
La = 10^(-0.2/10); % atmosåheric losses
c = 3e8; % speed of light
sigma_shadow = 3; % shadowing standard deviation
snr_th = 1; % SNR threshold
capture_th = 1; % capture threshold
aoiThreshold = 1; % aoi threshold
k_boltz = 1.380649e-23; % Boltzman constant
T_sys = 1000; % system temperature
B = 128e3; % bandwidth

N = 10; % number of sensors
numSatellites = 5; % number of satellites
activationProbability = 1/N;

noisePower = k_boltz * T_sys * B; % Noise

% Monter carlo interation number
MC = 1e5;

% Transmit power
Pt = 10;
delta = 1;

% elevation anagle
ElevAngle = 10:10:90;
finalAoIresult = zeros(size(ElevAngle));
P = 0.001:0.001:10;

for angleElevation = 1:length(ElevAngle)
    ElevDeg = ElevAngle(angleElevation);
    ElevRad = deg2rad(ElevDeg);

    % Calculate distance to satellite
    d = sqrt((Re + h_s)^2 - (Re + h_u)^2 * cos(ElevRad)^2) - ((Re + h_u) * sin(ElevRad));

    % satellites position, same elevation angle with the sensors, but the azimuth angle is different
    satellitePosition = linspace(0, 2 * pi, numSatellites+1);
    satellitePosition = satellitePosition(1:end-1);

    aoiUsers = ones(1, N);
    K = getKfromAngles(ElevDeg);  
    K_linear = 10.^(K/10); 
    
    % Calculate Rician PDF
    pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
    cdfRx = cumtrapz (P, pdf_RX);

    areaUnderPDF = trapz(P, pdf_RX);
    pdfNormalized = pdf_RX / areaUnderPDF;

    % Find peak and critical points
    [peak, peakPosition] = max(pdfNormalized);
    P_peak = P(peakPosition);

    aoiCalculation = ones(1, N);
    avgAoI = zeros(1, MC);

    % monte carlo simulation
    for mc = 1:MC
        FindActiveUsers = (aoiUsers >= delta) & (rand(1,N) < activationProbability);
        activeUsers = find(FindActiveUsers);
        receivePowerUsers = zeros(1,N);
        captureSensor = 0;

        if ~isempty(activeUsers)
            for j = 1:length(activeUsers)
                currentActiveUser = activeUsers(j);

                % random number for sampling
                u = rand();
                possibleValues = find(cdfRx < u);
                samplePower = P(possibleValues(end));

                for sat = 1:numSatellites
                    receivePowerUsers(sat, currentActiveUser) = samplePower * (0.95 + 0.1*rand());
                end
            end

            % For each active user, find which satellite receives the strongest signal
            [maxPowerPerUser, bestSatellitePerUser] = max(receivePowerUsers(:, activeUsers), [], 1);

            % capture effect per satellite
            for sat = 1:numSatellites
                userForThisSat = activeUsers(bestSatellitePerUser == sat);

                if ~isempty(userForThisSat)
                    activePower = receivePowerUsers(sat, userForThisSat);
                    [maxPower, positionOfMaxPower] = max(activePower);
                    strongestUser = userForThisSat(positionOfMaxPower);

                    % SINR
                    interference = sum(activePower) - maxPower;
                    SINR = maxPower / (interference + noisePower);

                if 10*log10(SINR) >= 3
                    captureSensor = strongestUser;
                    break;
                end
                end
            end
        end

        % AoI
        aoiUsers = aoiUsers + 1;

        % Reset AoI
        if captureSensor > 0
            aoiUsers(captureSensor) = 0;
        end

        % store results for MC iteration
        avgAoI(mc) = mean(aoiUsers);
    end

    % store final AoI
    finalAoIresult(angleElevation) = mean(avgAoI);
    fprintf('Elevation: %d° completed | Avg AoI: %.3f\n', ElevDeg, finalAoIresult(angleElevation));
end

%Plot AoI vs Elevation angle
figure;
plot(ElevAngle, finalAoIresult, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Elevation Angle (degrees)');
ylabel('Average AoI(in seconds)');
title( 'AoI Performance vs Elevation Angle');


function K = getKfromAngles(angleElevation)
    K_values = [0.9, 1.5, 2.2, 4.1, 8.9, 11.4, 13.5, 15.2, 18.6];
    K_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
    K = K_values(K_angles == angleElevation);
end