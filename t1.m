clear; clc; close all;

% Parameters
f_T = 2400e6; % frequency 2.4 GHz
Re = 6378e3; % earth radius
h_s = 400e3; % satellite altitude
h_u = 0; % ground sensor position
G_t = 1; % sensor antenna gain
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
T_sys = 1000;          % System temperature (K)
B = 128e3;             % Bandwidth (Hz)

% number of sensors
N = 10;
activationProbability = 1/N;

% noise
noise_power = k_boltz * T_sys * B;

% monte carlo
MC = 1e4;

% transmit power
Pt = 10;
delta_ = 1;
% elevation angel

ElevAngle = 10:10:90;
finalavgAoIresult = zeros(size(ElevAngle));
P = 0.001:0.001:5;  % Power range for PDF

for angleElevation = 1:length(ElevAngle)
    ElevDeg = ElevAngle(angleElevation);
    ElevRad = deg2rad(ElevDeg);

    aoiUsers = ones(1, N);
    K = getKfromAngles(ElevDeg);  
    K_linear = 10.^(K/10); 
    
    % Calculate Rician PDF
    pdf_RX = K_linear .* exp(-K_linear .* (P + 1)) .* besseli(0, 2 .* K_linear .* sqrt(P));
    
    % Normalize PDF
    areaUnderPDF = trapz(P, pdf_RX);
    pdfNormalized = pdf_RX / areaUnderPDF;
    
    % Find peak and critical points
    [peak, peakPosition] = max(pdfNormalized);
    P_peak = P(peakPosition);
    
    % first point of PDF
    P1 = P(1);
    y1 = pdfNormalized(1);
    
    % point 3 where y value of point 3 match with point 1
    postPeakValue = pdfNormalized(peakPosition+1:end);
    point3postion = find(postPeakValue <= y1, 1, 'first');
    if isempty(point3postion)
        P3 = P(end);
    else
        P3 = P(peakPosition + point3postion);
    end
    
    % regions before and after peak point
    region1 = find(P <= P_peak);       % before peak
    region2 = find(P > P_peak & P <= P3); % after peak until point 3
    
    % region probabilities
    probabilityRegion1 = trapz(P(region1), pdfNormalized(region1));
    probabilityRegion2 = trapz(P(region2), pdfNormalized(region2));
    totalProb = probabilityRegion1 + probabilityRegion2;
    
    % Tolerance for second check
    tolerance = 0.001 * peak;
    
    % Plot PDF and critical points only for first elevation
    if angleElevation == 1
        figure;
        plot(P, pdfNormalized, 'LineWidth', 2);
        hold on;
        plot(P_peak, peak, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        plot(P1, y1, 'go', 'MarkerSize', 10, 'LineWidth', 2);
        plot(P3, y1, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    % Calculate distances & path loss
    d = sqrt((Re + h_s)^2 - (Re + h_u)^2 * cos(ElevRad)^2) - ((Re + h_u) * sin(ElevRad)); 
    pathLoss = (c ./ (4 * pi * f_T * d)).^2;
    
    aoiCalculation = ones(1, N);
    avgAoI = zeros(1, MC);

    % monte carlo simulation
    for mc = 1:MC
        FindActiveUsers = (aoiUsers >= delta_) & (rand(1,N) < activationProbability);
        activeUsers = find(FindActiveUsers);
        receiverPowerUsers = zeros(1,N);
        captureSensor = 0;

        if ~isempty(activeUsers)
            for j = 1:length(activeUsers)
                currentActuveUser = activeUsers(j);
                
                % random number for sampling
                u = rand();

                region = (u <= probabilityRegion1/totalProb) * 1 + (u > probabilityRegion1/totalProb) * 2; % Determine region 1 or 2

                if region == 1
                    cdfRegion = cumtrapz(P(region1), pdfNormalized(region1)) / probabilityRegion1;
                    targetProb = u * totalProb / probabilityRegion1;
                    samplePower = P(region1(find(cdfRegion >= targetProb, 1)));
                else
                    cdfRegion = cumtrapz(P(region2), pdfNormalized(region2)) / probabilityRegion2;
                    targetProb = (u - probabilityRegion1 / totalProb) * totalProb / probabilityRegion2;
                    samplePower = P(region2(find(cdfRegion >= targetProb, 1)));                    
                end

                % Second check
                if samplePower > P_peak && samplePower <= P3
                    pdfValue = interp1(P, pdfNormalized, samplePower);
                    if abs(pdfValue - y1) < tolerance
                        samplePower = P_peak;  % Use peak value
                    end
                end

                % shadowing effect
                shadowing = 10.^(sigma_shadow * randn()/10);
                
                % actual received power
                receiverPowerUsers(currentActuveUser) = Pt * G_t * G_s * Lta * Lra * La * shadowing * pathLoss(currentActuveUser) * samplePower;
            end
        end
    end
end


function K = getKfromAngles(angleElevation)
    K_values = [0.9, 1.5, 2.2, 4.1, 8.9, 11.4, 13.5, 15.2, 18.6];
    K_angles = [10, 20, 30, 40, 50, 60, 70, 80, 90];
    K = K_values(K_angles == angleElevation);
end
