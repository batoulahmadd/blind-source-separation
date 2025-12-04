close all; clear all; clc;

rng default

T   = 5;        % periods for each signal
d   = 2;                % mixed observations
r   = 2;                % independent/principal components
numSimulations = 1000;  % Number of Monte Carlo simulations
numBits=1000;
EbNoVec = (-20:20)';    % Eb/No values (dB)

% resetting the values for the correlation functions 
corrfastICA_avg=zeros( 1, length(EbNoVec));
corrAMUSE_avg=zeros( 1, length(EbNoVec));
corrSOBI_avg=zeros( 1, length(EbNoVec));

best_corfastICA=zeros(1,numSimulations);
best_corAMUSE=zeros(1,numSimulations);
best_corSOBI=zeros(1,numSimulations);

% resetting the values for the time functions 
timefastICA_avg=zeros( 1, length(EbNoVec));
timeAMUSE_avg=zeros( 1, length(EbNoVec));
timeSOBI_avg=zeros( 1, length(EbNoVec));

for n = 1:length(EbNoVec)
% Convert Eb/No to SNR
snrdB = EbNoVec(n);
%functions
t = @(numBits,T) linspace(0,1,numBits) * 2 * pi * T; % generates n points. The spacing between the points is (x2-x1)/(n-1).
sigma = @(snrdB,X) exp(-snrdB / 20) * (norm(X(:)) / sqrt(numel(X)));
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
% Generate ground truth
Ztrue=zeros(r, numBits);
Ztrue(1,:)  = sin(t(numBits,T));
Ztrue(2,:)  = sigma(snrdB,Ztrue(1,:)) * randn(size(Ztrue(1,:)));

% resetting the values for the time functions for each snr value
best_timefastICA=zeros(1,numSimulations);
best_timeAMUSE=zeros(1,numSimulations);
best_timeSOBI=zeros(1,numSimulations);

    for i = 1:numSimulations
    
    % Generate mixed signals
    A = normRows(rand(d,2));
    Zmixed = A * Ztrue;
    
    % Perform the three algorithms
    tic;
    datafastICA = fastICA(Zmixed,r);
    best_timefastICA(i)=toc;
    
    tic;
    dataAMUSE = amuse(Zmixed,r);
    best_timeAMUSE(i)=toc;
    
    tic;
    [K, dataSOBI] = sobi(Zmixed,r);
    best_timeSOBI(i)=toc;
    
%         s=1:1000
%         figure();       
%         title('observed and estimated signals');
%         subplot(2,1,1)
%         plot(s,Zmixed(1,:) )
%         grid on;
%         xlabel('simulations')
%         ylabel('observed signal')
%         subplot(2,1,2)
%         plot(s,dataSOBI(2,:))
%         grid on;
%         xlabel('simulations')
%         ylabel('estimated signal')

   
    % Check sizes of the outputs
    numComponentsFastICA = size(datafastICA, 1);
    numComponentsAMUSE = size(dataAMUSE, 1);
    numComponentsSOBI = size(dataSOBI, 1);

   
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(1,:))', (datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(1,:))', (datafastICA(2,:))');   
    end
    best_corfastICA(i) = max(abs(corrfastICA(:, i)));
    
    corrAMUSE(1, i) = corr((Ztrue(1,:))',( dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(1,:))',( dataAMUSE(2,:))');
    end
    best_corAMUSE(i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(1,:)'),( dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(1,:))', (dataSOBI(2,:))');
    end
    best_corSOBI(i) = max(abs(corrSOBI(:, i)));
    
    end

    % Calculate average correlations
    corrfastICA_avg(n)=mean(best_corfastICA,2);
    corrAMUSE_avg(n)=mean(best_corAMUSE,2);
    corrSOBI_avg(n)=mean(best_corSOBI,2);
    
    % calculate avarage time
    timefastICA_avg(n)=mean(best_timefastICA,2);
    timeAMUSE_avg(n)=mean(best_timeAMUSE,2);
    timeSOBI_avg(n)=mean(best_timeSOBI,2);
    
 end
 
        figure(1);
        plot(EbNoVec,corrfastICA_avg, 'LineWidth', 2 )
        hold on;
        plot( EbNoVec,corrAMUSE_avg, 'LineWidth', 1.5)
        hold on
        plot(EbNoVec,corrSOBI_avg, 'LineWidth', 1)
        ylim([0.9 1.1]);
        legend('fastICA','AMUSE','SOBI')
        grid on;
        xlabel('snr')
        ylabel('correlation')
        title('correlation using three ICA algorithm');
        
        figure;
        plot(EbNoVec,timefastICA_avg)
        hold on;
        plot(EbNoVec,timeAMUSE_avg)
        hold on
        plot(EbNoVec,timeSOBI_avg)
        legend('fastICA','Amuse','SOBI')
        grid on;
        xlabel('snr')
        ylabel('running time')
        title('excution time of three ICA algorithms');

        
        

