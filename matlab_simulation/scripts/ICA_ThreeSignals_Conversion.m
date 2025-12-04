%% ICA demo (signals)
close all; clear all; clc;

rng default

T   = [5, 4, 3];        % periods for each signal
d   = 3;                % mixed observations
r   = 3;                % independent/principal components
numSimulations = 100;  % Number of Monte Carlo simulations
numBits=1000;
EbNoVec = (-2:3)';    % Eb/No values (dB)

%% correlation criterion for signals with noise

% resetting the values for the correlation functions 
corrfastICA_avg=zeros( r, length(EbNoVec));
corrAMUSE_avg=zeros( r, length(EbNoVec));
corrSOBI_avg=zeros( r, length(EbNoVec));

best_corfastICA=zeros(r,numSimulations);
best_corAMUSE=zeros(r,numSimulations);
best_corSOBI=zeros(r,numSimulations);

for n = 1:length(EbNoVec)
% Convert Eb/No to SNR
snrdB = EbNoVec(n);
%functions
t = @(numBits,T) linspace(0,1,numBits) * 2 * pi * T; % generates n points. The spacing between the points is (x2-x1)/(n-1).
sigma = @(snrdB,X) exp(-snrdB / 20) * (norm(X(:)) / sqrt(numel(X))); %scales ,  magnitude of the matrix X
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
% Generate ground truth
Ztrue=zeros(r, numBits);
Ztrue(1,:)  = sin(t(numBits,T(1)));
Ztrue(2,:)  = sign(sin(t(numBits,T(2))));
Ztrue(3,:) = sawtooth(t(numBits,T(3))); 

    for i = 1:numSimulations
    
    % Add some noise to make the signals "look" interesting
    Ztrue_with_noise = Ztrue + sigma(snrdB,Ztrue) * randn(size(Ztrue));
    % Generate mixed signals
    A = normRows(rand(d,3));
    Zmixed = A * Ztrue_with_noise;        
    
        
%         s=1:1000
%         figure();
%         plot(s,Zmixed(1,:) )
%         hold on;
%         plot(s,Zmixed(2,:) )
%         hold on
%         plot(s,Zmixed(3,:) )
%         legend('Sinusoid','Square','Sawtooth')
%         grid on;
%         xlabel('simulations')
%         ylabel('observed signals')
%         title('observed signals in case of three signals with noise');
    % Perform the three algorithms
    datafastICA = fastICA(Zmixed,r);
    dataAMUSE = amuse(Zmixed,3);
    [K, dataSOBI] = sobi(Zmixed,r);
   
    % Check sizes of the outputs
    numComponentsFastICA = size(datafastICA, 1);
    numComponentsAMUSE = size(dataAMUSE, 1);
    numComponentsSOBI = size(dataSOBI, 1);

    
    % For a Sinusoid signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(1,:))', (datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(1,:))', (datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(1,:))', (datafastICA(3,:))');    
    end
    best_corfastICA(1, i) = max(abs(corrfastICA(:, i)));
    
    corrAMUSE(1, i) = corr((Ztrue(1,:))',( dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(1,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(1,:))',( dataAMUSE(3,:))');
    end
    best_corAMUSE(1, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(1,:)'),( dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(1,:))', (dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(1,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(1, i) = max(abs(corrSOBI(:, i)));
    
    % Repeat the above for the Square and Sawtooth signals
    
    % For a Square signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(2,:))',( datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(2,:))',(datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(2,:))',( datafastICA(3,:))');    
    end
    best_corfastICA(2, i) = max(abs(corrfastICA(:, i)));
     
    corrAMUSE(1, i) = corr((Ztrue(2,:))', (dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(2,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(2,:))',( dataAMUSE(3,:))');
    end
    best_corAMUSE(2, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(2,:)'), (dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(2,:))',( dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(2,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(2, i) = max(abs(corrSOBI(:, i)));
     
    % For a Sawtooth signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(3,:))', (datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(3,:))', (datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(3,:))', (datafastICA(3,:))');    
    end
    best_corfastICA(3, i) = max(abs(corrfastICA(:, i)));
    
    corrAMUSE(1, i) = corr((Ztrue(3,:))',( dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(3,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(3,:))', (dataAMUSE(3,:))');
    end
    best_corAMUSE(3, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(3,:)'), (dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(3,:))', (dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(3,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(3, i) = max(abs(corrSOBI(:, i)));
    
    end

    % Calculate average correlations
    corrfastICA_avg( :, n)=mean(best_corfastICA,2);
    corrAMUSE_avg( :, n)=mean(best_corAMUSE,2);
    corrSOBI_avg( :, n)=mean(best_corSOBI,2);
    
 end
 
        figure(1);
        plot(EbNoVec,corrfastICA_avg(1,:) )
        hold on;
        plot( EbNoVec,corrfastICA_avg(2,:))
        hold on
        plot(EbNoVec,corrfastICA_avg(3,:))
        legend('Sinusoid','Square','Sawtooth')
        grid on;
        xlabel('snr')
        ylabel('correlation')
        title('correlation using ICA algorithm fastICA on three diffrent signals');

        
        
        figure(2);
        plot(EbNoVec,corrAMUSE_avg(1,:) )
        hold on;
        plot( EbNoVec,corrAMUSE_avg(2,:))
        hold on
        plot(EbNoVec,corrAMUSE_avg(3,:))
        legend('Sinusoid','Square','Sawtooth')
        grid on;
        xlabel('snr')
        ylabel('correlation')
        title('correlation using ICA algorithm AMUSE on three diffrent signals');
        
        
        
        figure(3);
        plot(EbNoVec,corrSOBI_avg(1,:) )
        hold on;
        plot( EbNoVec,corrSOBI_avg(2,:))
        hold on
        plot(EbNoVec,corrSOBI_avg(3,:))
        legend('Sinusoid','Square','Sawtooth')
        grid on;
        xlabel('snr')
        ylabel('correlation')
        title('correlation using ICA algorithm SOBI on three diffrent signals');
        
%% time criterion for signals with noise

% resetting the values for the time functions 
timefastICA_avg=zeros( 1, length(EbNoVec));
timeAMUSE_avg=zeros( 1, length(EbNoVec));
timeSOBI_avg=zeros( 1, length(EbNoVec));

for n = 1:length(EbNoVec)
% Convert Eb/No to SNR
snrdB = EbNoVec(n);
%functions
t = @(numBits,T) linspace(0,1,numBits) * 2 * pi * T;
sigma = @(snrdB,X) exp(-snrdB / 20) * (norm(X(:)) / sqrt(numel(X)));
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
% Generate ground truth
Ztrue=zeros(r, numBits);
Ztrue(1,:)  = sin(t(numBits,T(1)));
Ztrue(2,:)  = sign(sin(t(numBits,T(2))));
Ztrue(3,:) = sawtooth(t(numBits,T(3)));

% resetting the values for the time functions for each snr value
best_timefastICA=zeros(1,numSimulations);
best_timeAMUSE=zeros(1,numSimulations);
best_timeSOBI=zeros(1,numSimulations);

    for i = 1:numSimulations
    % Add some noise to make the signals "look" interesting
    Ztrue_with_noise = Ztrue + sigma(snrdB,Ztrue) * randn(size(Ztrue));
    % Generate mixed signals
    A = normRows(rand(d,3));
    Zmixed = A * Ztrue_with_noise;
            s=1:1000
        figure();       
        title('observed and estimated signals');
        subplot(2,1,1)
        plot(s,Zmixed(1,:) )
        hold on;
        plot(s,Zmixed(2,:) )
        hold on;
        plot(s,Zmixed(3,:) )
        grid on;
        xlabel('simulations')
        ylabel('observed signal')
        subplot(2,1,2)
        plot(s,dataSOBI(2,:))
        grid on;
        xlabel('simulations')
        ylabel('estimated signal')

    
    % Perform the three algorithms
    % computing running time for each algorithm
    tic;
    datafastICA = fastICA(Zmixed,r);
    best_timefastICA(i)=toc;
    
    tic;
    dataAMUSE = amuse(Zmixed,3);
    best_timeAMUSE(i)=toc;
    
    tic;
    [K, dataSOBI] = sobi(Zmixed,r);
    best_timeSOBI(i)=toc;
    
    end

    % Calculate average running time
    timefastICA_avg(n)=mean(best_timefastICA,2);
    timeAMUSE_avg(n)=mean(best_timeAMUSE,2);
    timeSOBI_avg(n)=mean(best_timeSOBI,2);
    
end

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
      
        
%% time and correlation criterion for signals without noise

% resetting the values for the correlation functions 
% for two specific snr values

corrfastICA_avg=zeros( r, 2);
corrAMUSE_avg=zeros( r, 2);
corrSOBI_avg=zeros( r, 2);

best_corfastICA=zeros(r,numSimulations);
best_corAMUSE=zeros(r,numSimulations);
best_corSOBI=zeros(r,numSimulations);

% resetting the values for the time functions 
timefastICA_avg=zeros( 1, 2);
timeAMUSE_avg=zeros( 1, 2);
timeSOBI_avg=zeros( 1, 2);


for q=1:2
    
% for two specific snr values
if q==1
    n=6;
else n=36;
end

% Convert Eb/No to SNR
snrdB = EbNoVec(n);
%functions
t = @(numBits,T) linspace(0,1,numBits) * 2 * pi * T;
sigma = @(snrdB,X) exp(-snrdB / 20) * (norm(X(:)) / sqrt(numel(X)));
normRows = @(X) bsxfun(@rdivide,X,sum(X,2));
% Generate ground truth
Ztrue=zeros(r, numBits);
Ztrue(1,:)  = sin(t(numBits,T(1)));
Ztrue(2,:)  = sign(sin(t(numBits,T(2))));
Ztrue(3,:) = sawtooth(t(numBits,T(3))); 

% resetting the values for the time functions for each snr value
best_timefastICA=zeros(1,numSimulations);
best_timeAMUSE=zeros(1,numSimulations);
best_timeSOBI=zeros(1,numSimulations);

    for i = 1:numSimulations
    
    % Add some noise to make the signals "look" interesting
    Ztrue_with_noise = Ztrue + sigma(snrdB,Ztrue) * randn(size(Ztrue));
    % Generate mixed signals
    A = normRows(rand(d,3));
    Zmixed = A * Ztrue;

    % computing running time for each algorithm
    tic;
    datafastICA = fastICA(Zmixed,r);
    best_timefastICA(i)=toc;
    
    tic;
    dataAMUSE = amuse(Zmixed,3);
    best_timeAMUSE(i)=toc;
    
    tic;
    [K, dataSOBI] = sobi(Zmixed,r);
    best_timeSOBI(i)=toc;
   
    % Check sizes of the outputs
    numComponentsFastICA = size(datafastICA, 1);
    numComponentsAMUSE = size(dataAMUSE, 1);
    numComponentsSOBI = size(dataSOBI, 1);

    % For a Sinusoid signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(1,:))', (datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(1,:))', (datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(1,:))', (datafastICA(3,:))');    
    end
    best_corfastICA(1, i) = max(abs(corrfastICA(:, i)));
    
    corrAMUSE(1, i) = corr((Ztrue(1,:))',( dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(1,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(1,:))',( dataAMUSE(3,:))');
    end
    best_corAMUSE(1, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(1,:)'),( dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(1,:))', (dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(1,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(1, i) = max(abs(corrSOBI(:, i)));
    
    % Repeat the above for the Square and Sawtooth signals
    
    % For a Square signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(2,:))',( datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(2,:))',(datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(2,:))',( datafastICA(3,:))');    
    end
    best_corfastICA(2, i) = max(abs(corrfastICA(:, i)));
     
    corrAMUSE(1, i) = corr((Ztrue(2,:))', (dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(2,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(2,:))',( dataAMUSE(3,:))');
    end
    best_corAMUSE(2, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(2,:)'), (dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(2,:))',( dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(2,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(2, i) = max(abs(corrSOBI(:, i)));
     
    % For a Sawtooth signal
    % Reset the correlation coefficients 
    corrfastICA=zeros(r,numSimulations);
    corrAMUSE=zeros(r,numSimulations);
    corrSOBI=zeros(r,numSimulations);
    %compute the correlation between the original and the estimated signals
    corrfastICA(1, i) = corr((Ztrue(3,:))', (datafastICA(1,:))');
    if numComponentsFastICA >= 2
        corrfastICA(2, i) = corr((Ztrue(3,:))', (datafastICA(2,:))');   
    end
    if numComponentsFastICA >= 3
        corrfastICA(3, i) = corr((Ztrue(3,:))', (datafastICA(3,:))');    
    end
    best_corfastICA(3, i) = max(abs(corrfastICA(:, i)));
    
    corrAMUSE(1, i) = corr((Ztrue(3,:))',( dataAMUSE(1,:))');
    if numComponentsAMUSE >= 2
        corrAMUSE(2, i) = corr((Ztrue(3,:))',( dataAMUSE(2,:))');
    end
    if numComponentsAMUSE >= 3
        corrAMUSE(3, i) = corr((Ztrue(3,:))', (dataAMUSE(3,:))');
    end
    best_corAMUSE(3, i) = max(abs(corrAMUSE(:, i)));
    
    corrSOBI(1, i) = corr((Ztrue(3,:)'), (dataSOBI(1,:))');
    if numComponentsSOBI >= 2
        corrSOBI(2, i) = corr((Ztrue(3,:))', (dataSOBI(2,:))');
    end
    if numComponentsSOBI >= 3
        corrSOBI(3, i) = corr((Ztrue(3,:))', (dataSOBI(3,:))');
    end
    best_corSOBI(3, i) = max(abs(corrSOBI(:, i)));
    
    end

    % Calculate average correlations
    corrfastICA_avg( :,q)=mean(best_corfastICA,2)
    corrAMUSE_avg( :,q)=mean(best_corAMUSE,2)
    corrSOBI_avg( :, q)=mean(best_corSOBI,2)
    % calculate avarage time
    timefastICA_avg(q)=mean(best_timefastICA,2)
    timeAMUSE_avg(q)=mean(best_timeAMUSE,2)
    timeSOBI_avg(q)=mean(best_timeSOBI,2)
    
end
