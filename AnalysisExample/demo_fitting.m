%% demonstration of the function fitting
% This program and dataset reproduces maximum-likelihood parameter 
% estimation used in the manuscript.
% The data contains one recording session that had 156 neurons.


%% adding the program folder to PATH.
addpath('Programs');

%% function fitting (only works under macOS)
DoAuditoryAnalysis(pwd);
% this line is for generating 'AuditorySpotSummary_1.mat'
% if you are using Windows or Linux, skip this part.

% this will take ~10 minutes using a laptop with 1.3 GHz core m7 processor.

%% showing example receptive fields
e = Experiment(pwd);
e.resetUseIDs;
examples = [9, 32, 144];

% This will show a raster plot associated with each virtual location at the
% top. and line plots and heat maps in the midle, and the 3D plot that
% represents the RF in spherical coordinate at the bottom.
% IW1-4 indicate different integration windows of time.
% IW1: 5-20 ms (used for the main analysis for the manuscript)
% IW2: 20-100 ms (slow timescale used for Fig. S3g)
% IW3: 105-120 ms (stimulus offset response)
% IW4: 120-200 ms (slow offset response)

% here you would see three neurons.
% All of them have a fast response.
% Neuron 9 and 32 have a nasal RF; neuron 144 has a temporal RF.
% You can go to the next neuron by pressing any button on the keyboard.

for i = 1:length(examples)
    e.plotAuditoryRF(examples(i));
    pause
end

%% showing a map from this recording. (subset of Figure 1h)
% A-P position is only relative in this dataset,
% not completely matched with Fig. 1h.

figure;e.plotAuditoryMap(3);
xlabel('A-P position (µm)');
ylabel('RF azimuth (°)');

