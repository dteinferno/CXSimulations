%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PB-EB integrate and fire model
% Dan Turner-Evans, 5/13/2016
% Based on code by Shaul Druckmann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial assumptions

% 1) The wedge and tile gall neurons get input from the gall when the fly
% moves to raise their activity (let's say rotational velocity + forward
% velocity)
% 2) There's a pre-existing bias in the wedge neurons for the bump to be
% in a predetermined location (stored by magic).
% 3) In the PB, the bump from the wedge neurons excites both populations
% of delta 7 neurons.
% 4) In the PB, delta 7 population #1 excites the gall-tile neurons
% everywhere but at the bump location.
% 5) In the EB, the gall-tile neurons inhibit the wedge neurons everywhere
% except at the bump location to help shape the bump
% 6) In the NO, one or the other side of NO-tile neurons receives uniform
% input (rotational velocity) when the fly is turning in one direction or
% the other.
% 7) In the PB, delta 7 population #2 inhibits the NO-tile neurons
% everywhere except at the bump location.
% 8) In the EB, the rotationally excited, bump uninhibited, shifted
% NO-tile neurons move the bump. 

%% Define connections and time constants

% All connections are done in the PB reference frame, the lower half of
% numbers coming from the left PB and the upper half coming from the right
% PB.

% THINK ABOUT INCLUDING GALL-TIP-WEDGE AND THREE GLOMERULI DELTA 7

% Specify the weights
NO_tile_wedge_weight = 1;
gall_tile_wedge_weight = 0;
wedge_d7_1_weight = 1;
wedge_d7_2_weight = 1;
d7_1_gall_tile_weight = -1/8;
d7_2_NO_tile_weight = -1/8;
inhibit_weight = -0.2;
wedge_inhibit_weight = 1/16;

% specify the time constants
tau_wedge = 1;
tau_NO_tile = 1;
tau_gall_tile = 1;
tau_d7_1 = 1;
tau_d7_2 = 1;
tau_inhibit = 1;

% The number of different neurons to be considered
num_wedges = 16;
wedge_offset = 0;
num_gall_tiles = 16;
gall_tile_offset = num_wedges;
num_NO_tiles = 16;
NO_tile_offset = gall_tile_offset+num_gall_tiles;
num_d7_1 = 8;
d7_1_offset = NO_tile_offset + num_NO_tiles;
num_d7_2 = 8;
d7_2_offset = d7_1_offset + num_d7_1;
num_inhibit = 1;
inhibit_offset = d7_2_offset + num_d7_2;
num_neurons = num_wedges+num_gall_tiles+ num_NO_tiles+num_d7_1+num_d7_2+num_inhibit;

tau_all = vertcat(...
    zeros(num_wedges,1) + tau_wedge,...
    zeros(num_gall_tiles,1) + tau_NO_tile,...
    zeros(num_NO_tiles,1) + tau_gall_tile,...
    zeros(num_d7_1,1) + tau_d7_1,...
    zeros(num_d7_2,1) + tau_d7_2,...
    zeros(num_inhibit,1) + tau_inhibit);


% Create an empty connectivity matrix
J = zeros(num_neurons,num_neurons);


% Define the connections in the EB

% NO-tile to wedge
for Ltile = 1:8
    J(NO_tile_offset+Ltile,wedge_offset+Ltile) = NO_tile_wedge_weight;
    J(NO_tile_offset+Ltile,wedge_offset+16-mod(6+Ltile,8)) = NO_tile_wedge_weight;
end

for Rtile = 9:16
    J(NO_tile_offset+Rtile,wedge_offset+Rtile) = NO_tile_wedge_weight;
    J(NO_tile_offset+Rtile,wedge_offset+8-mod(Rtile-2,8)) = NO_tile_wedge_weight;
end

% gall-tile to wedge
for Ltile = 1:8
    J(gall_tile_offset+Ltile,wedge_offset+Ltile) = gall_tile_wedge_weight;
    J(gall_tile_offset+Ltile,wedge_offset+16-mod(6+Ltile,8)) = gall_tile_wedge_weight;
end

for Rtile = 9:16
    J(gall_tile_offset+Rtile,wedge_offset+Rtile) = gall_tile_wedge_weight;
    J(gall_tile_offset+Rtile,wedge_offset+8-mod(Rtile-2,8)) = gall_tile_wedge_weight;
end

% Define the connections in the PB
% wedge to delta 7 #1 and #2
for LPB=1:8
    J(wedge_offset+LPB,d7_1_offset+LPB) = wedge_d7_1_weight;
    J(wedge_offset+LPB,d7_2_offset+LPB) = wedge_d7_2_weight;
end
for RPB=9:16
    J(wedge_offset+RPB,d7_1_offset+17-RPB) = wedge_d7_1_weight;
    J(wedge_offset+RPB,d7_2_offset+17-RPB) = wedge_d7_2_weight;
end

% Define the delta 7 #1 to NO-tile connections
for PB = 1:8
    J(d7_1_offset+PB,gall_tile_offset+[1:16]) = d7_1_gall_tile_weight;
    J(d7_1_offset+PB,gall_tile_offset+PB) = 0;
    J(d7_1_offset+PB,gall_tile_offset+17-PB) = 0;
    
    J(d7_2_offset+PB,NO_tile_offset+[1:16]) = d7_2_NO_tile_weight;
    if PB == 1
        J(d7_2_offset+PB,NO_tile_offset+8) = 0;
    else
        J(d7_2_offset+PB,NO_tile_offset+mod(PB-1,8)) = 0;
    end
    J(d7_2_offset+PB,NO_tile_offset+16-mod(16+PB,8)) = 0;
end

% Define the inhibitory connections
for wedge = 1:16
    J(inhibit_offset+1,wedge+wedge_offset) = inhibit_weight;
    J(wedge+wedge_offset,inhibit_offset+1) = wedge_inhibit_weight;
end

figure;
imagesc(J)

%% Define the initial conditions and run the model
% Specify the nonlinear function for the input currents to each neuron
fhandle = @ThresholdLinearSaturatingNeuron;
paramStruct = struct;
paramStruct.thresh = 0;
paramStruct.saturation = 10;
paramStruct.tau = tau_all;

% Make the time vector for the simulation
tVec = 0:0.1:20;
timePoints = length(tVec);

% Give the network some initial firing rates
initCond = randn(num_neurons,1);
initCond = zeros(num_neurons,1);
initCond(1:2) = 1;
initCond(15:16) = 1;

% Give the rotational and translational velocity inputs
flyMovement = zeros(num_neurons,1);
% flyMovement(1,1) = 1;
flyMovement(1:num_wedges+num_gall_tiles,1) = 0;
flyMovement(NO_tile_offset+1:NO_tile_offset+8,1) = 0.15;
% flyMovement(NO_tile_offset+15:NO_tile_offset+15,1) = 1;

% Solve the ODE's over this time vector
[t,y] = ode45(@(t, y) CxModelDynamicsForODESolver(t, y, fhandle, paramStruct,J',flyMovement), tVec, initCond);
inputActivity = y';
figure, imagesc(inputActivity)

