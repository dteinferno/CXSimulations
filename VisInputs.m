%% Model of the relevant visual inputs for bump movement

% input = movie of the visual scene over time

% output = plot of the dominant visual direction over time

% internal parameters: the RFs of the ring neurons specified as a three
% dimensional kernel (incorporating the time dependence)

% calculation = multiple the kernels by the visual display and sum them to
% get:
% 1) The visual area of maximum activity at any point in time
% 2) The most active glomerulus at any point in time

%% First generate the RFs

% The size of the visual scene
xPx = 360; % pixels
yPx = 120; % pixels
px2deg = 1; % pixel to degree conversion

% The center of the RF overlap and their spread
xRFcent = 65; % degrees
yRFcent = 0; % degrees
xRFspread = 30; % degrees
yRFspread = 15; % degrees

% The center and size of the contralateral inhibitory region
xIncent = 40; % degrees
radIn = 40; % degrees
yIncent = 0; % degrees

% The radial size of the RF and contralateral inhibitory region
RFrad = 30; % degrees

% The number of RFs per side
numRFs = 40;

% The time constants of the exhitatory and inhibitory regions
tauEx = 0.5; % sec
tauIn = 0.5; %3 sec
tStep = movieRate; % sec

% The strength of the excitatory and inhibitory peaks
sEx = 2;
sIn = -1;
sInCont = -0.25;

allRFs = zeros(numRFs*2, xPx, yPx, round(max(tauEx,tauIn)*2/tStep));

% Generate the RF kernel at t = 0 for one side of the visual field
for RF = 1:numRFs
    % add some variance to the RF center and radius
    RFcent = [xRFcent+randn*xRFspread yRFcent + randn*yRFspread];
    RFradNow = RFrad*(1+randn/3);
    
    if RFcent(1)-RFradNow < 0
        RFcent(1) = RFradNow;
    end
    if RFcent(1) + RFradNow > xPx/2*px2deg
        RFcent(1) = xPx/2*px2deg-RFradNow-1;
    end
    if RFcent(2)-RFradNow < -yPx/2*px2deg
        RFcent(2) = -yPx/2*px2deg+RFradNow+1;
    end
    if RFcent(2)+RFradNow > yPx/2*px2deg
        RFcent(2) = yPx/2*px2deg-RFradNow-1;
    end

    % create the local exhitatory region (assume that 1 pixel = 1 degree,
    % square fields, and that the exhitatory portion is 1/3 of the 
    % horizontal field, for now)
    allRFs(RF,...
        round(xPx/2-RFcent(1)-RFradNow/3):round(xPx/2-RFcent(1)+RFradNow/3),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sEx;
    
    % create the farther out inhibitory regions (see above assumptions
    allRFs(RF,...
        round(xPx/2-RFcent(1)-RFradNow):round(xPx/2-RFcent(1)-RFradNow/3),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sIn;
    allRFs(RF,...
        round(xPx/2-RFcent(1)+RFradNow/3):round(xPx/2-RFcent(1)+RFradNow),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sIn;
    
    for xCirc = 1:xPx
        for yCirc = 1:yPx
            if ((xCirc-xPx/2+RFcent(1))/ RFradNow)^2 + ((yCirc-yPx/2-RFcent(2)) /RFradNow)^2 > 1
                allRFs(RF, xCirc, yCirc, 1) = 0;
            end
        end
    end
    
    % create the contralateral inhibitory region (assume that 1 pixel = 1
    % degree and square fields, for now)
    for xCirc = 1:xPx
        for yCirc = 1:yPx
            if ((xCirc-xPx/2-xIncent)/ radIn)^2 + ((yCirc-yPx/2-yIncent) /radIn)^2 <= 1
                allRFs(RF, xCirc, yCirc, 1) = sInCont;
            end
        end
    end
end

% Generate the RF kernel at t = 0 for the other side of the visual field
for RF = numRFs+1:numRFs*2
    % add some variance to the RF center and radius
    RFcent = [xRFcent+randn*xRFspread yRFcent + randn*yRFspread];
    RFradNow = RFrad*(1+randn/3);
    
    if RFcent(1)-RFradNow < 0
        RFcent(1) = RFradNow;
    end
    if RFcent(1) + RFradNow > xPx/2*px2deg
        RFcent(1) = xPx/2*px2deg-RFradNow-1;
    end
    if RFcent(2)-RFradNow < -yPx/2*px2deg
        RFcent(2) = -yPx/2*px2deg+RFradNow+1;
    end
    if RFcent(2)+RFradNow > yPx/2*px2deg
        RFcent(2) = yPx/2*px2deg-RFradNow-1;
    end

    % create the local exhitatory region (assume that 1 pixel = 1 degree,
    % square fields, and that the exhitatory portion is 1/3 of the 
    % horizontal field, for now)
    allRFs(RF,...
        round(xPx/2+RFcent(1)-RFradNow/3):round(xPx/2+RFcent(1)+RFradNow/3),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sEx;
    
    % create the farther out inhibitory regions (see above assumptions
    allRFs(RF,...
        round(xPx/2+RFcent(1)-RFradNow):round(xPx/2+RFcent(1)-RFradNow/3),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sIn;
    allRFs(RF,...
        round(xPx/2+RFcent(1)+RFradNow/3):round(xPx/2+RFcent(1)+RFradNow),...
        round(yPx/2+RFcent(2)-RFradNow):round(yPx/2+RFcent(2)+RFradNow),...
        1) = sIn;
    
    for xCirc = 1:xPx
        for yCirc = 1:yPx
            if ((xCirc-xPx/2-RFcent(1))/ RFradNow)^2 + ((yCirc-yPx/2-RFcent(2)) /RFradNow)^2 > 1
                allRFs(RF, xCirc, yCirc, 1) = 0;
            end
        end
    end
    
    % create the contralateral inhibitory region (assume that 1 pixel = 1
    % degree and square fields, for now)
    for xCirc = 1:xPx
        for yCirc = 1:yPx
            if ((xCirc-xPx/2+xIncent)/ radIn)^2 + ((yCirc-yPx/2+yIncent) /radIn)^2 <= 1
                allRFs(RF, xCirc, yCirc, 1) = sInCont;
            end
        end
    end
end

% Propagate the kernel over time
for RF = 1:numRFs*2
    % Make masks for the excitatory and inhibitory regions
    RF0 = squeeze(allRFs(RF,:,:,1));
    exMask = 0.5*(RF0+abs(RF0));
    inMask = 0.5*(RF0-abs(RF0));
    for tIt = 1:round(max(tauEx,tauIn)*2/tStep)
        % propagate for the excitatory stages
        allRFs(RF,:,:,tIt) = exMask.*exp(-tIt*tStep/tauEx) + inMask.*exp(-tIt*tStep/tauIn);
    end
end

%% Now convolute the RFs with the movie

% Calculate the mean RF responses and the ROI with the max response
% meanResp = zeros(xPx,yPx,size(movie,3)-size(allRFs,4)+1);
RFMax = zeros(size(movie,3)-size(allRFs,4)+1,1);
RFMaxVals = zeros(size(allRFs,1),size(movie,3)-size(allRFs,4)+1);

for i = size(allRFs,4):length(RFMax)
    i
    tempRFcomp = zeros(size(allRFs,1), size(allRFs,2), size(allRFs,3));
    for RF = 1:size(allRFs,1)
        for tIt = 1:size(allRFs,4)
            tempRFcomp(RF,:,:) = squeeze(tempRFcomp(RF,:,:)) + squeeze(allRFs(RF,:,:,tIt)).*squeeze(movie(:,:,i+tIt-1));
        end
    end
%     meanResp(:,:,i) = squeeze(mean(tempRFcomp,1));
    tempRFcompMean = squeeze(mean(mean(tempRFcomp,2),3));
    RFmaxVals(:,i) = tempRFcompMean;
    if max(tempRFcompMean) <= 0
        RFMax(i) = 0;
    else
        fndMaxRFs = find(tempRFcompMean == max(tempRFcompMean));
        RFMax(i) = fndMaxRFs(1);
    end
end

%% Plot the maximally excited RF center vs. the feature positions
ExAng = zeros(length(RFMax),1);
for i=1:length(RFMax)
    if RFMax(i) > 0
        exRFNow = 0.5*(squeeze(allRFs(RFMax(i),:,:,1))+abs(squeeze(allRFs(RFMax(i),:,:,1))));
        ExAng(i) = mean(find(mean(exRFNow,2)))-180;
    end
end

OffsetRotUnwrap = UnWrap(OffsetRot,2,0);
figure;
plot(t(1:40:length(ExAng)*40),ExAng*pi/180,'k','LineWidth',2);
hold on;
plot(t,(mod(OffsetRotUnwrap-25,360)-180)*pi/180,'r','LineWidth',2);
plot(t,(mod(OffsetRotUnwrap+95,360)-180)*pi/180,'b','LineWidth',2);
plot(t,(mod(OffsetRotUnwrap-115,360)-180)*pi/180,'g','LineWidth',2);

set(gca,'FontSize',16);
xlabel('Time (sec)');
ylabel('Obj. Orientation (rad)');
xlim([25 190]);
ylim([-pi pi]);

rectangle('Position', [t(1) -pi t(end)-t(1) pi/3]);
rectangle('Position', [t(1) 2*pi/3 t(end)-t(1) pi/3]);

RFcentVals = zeros(size(allRFs,1),1);
for i=1:length(RFcentVals)
    exRFNow = 0.5*(squeeze(allRFs(i,:,:,1))+abs(squeeze(allRFs(i,:,:,1))));
    RFcentVals(i) = mean(find(mean(exRFNow,2)))-180;
end

angs = [-180:1:180];
allMaxs = zeros(length(angs),length(RFmaxVals));

for i=1:length(RFcentVals)
    allMaxs(round(RFcentVals(i))+180,:) = squeeze(RFmaxVals(i,:));
end

figure;
imagesc(t(1:40:length(ExAng)*40),pi/180*angs,flipud(allMaxs));
caxis([0 max(max(allMaxs))]);
xlim([25 190]);
ylim([-pi pi]);
set(gca,'FontSize',16);
xlabel('Time (sec)');
ylabel('Obj. Orientation (rad)');

%% Plot a few representative RFs and the sum of all RFs
figure;

for i=1:6
    subplot(4,3,i)
    RF2plt = round(rand*numRFs*2);
    imagesc([1:360]-180,[1:120]-60,squeeze(allRFs(RF2plt,:,:,1))')
    axis equal;
    xlim([-180 180]);
    ylim([-60 60]);
    xlabel('Visual field (deg)');
end

subplot(4,3,7:12)
imagesc([1:360]-180,[1:120]-60,squeeze(mean(allRFs(:,:,:,1),1))')
axis equal;
xlim([-180 180]);
ylim([-60 60]);
xlabel('Visual field (deg)');

%% Show the movie template
figure;
imagesc([1:360]-180,[1:120]-60,flipud(movTemplate'));
axis equal;
xlim([-180 180]);
ylim([-60 60]);
xlabel('Visual field (deg)');
set(gca,'FontSize',16);
colormap('gray');

%% Show the fly's rotation
OffsetRotUnwrap = UnWrap(OffsetRot,2,0);
figure;
plot(t,OffsetRotUnwrap*pi/180,'k','LineWidth',2);

set(gca,'FontSize',16);
xlabel('Time (sec)');
ylabel('Obj. Orientation (rad)');
xlim([25 190]);
