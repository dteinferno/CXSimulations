%% Generate a fake movie with three objects

%% Create a movie template of three stripes
movTemplate = zeros(round(movieResX*1.5),movieResY);
stripeWidth = 15;
stripeHeight = 30;
movTemplate(51:51+stripeWidth,movieResY/2:movieResY/2+stripeHeight) = 1;
movTemplate(50+121:50+121+stripeWidth,movieResY/2:movieResY/2+stripeHeight) = 1;
movTemplate(170+91:170+91+stripeWidth,movieResY/2:movieResY/2+stripeHeight) = 1;

%% Create a movie template of three, distinct objects
% Set the image size and generate the blank matrix
yPx = 120;
xPx = yPx * 3;
texIm = zeros(xPx,yPx);

% Specify the centroids of the objects
numObjs = 3;
objSpac1 = xPx*1/8;
objSpac2 = xPx*1/3;
objSpac3 = xPx*5/12;
bottomD = 70;

% Draw a diamond at the first point
dWidth = 20;
texIm(objSpac1-dWidth/2:objSpac1+dWidth/2,bottomD) = 255;
for dSide = 1:dWidth-1
    texIm(objSpac1-dWidth/2+dSide:objSpac1+dWidth/2-dSide,bottomD+dSide) = 255;
    texIm(objSpac1-dWidth/2+dSide:objSpac1+dWidth/2-dSide,bottomD-dSide) = 255;
end

% Draw a triangle at the second point
tAng = pi/4;
botT = round(2*dWidth/sqrt(2*tan(tAng)));
hT = round(dWidth*sqrt(tan(tAng)/2));
texIm((objSpac1+objSpac2)-1:(objSpac1+objSpac2)+1,bottomD-round(hT*2/3)) = 255;
for tSide = 1:hT
   texIm((objSpac1+objSpac2)-round(tSide/tan(tAng)):(objSpac1+objSpac2)+round(tSide/tan(tAng)),bottomD-round(hT*2/3)+tSide) = 255;
end

% Draw a star at the third point
sqSide = 8;
tH = round((dWidth^2/2-sqSide^2)/(2*sqSide));
for sqDraw = 1:sqSide+1
    texIm((objSpac1+objSpac2+objSpac3)-sqSide/2:(objSpac1+objSpac2+objSpac3)+sqSide/2,bottomD+sqSide/2+1-sqDraw) = 255;
end
for tDraw = 1:tH
    texIm((objSpac1+objSpac2+objSpac3)-round(tDraw*sqSide/tH/2):(objSpac1+objSpac2+objSpac3)+round(tDraw*sqSide/tH/2),bottomD+sqSide/2+tH+1-tDraw) = 255;
    texIm((objSpac1+objSpac2+objSpac3)-round(tDraw*sqSide/tH/2):(objSpac1+objSpac2+objSpac3)+round(tDraw*sqSide/tH/2),bottomD-sqSide/2-tH-1+tDraw) = 255;
    texIm((objSpac1+objSpac2+objSpac3)+sqSide/2+tH+1-tDraw,bottomD-round(tDraw*sqSide/tH/2):bottomD+round(tDraw*sqSide/tH/2)) = 255;
    texIm((objSpac1+objSpac2+objSpac3)-sqSide/2-tH-1+tDraw,bottomD-round(tDraw*sqSide/tH/2):bottomD+round(tDraw*sqSide/tH/2)) = 255;
end

movTemplate = fliplr(circshift(texIm,140)')';

%% Generate the movie
% Load a real fly's trajectory
fName = 'D:\Imaging\SS131\3Obj\20160311\Fly1_4-6day_6fxSS131_04.TXT';
fileID = fopen(fName);
tstamp = fgetl(fileID);
exTypeStr = strsplit(tstamp,'_');
exType = exTypeStr{end}(1:end-4);
formatSpec = '%s %f %s %f %s %f %s %f %s %f %s %f %s %d %s %d %s %d %s %d %s %d %s %f';
N=400000;
C = textscan(fileID,formatSpec,N,'CommentStyle','Current','Delimiter','\t');
t = C{1,2}; % Time
OffsetRot = C{1,4}; % Stripe rotational offset
OffsetRot = mod(OffsetRot+180, 360)-180;
OffsetFor = C{1,6}; % Stripe forward offset
OffsetLat = C{1,8}; % Stripe lateral offset
dx0 = C{1,10}; % X position of the ball from camera 1 
dx1 = C{1,12}; % X position of the ball from camera 2
dy0 = C{1,14};
dy1 = C{1,16};
closed = C{1,18};
direction = C{1,20};
trans = C{1,22};
gain = C{1,24};
fclose(fileID);

% specify temporal parameters of the movie
movieRate = mean(diff(t))*40; % sec
movieResX = 240; % Px
movieResY = 120; % Px

movie = zeros(movieResX,movieResY,round(length(t)/40));

for i=1:floor(length(t)/40)
    frameNow = circshift(movTemplate,round(OffsetRot(i*40)));
    movie(:,:,i) = frameNow(1:movieResX,:);
end

%% Show the fake movie
figure;
for i=1:length(movie)
    imagesc(squeeze(movie(:,:,i))');
    axis equal;
    text(0,-5,num2str(i))
    drawnow;
end