%% Load and process data
tic
% close all
% clear variables
clc

% FILL IN CAMERA NRs, TRIAL NR AND DIRECTION, AND PLAYER ID
% camNums = [4 3];            % [1 2] for agility FCT or [4 3] for arrowhead
trialName = 'arrowhead';    %arrowhead or agilityFCT
% trialDir = 'L';             %L for left or R for right
% subjectNum = 1;            %Player ID (PID)


% PID = PID';

%%
% f = figure;
% f.WindowState = "maximized";
% ax1 = subplot(2,1,1);
% hold on
% ylim([0 12])
% xlabel('time (s)','FontSize',16)
% ylabel('Velocity (m/s)','FontSize',16)
% yline(9,'m--','LineWidth', 2)
% title('Velocity of all players over time','FontSize',16)
% ax = gca;
% ax.FontSize = 17; 
% 
% ax2 = subplot(2,1,2);
% hold on
% ylim([-10 10])
% xlabel('time (s)','FontSize',16)
% ylabel('Acceleration (m/s^2)','FontSize',16)
% yline(6,'m--','LineWidth', 2)
% yline(-8,'m--','LineWidth', 2)
% title('Acceleration of all players over time','FontSize',16)
% ax = gca;
% ax.FontSize = 17; 
%
%LEFT
% PID = 1:50; PID(33:34) = []; %PID 33 and 34 have not been recorded
% for pIdx = 1:48 %11:32 %

%RIGHT 30 RIGHT IS ALSO UNUSABLE DUE TO RAIN ON LENS
% PID = 1:50; PID(1:11,33:34) = [];% PID(33:34) = []; %PID 33 and 34 have not been recorded
for dir = 1:2
    if dir == 1 %LEFT  %strcmp(trialDir,'L') == 1
        trialDir = 'L';
        PID = 1:50; 
        PID(PID==26) = []; %druppel op de lens
        PID(PID==27) = []; %druppel op de lens
        PID(PID==33) = []; %niet gemeten wegens technische problemen
        PID(PID==34) = []; %niet gemeten wegens technische problemen
    elseif dir == 2 %RIGHT  %strcmp(trialDir,'R') == 1
        trialDir = 'R';
        PID = 11:50;
        PID(PID==25) = []; %speler buiten beeld door te ruime bocht
        PID(PID==30) = []; %druppel op de lens
        PID(PID==31) = []; %speler buiten beeld door te ruime bocht
        PID(PID==32) = []; %speler buiten beeld door te ruime bocht
        PID(PID==33) = []; %niet gemeten wegens technische problemen
        PID(PID==34) = []; %niet gemeten wegens technische problemen
    end
% 
% PID 25 L is traag maar wel heel snel bij tweede bocht?
for pIdx = 1:size(PID,2) %11:32 %
% % subjectNum = 47;
% 
% try
subjectNum = PID(pIdx);
disp(subjectNum)
% for subjectNum = 1:1
clearvars -except trialName trialDir subjectNum PID f
if strcmp(trialName,'arrowhead') == 1 %trialName == 'arrowhead'
    trialAbbrev = 'AH';
    if subjectNum >= 1 && subjectNum <= 10 %monday morning
        load('stereoCalibrationSessionMondayMorningAH.mat');
        camNums = [3 4];
        % lineColor = 'm';
    elseif subjectNum >= 11 && subjectNum <= 21 %monday afternoon
        load('stereoCalibrationSessionMondayAfternoonAH2.mat');
        camNums = [4 3];
        % lineColor = 'c';
    elseif subjectNum > 21
        load('stereoCalibrationSessionTuesdayAH2.mat');
        camNums = [4 3];
        % lineColor = 'g';
    end
elseif strcmp(trialName,'agilityFCT') == 1 %trialName == 'agilityFCT'
    trialAbbrev = 'AFCT';
    camNums = [1 2];
    if subjectNum >= 1 && subjectNum <= 21 %monday
        load('stereoCalibrationSessionMondayAFCT.mat');
    elseif subjectNum > 21
        load('stereoCalibrationSessionTuesdayAFCT.mat');
    end
end
stereoCalibrationSession = calibrationSession; clearvars calibrationSession;

vFile = strcat(['D:\GP\Output\output_videos\testing-days\',trialName,'\PID-',sprintf('%02d', subjectNum),'-' ...
    ,trialAbbrev,'-',trialDir,'-cam',num2str(camNums(1)),'.avi']);
vid = VideoReader(vFile);
% vid_duration = vid.Duration;

% if subjectNum >= 1 && subjectNum <= 21 %if monday
%     if trialName == 1 %if arrowhead

% strcat('D:\GP\Output\output_jsons\testing-days\',test,'\cam',num2str(camNums(k)),'\PID-',num2str(subjectNum),'-AH-',trialDir)

% For all camera's k, do:
for k = 1:length(camNums)
   initCals(:,k) = load(strcat('initCalibrationSessionCam',num2str(camNums(k)),'.mat'));

   % Load specified camera and trial data into workspace
   % cameraParams(:,k) = load(strcat('initialCalibrationSessionCam',num2str(camNums(k)),'.mat'));
   inputFiles(:,k) = dir(strcat('D:\GP\Output\output_jsons\testing-days\',trialName,'\cam',num2str(camNums(k)),'\PID-',sprintf('%02d', subjectNum),'-',trialAbbrev,'-',trialDir));
   
   % For each JSON file (containing the pose keypoints within one frame of a video), load the path and file names
   for i = 3:length(inputFiles)
       jsonName{i-2,k} = inputFiles(i,k).name;
       path{i-2,k} = inputFiles(i,k).folder;
   end

   % For each JSON file, load and decode the data into useable variables
   for i = 1:length(jsonName)
       str{i,k} = fileread(fullfile(path{i,k},jsonName{i,k}));
       data{i,k} = jsondecode(str{i,k});
       try
           pos{i,k} = data{i,k}.people.pose_keypoints_2d;
       catch
           continue
       end
   end
end

% for i = 1:size(pos,1)
%     conf1(i,:) = pos{i,1}(3:3:end);
% end
% 
% conf1(:,1:5) = [];
% 
% figure;
% hold on;
% xlim([0 size(conf1,1)])
% for i = 1:size(conf1,2)
%     plot(conf1(:,i))
%     while waitforbuttonpress ~= 1
%     end
% end
%%

fps = 60;
frames = length(inputFiles)-2;
time = 0:1/fps:(frames-1)/fps;

% Remove leading frames in which no pose was detected yet

% Find the index of the first non-empty element
firstNonEmptyIdx1 = find(~cellfun('isempty', pos(:,1)), 1);
firstNonEmptyIdx2 = find(~cellfun('isempty', pos(:,2)), 1);
firstNonEmptyIdx = max([firstNonEmptyIdx1 firstNonEmptyIdx2]);

% Remove the leading empty elements (if any)
if ~isempty(firstNonEmptyIdx)
    pos = pos(firstNonEmptyIdx:end,:);
end

playerData = readtable("D:\GP\test-results.xlsx",'Sheet','Test Results','ReadVariableNames',true,'VariableNamingRule','preserve');
rowIdx = playerData.PID == subjectNum;

%
if strcmp(trialName,'arrowhead') == 1 && strcmp(trialDir,'R') == 1 %trialName == 'arrowhead'
    endTime = playerData.ArrowR(rowIdx);
elseif strcmp(trialName,'arrowhead') == 1 && strcmp(trialDir,'L') == 1
    endTime = playerData.ArrowL(rowIdx);
elseif strcmp(trialName,'agilityFCT') == 1 && strcmp(trialDir,'R') == 1
    endTime = playerData.COD_ArrowR(rowIdx);
elseif strcmp(trialName,'agilityFCT') == 1 && strcmp(trialDir,'L') == 1
    endTime = playerData.COD_ArrowL(rowIdx);
end
%%
% size(pos)
% pos(1:,:)
% posNew = min(time(time>endTime));
% time(time == posNew)
% pos = pos(1:find(time == min(time(time>endTime))),:);
% % find(time == min(time(time>endTime)))
% time = time(1:size(pos,1));
% find(time == min(time(time>endTime)))
%%
clearvars k i inputFiles jsonName path str data fps frames vFile %trialAbbrev


% Load stereo parameters
% load(strcat('stereoCalibrationSessionCam',num2str(camNums(1)),num2str(camNums(2)),'.mat'));
% stereoCalibrationSession = calibrationSession;

% clearvars calibrationSession

% completionTimes = readtable("D:\GP\Trial-Data.xlsx",'ReadVariableNames',true,'VariableNamingRule','preserve');
% rowIdx = completionTimes.trial == trialName;

% testResults.completionTimeAFCT505 = completionTimes.agilityFCT_505(rowIdx);
% testResults.completionTimeAFCTTriangle = completionTimes.agilityFCT_triangle(rowIdx);
% testResults.completionTimeAFCTTotal = completionTimes.agilityFCT_total(rowIdx);
% testResults.completionTime1000m = completionTimes.("1000m")(rowIdx);
% testResults.completionTime400m = completionTimes.("400m")(rowIdx);
% testResults.completionTimeSprint = completionTimes.("sprint")(rowIdx);
% testResults.scoreJump = completionTimes.("jump")(rowIdx);
% testResults.scorePower = completionTimes.("power")(rowIdx);
% testResults.scoreMobility = completionTimes.("mobility")(rowIdx);

% clearvars completionTimes rowidx
%


% Anthropometric and player data
anthropoData.weight = playerData.Gewicht(rowIdx);
anthropoData.height = playerData.Lengte(rowIdx);
anthropoData.sittingHeight = playerData.Zithoogte(rowIdx);
anthropoData.age = playerData.Leeftijd(rowIdx);

% Sprint test results
testResults.sprintTest10m1st = playerData.("10m1")(rowIdx);
testResults.sprintTest10m2nd = playerData.("10m2")(rowIdx);
testResults.sprintTest10mAvg = mean([testResults.sprintTest10m1st testResults.sprintTest10m2nd]);
testResults.sprintTest30m1st = playerData.("30m1")(rowIdx);
testResults.sprintTest30m2nd = playerData.("30m2")(rowIdx);
testResults.sprintTest30mAvg = mean([testResults.sprintTest30m1st testResults.sprintTest30m2nd]);

% Agility FCT test results
testResults.AFCT505Left = playerData.COD_505L(rowIdx);
testResults.AFCT505Right = playerData.COD_505R(rowIdx);
testResults.AFCT505Avg = mean([testResults.AFCT505Left testResults.AFCT505Right]);
testResults.AFCTFullLeft = playerData.COD_ArrowL(rowIdx);
testResults.AFCTFullRight = playerData.COD_ArrowR(rowIdx);
testResults.AFCTFullAvg = mean([testResults.AFCTFullLeft testResults.AFCTFullRight]);

% 505 test results
testResults.Left505 = playerData.("505L")(rowIdx);
testResults.Right505 = playerData.("505R")(rowIdx);
testResults.Avg505 = mean([testResults.Left505 testResults.Right505]);

% Arrowhead test results
testResults.arrowheadLeft = playerData.ArrowL(rowIdx);
testResults.arrowheadRight = playerData.ArrowR(rowIdx);
testResults.arrowheadAvg = mean([testResults.arrowheadLeft testResults.arrowheadRight]);

% Endurance test results
testResults.endurance1000m = playerData.("1000_sec")(rowIdx);
testResults.endurance400m = playerData.("400_sec")(rowIdx);
%%
% clearvars rowIdx playerData

% Butterworth filter design:
% Define filter parameters
fs = 60;           % Sampling rate (Hz)
fc = 4;            % Cutoff frequency (Hz)
order = 2;         % Filter order

% Normalize the cutoff frequency
Wn = fc / (fs/2);

% Design the filter (low pass, to filter out high frequency outliers)
[b,a] = butter(order, Wn, 'low');

%Define variable for frame number plotting
frame = (1:1:length(pos));

% For all camera's, extract joints per frame from the pos variable, filter using 1-D median (medfilt1) and
% Butterworth filtering, and insert the filtered data back into the pos variable
for k = 1:length(camNums) 
    for i = 1:75
        for j = 1:length(pos)
            try
                joint(j,i) = pos{j,k}(i);
            catch
                continue
            end
        end
    end

    for i = 1:75
        joint_medf(:,i) = medfilt1(joint(:,i),5,'omitnan');
        joint_smooth(:,i) = filtfilt(b, a, joint_medf(:,i)')';
        % joint_smooth(:,i) = movmedian(joint(:,i),3);
        % joint_smooth(:,i) = filtfilt(b, a, joint_smooth(:,i)')';
        % joint_spline(:,i) = filloutliers(joint_smooth(:,i),'linear');
    end

    for i = 1:75
        for j = 1:length(pos)
            % pos_{j,k}(i) = joint_smooth(j,i);
            if ~isempty(pos(j,k))
                pos{j,k}(i,:) = joint_smooth(j,i);
            end
        end
    end
end

%
clearvars fs fc order Wn k i joint joint_medf joint_smooth b a

% Extract X and Y data from pos variable and load into separate variables
for i = 1:size(pos,1)
    posX1(i,:) = pos{i,1}(1:3:end);
    posY1(i,:) = pos{i,1}(2:3:end);
    posX2(i,:) = pos{i,2}(1:3:end);
    posY2(i,:) = pos{i,2}(2:3:end);
    % conf1(i,:) = pos{i,1}(3:3:end);
    % 
    % posX1_(i,:) = pos_{i,1}(1:3:end);
    % posY1_(i,:) = pos_{i,1}(2:3:end);
    % posX2_(i,:) = pos_{i,2}(1:3:end);
    % posY2_(i,:) = pos_{i,2}(2:3:end);
end

% Remove face joints as they are not needed (nrs 0-4)
del = 5; %5 joints + 1 for translation to Matlab numbering
posX1(:,1:del) = [];
posY1(:,1:del) = [];
posX2(:,1:del) = [];
posY2(:,1:del) = [];
% conf1(:,1:del) = [];

%%
% 
% figure;
% plot(conf1(:,10))

% posX1_(:,1:del) = [];
% posY1_(:,1:del) = [];
% posX2_(:,1:del) = [];
% posY2_(:,1:del) = [];

% for i = 1:size(pos,1)
%     pos1undistorted = undistortPoints([posX1(i,:)' posY1(i,:)'],initCals(1).calibrationSession.CameraParameters);
%     pos2undistorted = undistortPoints([posX2(i,:)' posY2(i,:)'],initCals(2).calibrationSession.CameraParameters);
%     posX1(i,:) = pos1undistorted(:,1)';
%     posY1(i,:) = pos1undistorted(:,2)';
%     posX2(i,:) = pos2undistorted(:,1)';
%     posY2(i,:) = pos2undistorted(:,2)';
% end

%%

%Triangulate 3D positional data per joint per frame based on stereo
%parameters. Load into separate X,Y,Z world coordinate variables
Rx = rotx(-90);
for i = 1:size(posX1,2)
    wPts = triangulate([posX1(:,i) posY1(:,i)], [posX2(:,i) posY2(:,i)],stereoCalibrationSession.CameraParameters);
    aPts = [wPts(:,1) wPts(:,2) wPts(:,3)];
    for k = 1:length(aPts)
        RPts(:,k) = Rx*aPts(k,:)';
    end
    WX(:,i) = RPts(1,:)';
    WY(:,i) = RPts(2,:)';
    WZ(:,i) = RPts(3,:)';
end

%%
% -1/(sqrt(2)*erfcinv(3/2))
%%
% clearvars posX1 posX2 posY1 posY2 Rx i wPts aPts RPts

SIP_table = readtable("C:\Users\jaspe\OneDrive - University of Twente\Documenten\University of Twente\" + ...
    "Year 2\Graduation Project\OpenPose\Segment Inertia Parameters - De Leva 1996.xlsx");

SIP_thigh = strcmp(SIP_table.segment, 'thigh');
real_thigh = ((playerData.Lengte(rowIdx)*10) * table2array(SIP_table(SIP_thigh, 3)))/1741;

for i = 1:size(pos,1) %for each frame
    virtual_thigh_len(i,:) = norm([WX(i,10) WY(i,10) WZ(i,10)] - [WX(i,8) WY(i,8) WZ(i,8)]);
end

virt_thigh_avg = mean(rmoutliers(virtual_thigh_len));
scale_factor = real_thigh/virt_thigh_avg;
% scale_factor = 1.5;

for j = 1:size(WX,2) %for each openpose keypoint (body_25B with head joints removed except top of head)
    WX(:,j) = (WX(:,j)) * scale_factor;
    WY(:,j) = (WY(:,j)) * scale_factor;
    WZ(:,j) = (WZ(:,j)) * scale_factor;
end

clearvars SIP_thigh 
% clearvars virt_thigh_avg virtual_thigh_len scale_factor real_thigh

% Calculate joint angles
for i = 1:size(pos,1)
    %right knee
    knee_vec_r = [WX(i,14-del+1) WY(i,14-del+1) WZ(i,14-del+1)];
    init_vec1 = [WX(i,12-del+1) WY(i,12-del+1) WZ(i,12-del+1)];    
    init_vec2 = [WX(i,16-del+1) WY(i,16-del+1) WZ(i,16-del+1)];
    vec1 = knee_vec_r - init_vec1;
    vec2 = knee_vec_r - init_vec2;
    % vec3 = init_vec2 - knee_vec_r;
    jointAngles.kneeAngle(i,1) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));
    % jointAngles.kneeAngle(i,3) = atan2d(norm(cross(vec1,vec3)),dot(vec1,vec3));

    %left knee
    knee_vec_l = [WX(i,13-del+1) WY(i,13-del+1) WZ(i,13-del+1)];
    init_vec1 = [WX(i,11-del+1) WY(i,11-del+1) WZ(i,11-del+1)];
    init_vec2 = [WX(i,15-del+1) WY(i,15-del+1) WZ(i,15-del+1)];
    vec1 = knee_vec_l - init_vec1;
    vec2 = knee_vec_l - init_vec2;
    % vec3 = init_vec2 - knee_vec_l;
    jointAngles.kneeAngle(i,2) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));
    % jointAngles.kneeAngle(i,4) = atan2d(norm(cross(vec1,vec3)),dot(vec1,vec3));

    %right hip
    hip_vec_r = [WX(i,12-del+1) WY(i,12-del+1) WZ(i,12-del+1)];
    init_vec1 = [WX(i,6-del+1) WY(i,6-del+1) WZ(i,6-del+1)];    
    init_vec2 = [WX(i,14-del+1) WY(i,14-del+1) WZ(i,14-del+1)];
    vec1 = hip_vec_r - init_vec1;
    vec2 = hip_vec_r - init_vec2;
    jointAngles.hipAngle(i,1) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));

    %left hip
    hip_vec_l = [WX(i,11-del+1) WY(i,11-del+1) WZ(i,11-del+1)];
    init_vec1 = [WX(i,5-del+1) WY(i,5-del+1) WZ(i,5-del+1)];    
    init_vec2 = [WX(i,13-del+1) WY(i,13-del+1) WZ(i,13-del+1)];
    vec1 = hip_vec_l - init_vec1;
    vec2 = hip_vec_l - init_vec2;
    jointAngles.hipAngle(i,2) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));

    %right ankle
    ankle_vec_r = [WX(i,16-del+1) WY(i,16-del+1) WZ(i,16-del+1)];
    init_vec1 = [WX(i,14-del+1) WY(i,14-del+1) WZ(i,14-del+1)];    
    init_vec2 = [WX(i,23-del+1) WY(i,23-del+1) WZ(i,23-del+1)];
    vec1 = ankle_vec_r - init_vec1;
    vec2 = ankle_vec_r - init_vec2;
    jointAngles.ankleAngle(i,1) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));

    %left ankle
    ankle_vec_l = [WX(i,15-del+1) WY(i,15-del+1) WZ(i,15-del+1)];
    init_vec1 = [WX(i,13-del+1) WY(i,13-del+1) WZ(i,13-del+1)];    
    init_vec2 = [WX(i,19-del+1) WY(i,19-del+1) WZ(i,19-del+1)];
    vec1 = ankle_vec_l - init_vec1;
    vec2 = ankle_vec_l - init_vec2;
    jointAngles.ankleAngle(i,2) = atan2d(norm(cross(vec1,vec2)),dot(vec1,vec2));
end

clearvars knee_vec_l knee_vec_r hip_vec_l hip_vec_r ankle_vec_l ankle_vec_r init_vec1 init_vec2 vec1 vec2 i

% jointAngles.kneeAngle = rmoutliers(jointAngles.kneeAngle);
% jointAngles.hipAngle = rmoutliers(jointAngles.hipAngle);
% jointAngles.ankleAngle = rmoutliers(jointAngles.ankleAngle);

% Filter both knee joint angles
% for k = 1:2
%     jointAngles.kneeAngle(:,k) = medfilt1(jointAngles.kneeAngle(:,k),2,'omitnan');
%     jointAngles.hipAngle(:,k) = medfilt1(jointAngles.hipAngle(:,k),2,'omitnan');
%     jointAngles.ankleAngle(:,k) = medfilt1(jointAngles.ankleAngle(:,k),2,'omitnan');
% end

% Define the edges of the angle bins and bin the angles
edges = 0:10:180;
jointAngles.kneeAngleCounter(:,1) = histcounts(jointAngles.kneeAngle(:,1), edges)'; %right knee
jointAngles.kneeAngleCounter(:,2) = histcounts(jointAngles.kneeAngle(:,2), edges)'; %left knee
jointAngles.hipAngleCounter(:,1) = histcounts(jointAngles.hipAngle(:,1), edges)'; %right hip
jointAngles.hipAngleCounter(:,2) = histcounts(jointAngles.hipAngle(:,2), edges)'; %left hip
jointAngles.ankleAngleCounter(:,1) = histcounts(jointAngles.ankleAngle(:,1), edges)'; %right ankle
jointAngles.ankleAngleCounter(:,2) = histcounts(jointAngles.ankleAngle(:,2), edges)'; %left ankle

%indicate start and finish frames
% start = 150;
% finish = 550;

% Center of mass
% The CoM position is calculated using the segmentation method, where the CoM is estimated 
% by summing the moments of masses of individual body segments (= weighted segment CoM positions).
% CoM = sum of all segments( segment mass percentage of body * coordinates of individual segments )

for i = 1:size(pos,1)
    
    cm = [];
    %start + l(end-start)
    %calculate center of mass of left upper arm (male)
    segment = strcmp(SIP_table.segment, 'upper arm');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,1) WY(i,1) WZ(i,1)] + lcmp * ([WX(i,3) WY(i,3) WZ(i,3)] - [WX(i,1) WY(i,1) WZ(i,1)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of right upper arm (male)
    segment_cm = smass * ([WX(i,2) WY(i,2) WZ(i,2)] + lcmp * ([WX(i,4) WY(i,4) WZ(i,4)] - [WX(i,2) WY(i,2) WZ(i,2)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of left forearm (male)
    segment = strcmp(SIP_table.segment, 'forearm');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,3) WY(i,3) WZ(i,3)] + lcmp * ([WX(i,5) WY(i,5) WZ(i,5)] - [WX(i,3) WY(i,3) WZ(i,3)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of right forearm (male)
    segment_cm = smass * ([WX(i,4) WY(i,4) WZ(i,4)] + lcmp * ([WX(i,6) WY(i,6) WZ(i,6)] - [WX(i,4) WY(i,4) WZ(i,4)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of left thigh (male)
    segment = strcmp(SIP_table.segment, 'thigh');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,7) WY(i,7) WZ(i,7)] + lcmp * ([WX(i,9) WY(i,9) WZ(i,9)] - [WX(i,7) WY(i,7) WZ(i,7)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of right thigh (male)
    segment_cm = smass * ([WX(i,8) WY(i,8) WZ(i,8)] + lcmp * ([WX(i,10) WY(i,10) WZ(i,10)] - [WX(i,8) WY(i,8) WZ(i,8)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of left shank (male)
    segment = strcmp(SIP_table.segment, 'shank');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,9) WY(i,9) WZ(i,9)] + lcmp * ([WX(i,11) WY(i,11) WZ(i,11)] - [WX(i,9) WY(i,9) WZ(i,9)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of right shank (male)
    segment_cm = smass * ([WX(i,10) WY(i,10) WZ(i,10)] + lcmp * ([WX(i,12) WY(i,12) WZ(i,12)] - [WX(i,10) WY(i,10) WZ(i,10)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table
    
    %calculate center of mass of head (male)
    segment = strcmp(SIP_table.segment, 'head');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,14) WY(i,14) WZ(i,14)] + lcmp * ([WX(i,13) WY(i,13) WZ(i,13)] - [WX(i,14) WY(i,14) WZ(i,14)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

     %calculate center of mass of trunk (male)
    segment = strcmp(SIP_table.segment, 'trunk');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    MIDH = mean([[WX(i,7) WY(i,7) WZ(i,7)];[WX(i,8) WY(i,8) WZ(i,8)]]); %middle point between hips (Mid Hip)
    segment_cm = smass * ([WX(i,13) WY(i,13) WZ(i,13)] + lcmp * (MIDH - [WX(i,13) WY(i,13) WZ(i,13)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table
    
    %calculate center of mass of left foot (male)
    segment = strcmp(SIP_table.segment, 'foot');
    lcmp = table2array(SIP_table(segment, 7))/100; %longitudinal center of mass position (%)
    smass = table2array(SIP_table(segment, 5))/100; %segment mass, percentage of full body weight
    segment_cm = smass * ([WX(i,17) WY(i,17) WZ(i,17)] + lcmp * ([WX(i,15) WY(i,15) WZ(i,15)] - [WX(i,17) WY(i,17) WZ(i,17)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    %calculate center of mass of right foot (male)
    segment_cm = smass * ([WX(i,20) WY(i,20) WZ(i,20)] + lcmp * ([WX(i,19) WY(i,19) WZ(i,19)] - [WX(i,20) WY(i,20) WZ(i,20)]));
    cm = vertcat(cm, segment_cm); %add segment center of mass to table

    % com_forearm = [WX(frame,3) WY(frame,3) WZ(frame,3)] + lcpm_forearm * ([WX(frame,5) WY(frame,5) WZ(frame,5)] - [WX(frame,3) WY(frame,3) WZ(frame,3)]);
    %benadering waarbij m_i bij beide 50% is, in realiteit beide segments vermenigvuldigen met waarde uit tabel!
    % com_arm = (com_upperarm + com_forearm)/2; 
    %com_arm = (0.5 * com_upperarm) + (0.5 * com_forearm)
    % segment = strcmp(SIP_table.segment, 'foot');
    % footmassperc = table2array(SIP_table(strcmp(SIP_table.segment, 'foot'), 5))/100;
    handmassperc = table2array(SIP_table(strcmp(SIP_table.segment, 'hand'), 5))/100;
    centerofmass(i,:) = sum(cm)/(1-handmassperc);
end

clearvars k i segment lcmp smass segment_cm cm handmassperc

%%
% close all
clc 

%designing Butterworth filter for kinematic data of CoM, using 60 fps as
%sampling rate and a cutoff frequency of 1 Hz based on an iterative
%approach to reach the desired effect
fs = 60;   % Sampling rate (Hz)
fc = 1;            % Cutoff frequency (Hz)
order = 2;         % Filter order
Wn = fc / (fs/2);
[b,a] = butter(order, Wn, 'low');
centerofmass = filtfilt(b,a,centerofmass); 

gap = 1; % equals to measuring distance every 1/12th of a second (when using 60fps)
for i = 1+gap:size(pos,1)
    currentComLoc = centerofmass(i,:);
    previousComLoc = centerofmass(i-gap,:);
    displacement(i) = norm(currentComLoc - previousComLoc);
    % direction = displacement / norm(displacement);
    % dirDistCoM = dot(displacement,direction);
    % speedCoM(i-gap,:) = displacement / time(gap+1); % in mm per second!
end

%calculate displacement of center of mass between each frame
dcom = diff(centerofmass);
freq = 60;
velocityCoM = (sqrt(dcom(:,1).^2 + dcom(:,2).^2 + dcom(:,3).^2) * freq)/1000;
accelerationCoM = (gradient(velocityCoM, 1/freq));


% ax1 = subplot(2,1,1);
% hold on
% xlabel('time (s)')
% ylabel('Speed (m/2)')
% plot(time(idx_start+1:idx_end),velocityCoM,'b-')
% 
% ax1 = subplot(2,1,2);
% hold on
% xlabel('time (s)')
% ylabel('Acceleration (m/s^2)')
% plot(time(idx_start+1:idx_end),accelerationCoM,'r-')

if trialDir == 'L' && strcmp(trialName,'arrowhead') == 1
    completionTime = testResults.arrowheadLeft;
elseif trialDir == 'L' && strcmp(trialName,'agilityFCT') == 1
    completionTime = testResults.AFCTFullLeft;
elseif trialDir == 'R' && strcmp(trialName,'arrowhead') == 1
    completionTime = testResults.arrowheadRight;
elseif trialDir == 'R' && strcmp(trialName,'agilityFCT') == 1
    completionTime = testResults.AFCTFullRight;
end

if completionTime >= 9
    lineColor = 'r';
elseif completionTime <= 8.5
    lineColor = 'g';
else
    lineColor = 'c';
end

removeFrames = 90;
limit = 20000;
yCOM = centerofmass(:,2);
yCOM(1:removeFrames) = limit; yCOM(end-removeFrames:end) = limit;

xCOM = centerofmass(:,1);
xCOM(1:removeFrames) = limit; xCOM(end-removeFrames:end) = limit;

if strcmp(trialDir,'L') == 1
    [pkX,locX] = (findpeaks(xCOM));
elseif strcmp(trialDir,'R') == 1
    [pkX,locX] = (findpeaks(-xCOM));
end

maxIdX = find(pkX == max(pkX));
threshIdX = locX(maxIdX);
thresholdX = yCOM(threshIdX);

newdata = yCOM(yCOM >= thresholdX);
newdata_idxs = find(yCOM >= thresholdX);
[pks,locs] = findpeaks(diff(newdata_idxs));

while size(locs,1) > 1
    % i = i - 10;
    thresholdX = thresholdX - 10;
    newdata = yCOM(yCOM >= thresholdX);
    newdata_idxs = find(yCOM >= thresholdX);
    [pks,locs] = findpeaks(diff(newdata_idxs));
    % disp(threshold)
end

if strcmp(trialDir,'L') == 1
    threshPart = min(xCOM(xCOM(1:locs) > max(pkX)/2));
elseif strcmp(trialDir,'R') == 1
    threshPart = max(xCOM(xCOM(1:locs) < -max(pkX)/2));
end

idx_start = find(xCOM == max(threshPart));
idx_end = locs+pks;

% figure;
% hold on
% axis equal
% % plot(xCOM)
% plot(centerofmass(:,1),centerofmass(:,2))
% xline(max(threshPart))
% plot(centerofmass(idx_start,1),centerofmass(idx_start,2),'r.')

% ax1 = subplot(2,1,1);
% hold on
% % plot(time(1:size(velocityCoM)),velocityCoM,'Color','#1C75BC','LineStyle','-','LineWidth',1)
% plot(time(idx_start:idx_end),velocityCoM(idx_start:idx_end),'Color','#1C75BC','LineStyle','-','LineWidth',1)
% % xline([time(idx_start) time(idx_end)])
% 
% ax2 = subplot(2,1,2);
% hold on
% % plot(time(1:size(accelerationCoM)),accelerationCoM,'Color','#EF4136','LineStyle','-','LineWidth',1)
% plot(time(idx_start:idx_end),accelerationCoM(idx_start:idx_end),'Color','#EF4136','LineStyle','-','LineWidth',1)
% % xline([time(idx_start) time(idx_end)])

WX = WX(idx_start:idx_end,:);
WY = WY(idx_start:idx_end,:);
WZ = WZ(idx_start:idx_end,:);
pos = pos(idx_start:idx_end,:);
centerofmass = centerofmass(idx_start:idx_end,:);
posX1 = posX1(idx_start:idx_end,:);
posY1 = posY1(idx_start:idx_end,:);
posX2 = posX2(idx_start:idx_end,:);
posY2 = posY2(idx_start:idx_end,:);

% Calculating speed and acceleration of joints with the purpose of foot ground contact time

heelLeft = 17;
toeLeft = 15;
heelRight = 20;
toeRight = 19;

gap = 1;
for i = 1+gap:size(pos,1)
    %left new
    currentHeelLeft = WZ(i,heelLeft); %vertical
    previousHeelLeft = WZ(i-gap,heelLeft);
    displacementL = norm(currentHeelLeft - previousHeelLeft);
    vertSpeedHeelL(i-gap,:) = (displacementL / time(gap+1))/1000; % in mm per second!

    currentToeLeft = WZ(i,toeLeft); %vertical
    previousToeLeft = WZ(i-gap,toeLeft);
    displacementL = norm(currentToeLeft - previousToeLeft);
    vertSpeedToeL(i-gap,:) = (displacementL / time(gap+1))/1000; % in mm per second!

    currentHeelRight = WZ(i,heelRight); %vertical
    previousHeelRight = WZ(i-gap,heelRight);
    displacementR = norm(currentHeelRight - previousHeelRight);
    vertSpeedHeelR(i-gap,:) = (displacementR / time(gap+1))/1000; % in mm per second!

    currentToeRight = WZ(i,toeRight); %vertical
    previousToeRight = WZ(i-gap,toeRight);
    displacementR = norm(currentToeRight - previousToeRight);
    vertSpeedToeR(i-gap,:) = (displacementR / time(gap+1))/1000; % in mm per second!
end

[pks_lhs,locs_lhs] = findpeaks(-vertSpeedHeelL); 
[pks_lto,locs_lto] = findpeaks(vertSpeedToeL);
[pks_rhs,locs_rhs] = findpeaks(-vertSpeedHeelR); 
[pks_rto,locs_rto] = findpeaks(vertSpeedToeR);

events_openpose.locs_lhs_filt = []; events_openpose.pks_lhs_filt = [];
for j = 2:length(pks_lhs)-1
    if WZ(locs_lhs(j),heelLeft) < WZ(locs_lhs(j-1),heelLeft) && WZ(locs_lhs(j),heelLeft) < WZ(locs_lhs(j+1),heelLeft)
        events_openpose.locs_lhs_filt = [events_openpose.locs_lhs_filt; locs_lhs(j)];
        events_openpose.pks_lhs_filt = [events_openpose.pks_lhs_filt; pks_lhs(j)];
    end
end

events_openpose.locs_rhs_filt = []; events_openpose.pks_rhs_filt = [];
for j = 2:length(pks_rhs)-1
    if WZ(locs_rhs(j),heelRight) < WZ(locs_rhs(j-1),heelRight) && WZ(locs_rhs(j),heelRight) < WZ(locs_rhs(j+1),heelRight)
        events_openpose.locs_rhs_filt = [events_openpose.locs_rhs_filt; locs_rhs(j)];
        events_openpose.pks_rhs_filt = [events_openpose.pks_rhs_filt; pks_rhs(j)];
    end
end

events_openpose.locs_lto_filt = []; events_openpose.pks_lto_filt = [];
for j = 1:length(events_openpose.locs_lhs_filt)
    %find toe-off event between current heel strike and subsequent knee
    %raise moment
    currentLhs = events_openpose.locs_lhs_filt(j);
    nextLeftKA = min(locs_lhs(locs_lhs > events_openpose.locs_lhs_filt(j))); %knee adduction
    toeoffs = locs_lto(locs_lto > currentLhs & locs_lto < nextLeftKA);

    %if more than 1 toe-off is found, take the highest peak in vertical
    %speed.
    if length(toeoffs) > 1
        highestTO = 0;
        for i = 1:length(toeoffs)
            temp_idx = find(locs_lto == toeoffs(i));
            if pks_lto(temp_idx) > highestTO
                highestTO = pks_lto(temp_idx);
                temp_loc = locs_lto(temp_idx);
                temp_pks = pks_lto(temp_idx);
            end
        end
        events_openpose.locs_lto_filt = [events_openpose.locs_lto_filt; temp_loc];
        events_openpose.pks_lto_filt = [events_openpose.pks_lto_filt; temp_pks];
    elseif length(toeoffs) == 1
        temp_loc = min(locs_lto(locs_lto > events_openpose.locs_lhs_filt(j)));
        temp_idx = find(locs_lto == temp_loc);
        temp_pks = pks_lto(temp_idx);
        events_openpose.locs_lto_filt = [events_openpose.locs_lto_filt; temp_loc];
        events_openpose.pks_lto_filt = [events_openpose.pks_lto_filt; temp_pks];
    else
        if ~isempty(locs_lto(locs_lto > events_openpose.locs_lhs_filt(j)))
            temp_loc = min(locs_lto(locs_lto > events_openpose.locs_lhs_filt(j)));
            temp_idx = find(locs_lto == temp_loc);
            temp_pks = pks_lto(temp_idx);
            events_openpose.locs_lto_filt = [events_openpose.locs_lto_filt; temp_loc];
            events_openpose.pks_lto_filt = [events_openpose.pks_lto_filt; temp_pks];
        else
            continue
        end
    end
end

events_openpose.locs_rto_filt = []; events_openpose.pks_rto_filt = [];
for j = 1:length(events_openpose.locs_rhs_filt)
    currentRhs = events_openpose.locs_rhs_filt(j);
    nextRightKA = min(locs_rhs(locs_rhs > events_openpose.locs_rhs_filt(j))); %knee adduction
    toeoffs = locs_rto(locs_rto > currentRhs & locs_rto < nextRightKA);

    if length(toeoffs) > 1
        highestTO = 0;
        for i = 1:length(toeoffs)
            temp_idx = find(locs_rto == toeoffs(i));
            if pks_rto(temp_idx) > highestTO
                highestTO = pks_rto(temp_idx);
                temp_loc = locs_rto(temp_idx);
                temp_pks = pks_rto(temp_idx);
            end
        end
        events_openpose.locs_rto_filt = [events_openpose.locs_rto_filt; temp_loc];
        events_openpose.pks_rto_filt = [events_openpose.pks_rto_filt; temp_pks];
    elseif length(toeoffs) == 1
        temp_loc = min(locs_rto(locs_rto > events_openpose.locs_rhs_filt(j)));
        temp_idx = find(locs_rto == temp_loc);
        temp_pks = pks_rto(temp_idx);
        events_openpose.locs_rto_filt = [events_openpose.locs_rto_filt; temp_loc];
        events_openpose.pks_rto_filt = [events_openpose.pks_rto_filt; temp_pks];
    else
        if ~isempty(locs_rto(locs_rto > events_openpose.locs_rhs_filt(j)))
            temp_loc = min(locs_rto(locs_rto > events_openpose.locs_rhs_filt(j)));
            temp_idx = find(locs_rto == temp_loc);
            temp_pks = pks_rto(temp_idx);
            events_openpose.locs_rto_filt = [events_openpose.locs_rto_filt; temp_loc];
            events_openpose.pks_rto_filt = [events_openpose.pks_rto_filt; temp_pks];
        else
            continue
        end
        
    end
end

events_openpose.lhs_frames = events_openpose.locs_lhs_filt'; events_openpose.lto_frames = events_openpose.locs_lto_filt';
events_openpose.rhs_frames = events_openpose.locs_rhs_filt'; events_openpose.rto_frames = events_openpose.locs_rto_filt';

clearvars locs_rto locs_rhs locs_lto locs_lhs pks_rto pks_rhs pks_lto pks_lhs temp_idx temp_loc temp_pks ...
highestTO currentRhs currentLhs nextRightKA nextLeftKA
 

%%
% create bounds to ignore the first and last 1/8 of total frames, which are too far away
% from the camera for accurate pose detection
% i = 10;
% bound_lower = round(size(centerofmass,1)/i);
% bound_upper = round(size(centerofmass,1)/i*(i-1));
% acc_withinBounds = accelerationCoM(bound_lower:bound_upper);
% 
% while max(acc_withinBounds) > 9.8 || min(acc_withinBounds) < -9.8
%     i = i - 1;
%     bound_lower = round(size(centerofmass,1)/i);
%     bound_upper = round(size(centerofmass,1)/i*(i-1));
%     acc_withinBounds = accelerationCoM(bound_lower:bound_upper);
% end

% gap = 1; % equals to measuring distance every 1/6th of a second (when using 60fps)
% for i = 1+gap:size(pos,1)
%     currentComLoc = centerofmass(i,:);
%     previousComLoc = centerofmass(i-gap,:);
%     displacement = norm(currentComLoc - previousComLoc); % in mm
% 
%     currentTime = time(i);
%     previousTime = time(i-gap);
%     duration = currentTime - previousTime; % in seconds
%     % direction = displacement / norm(displacement);
%     % dirDistCoM = dot(displacement,direction);
%     speedCoM(i,:) = displacement / duration; % in mm per second!
% end
% close all
% clc 

% f1 = figure("WindowState","maximized");
% 
% tiledlayout(2,1)
% nexttile
% hold on;
% ylim([0 20000])
% plot(speedCoM,'b')
% plot(rmoutliers(speedCoM),'r')
% title("old method")

% nexttile
% hold on;
% ylim([0 20000])
% plot(speedCoM_new,'b')
% plot(rmoutliers(speedCoM_new),'r')
% % plot(filloutliers(displacement,"spline"),'g')
% ylabel("Distance")
% xlabel("Time (in frames)")
% title("new method")

% f2 = figure;
% hold on
% grid on
% ylabel('Y')
% xlabel('X')
% zlabel('Z')
% xlim([-2000 9000]);
% ylim([0 25000]);
% zlim([-2000 3000]);
% view([90 0]); %plot ZY plane (side view)
% view([90 90]); %plot XY plane (topview)
% title('Running pattern')
% 
% for i = 1:size(centerofmass,1)
%     plot3(centerofmass(i,1),centerofmass(i,2),centerofmass(i,3),'Color',[i/size(pos,1) (size(pos,1)-i)/size(pos,1) 1],'Marker','.');
% end
% plot(1500,4500,'ro')
% plot(1500,9500,'ro')
% plot(1500,14500,'ro')
% plot(1500,23500,'ro')
% plot(6500,9500,'ro')
% plot(0,9500,'ro')
%%
% curve = fit(time(1:size(rmoutliers(displacement),2))',rmoutliers(displacement)','smoothingspline');
% diffCurve = differentiate(curve,time(1:size(rmoutliers(displacement),2)));
% 
% figure;
% plot(curve,time(1:size(rmoutliers(displacement),2)),rmoutliers(displacement));
% 
% figure;
% ylim([0 1000])
% plot(time(1:size(rmoutliers(displacement),2)),diffCurve);

% size(displacement,2)
% time(1:size)

% diffCurve = differentiate(curve,time(start:fin));
%%
% speedCoM_ms = speedCoM/1000; %speed in m/s
% speedCoM_ms_filtered = filloutliers(speedCoM,"center")/1000; %filtered speed in m/s
% speedCoM_ms_filt_smooth = smoothdata(speedCoM_ms_filtered,'gaussian',60);
% acc_ms2 = (diff(speedCoM_ms) / time(2));
% acc_ms2_filtered = (diff(speedCoM_ms_filtered) / time(2));
% acc_ms2_filt_smooth = (diff(speedCoM_ms_filt_smooth)/time(2));
% acc_ms2_filt_smooth2 = smoothdata(acc_ms2_filt_smooth,'gaussian',10);

% speedCoM_ms = speedCoM/1000; %speed in m/s
% speedCoM_filt = smoothdata(filloutliers(speedCoM,"center")/1000,'gaussian',60);
% accCoM_filt = smoothdata((diff(speedCoM_filt)/time(2)),'gaussian',10);

% fs = size(pos,1);   % Sampling rate (Hz)
% fc = 12;            % Cutoff frequency (Hz)
% order = 4;         % Filter order
% 
% % Normalize the cutoff frequency
% Wn = fc / (fs/2);
% 
% % Design the filter (low pass, to filter out high frequency outliers)
% [b,a] = butter(order, Wn, 'low');
% 
% speedCoM_filt = filtfilt(b, a, filloutliers(speedCoM,"center")/1000);
% accCoM = (diff(speedCoM_filt)/time(2));
% accCoM_filt = filtfilt(b,a,accCoM);

%%

% speedCoM_ = medfilt1(speedCoM,10,'omitnan');
% accelerationCoM = (diff(speedCoM,2) / time(gap+1));

clearvars i gap currentComLoc previousComLoc displacement direction dirDistCoM meanSpeedRaw meanAccelerationRaw

% save('data_processing')

%% Calculate Agility Parameters

% clear variables
% close all
% 
% load("data_processing.mat");

% calculate takeoff distance (i.e. distance from trail leg to center of mass at toe off)
for j = 1:length(events_openpose.lto_frames)
    temp_L(j) = norm(centerofmass(events_openpose.lto_frames(:,j),:) - ...
        [WX(events_openpose.lto_frames(:,j),16) WY(events_openpose.lto_frames(:,j),16) ...
        WZ(events_openpose.lto_frames(:,j),16)])/1000;
end
agilityParams.takeoffDistance.left = rmoutliers(temp_L); clearvars temp_L

for j = 1:length(events_openpose.rto_frames)
    temp_R(j) = norm(centerofmass(events_openpose.rto_frames(:,j),:) - ...
        [WX(events_openpose.rto_frames(:,j),18) WY(events_openpose.rto_frames(:,j),18) ...
        WZ(events_openpose.rto_frames(:,j),18)])/1000;
end
agilityParams.takeoffDistance.right = rmoutliers(temp_R); clearvars temp_R

% calculate stance times left and right
for j = 1:length(events_openpose.lhs_frames)
    if ~isempty(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.lhs_frames(j)))
        temp_L(j) = (time(min(events_openpose.lto_frames(events_openpose.lto_frames > events_openpose.lhs_frames(j))) - ...
            events_openpose.lhs_frames(j)));
    end
end
agilityParams.stanceTime.left = rmoutliers(temp_L); clearvars j temp_L 
agilityParams.stanceTime.left(agilityParams.stanceTime.left <= time(2)) = [];

for j = 1:length(events_openpose.rhs_frames)
    if ~isempty(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.rhs_frames(j)))
        temp_R(j) = time(min(events_openpose.rto_frames(events_openpose.rto_frames > events_openpose.rhs_frames(j))) - ...
            events_openpose.rhs_frames(j));
    end
end
agilityParams.stanceTime.right = rmoutliers(temp_R); clearvars j temp_R 
agilityParams.stanceTime.right(agilityParams.stanceTime.right <= time(2)) = [];
% agilityParams.stanceTime.right <= 0.0167
% agilityParams.stanceTime.right(agilityParams.stanceTime.right <= 0.0167) = []

% calculate step times
for j = 1:length(events_openpose.rhs_frames)
    if ~isempty(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.rhs_frames(j)))
        temp_L(j) = time(min(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.rhs_frames(j))) - ...
            events_openpose.rhs_frames(j));
    end
end
agilityParams.stepTime.left = rmoutliers(temp_L); clearvars j temp_L
agilityParams.stepTime.left(agilityParams.stepTime.left <= time(2)) = [];

for j = 1:length(events_openpose.lhs_frames)
    if ~isempty(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.lhs_frames(j)))
        temp_R(j) = time(min(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.lhs_frames(j))) - ...
            events_openpose.lhs_frames(j));
    end
end
agilityParams.stepTime.right = rmoutliers(temp_R); clearvars j temp_R
agilityParams.stepTime.right(agilityParams.stepTime.right <= time(2)) = [];

% calculate swing times
for j = 1:length(events_openpose.lhs_frames)
    if ~isempty(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.lhs_frames(j)))
        temp_L(j) = time(min(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.lhs_frames(j)))) - ...
            time(min(events_openpose.lto_frames(events_openpose.lto_frames > events_openpose.lhs_frames(j))));
    end
end
agilityParams.swingTime.left = rmoutliers(temp_L); clearvars j temp_L 
agilityParams.swingTime.left(agilityParams.swingTime.left <= time(2)) = [];
    
for j = 1:length(events_openpose.rhs_frames)
    if ~isempty(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.rhs_frames(j)))
        temp_R(j) = time(min(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.rhs_frames(j)))) - ...
            time(min(events_openpose.rto_frames(events_openpose.rto_frames > events_openpose.rhs_frames(j))));
    end
end
agilityParams.swingTime.right = rmoutliers(temp_R); clearvars j temp_R
agilityParams.swingTime.right(agilityParams.swingTime.right <= time(2)) = [];

%%
% remove last heel strike if there is no toe off to match with
% events_openpose.lhs_frames(2)
%%
while length(events_openpose.lhs_frames) ~= length(events_openpose.lto_frames)
    if events_openpose.lhs_frames(1) > events_openpose.lto_frames(1) 
        events_openpose.lto_frames(1) = [];
    else
        events_openpose.lhs_frames(end) = [];
    end
end

while length(events_openpose.rhs_frames) ~= length(events_openpose.rto_frames)
    if events_openpose.rhs_frames(1) > events_openpose.rto_frames(1) 
        events_openpose.rto_frames(1) = [];
    else
        events_openpose.rhs_frames(end) = [];
    end
end

% calculate double support times, in which negative values indicate that
% the heel strike of the next step came after the toe off of the current
% step, meaning during the resulting time no feet were touching the ground.

for j = 1:length(events_openpose.lhs_frames)
    if ~isempty(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.lhs_frames(j)))
        try
            temp_LR(j) = time(min(events_openpose.lto_frames(events_openpose.lto_frames > events_openpose.lhs_frames(j)))) - ...
                time(min(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.lhs_frames(j))));
        catch
            disp('agilityParams.dsTime.left_to_right failed')
        end
    end
end
agilityParams.dsTime.left_to_right = rmoutliers(temp_LR); clearvars j temp_LR   
   
for j = 1:length(events_openpose.rhs_frames)
    if ~isempty(events_openpose.rhs_frames(events_openpose.rhs_frames > events_openpose.rhs_frames(j)))
        try
            temp_RL(j) = time(min(events_openpose.rto_frames(events_openpose.rto_frames > events_openpose.rhs_frames(j)))) - ...
                time(min(events_openpose.lhs_frames(events_openpose.lhs_frames > events_openpose.rhs_frames(j))));
        catch
            disp('agilityParams.dsTime.right_to_left failed')
        end
    end
end
agilityParams.dsTime.right_to_left = rmoutliers(temp_RL); clearvars j temp_RL

% calculate step lengths
for j = 1:length(events_openpose.rhs_frames)
    temp_R(j) = norm([WX(events_openpose.rhs_frames(:,j),12) WY(events_openpose.rhs_frames(:,j),12) WZ(events_openpose.rhs_frames(:,j),12)] - ...
    [WX(events_openpose.rhs_frames(:,j),11) WY(events_openpose.rhs_frames(:,j),11) WZ(events_openpose.rhs_frames(:,j),11)])/1000;
end
agilityParams.stepLength.right = rmoutliers(temp_R); clearvars j temp_R %remove step lengths greater than 3 median absolute deviations

for j = 1:length(events_openpose.lhs_frames)
    temp_L(j) = norm([WX(events_openpose.lhs_frames(:,j),11) WY(events_openpose.lhs_frames(:,j),11) WZ(events_openpose.lhs_frames(:,j),11)] - ...
    [WX(events_openpose.lhs_frames(:,j),12) WY(events_openpose.lhs_frames(:,j),12) WZ(events_openpose.lhs_frames(:,j),12)])/1000;
end
agilityParams.stepLength.left = rmoutliers(temp_L); clearvars j temp_L %remove step lengths greater than 3 median absolute deviations

% Calculate agility predictors to feed into the model

% Takeoff Distance (i.e. distance between trail leg and center of mass at toe off)
agilityPredictors.takeoffDistance.meanLeftTakeoffDistance = mean(agilityParams.takeoffDistance.left);
agilityPredictors.takeoffDistance.meanRightTakeoffDistance = mean(agilityParams.takeoffDistance.right);
agilityPredictors.takeoffDistance.meanBilateralTakeoffDistance = mean([agilityParams.takeoffDistance.left agilityParams.takeoffDistance.right]);
agilityPredictors.takeoffDistance.maxLeftTakeoffDistance = max(agilityParams.takeoffDistance.left);
agilityPredictors.takeoffDistance.maxRightTakeoffDistance = max(agilityParams.takeoffDistance.right);
agilityPredictors.takeoffDistance.maxBilateralTakeoffDistance = max([agilityParams.takeoffDistance.left agilityParams.takeoffDistance.right]);


% Stance Time (i.e. foot-ground contact time)
agilityPredictors.stanceTime.meanLeftStanceTime = mean(agilityParams.stanceTime.left);
agilityPredictors.stanceTime.meanRightStanceTime = mean(agilityParams.stanceTime.right);
agilityPredictors.stanceTime.meanBilateralStanceTime = mean([agilityParams.stanceTime.left agilityParams.stanceTime.right]);
agilityPredictors.stanceTime.minLeftStanceTime = min(agilityParams.stanceTime.left);
agilityPredictors.stanceTime.minRightStanceTime = min(agilityParams.stanceTime.right);
agilityPredictors.stanceTime.minBilateralStanceTime = min([agilityParams.stanceTime.left agilityParams.stanceTime.right]);
agilityPredictors.stanceTime.maxLeftStanceTime = max(agilityParams.stanceTime.left);
agilityPredictors.stanceTime.maxRightStanceTime = max(agilityParams.stanceTime.right);
agilityPredictors.stanceTime.maxBilateralStanceTime = max([agilityParams.stanceTime.left agilityParams.stanceTime.right]);
agilityPredictors.stanceTime.totalLeftStanceTime = sum(agilityParams.stanceTime.left);
agilityPredictors.stanceTime.totalRightStanceTime = sum(agilityParams.stanceTime.right);

% Step time (i.e. time in between left and right steps)
agilityPredictors.stepTime.meanL2RStepTime = mean(agilityParams.stepTime.left);
agilityPredictors.stepTime.meanR2LStepTime = mean(agilityParams.stepTime.right);
agilityPredictors.stepTime.meanBilateralStepTime = mean([agilityParams.stepTime.left agilityParams.stepTime.right]);
agilityPredictors.stepTime.minL2RStepTime = min(agilityParams.stepTime.left);
agilityPredictors.stepTime.minR2LStepTime = min(agilityParams.stepTime.right);
agilityPredictors.stepTime.minBilateralStepTime = min([agilityParams.stepTime.left agilityParams.stepTime.right]);
agilityPredictors.stepTime.maxL2RStepTime = max(agilityParams.stepTime.left);
agilityPredictors.stepTime.maxR2LStepTime = max(agilityParams.stepTime.right);
agilityPredictors.stepTime.maxBilateralStepTime = max([agilityParams.stepTime.left agilityParams.stepTime.right]);

% Swing time (i.e. time in between heel strikes of the same foot)
agilityPredictors.swingTime.meanLeftSwingTime = mean(agilityParams.swingTime.left);
agilityPredictors.swingTime.meanRightSwingTime = mean(agilityParams.swingTime.right);
agilityPredictors.swingTime.meanBilateralSwingTime = mean([agilityParams.swingTime.left agilityParams.swingTime.right]);
agilityPredictors.swingTime.minLeftSwingTime = min(agilityParams.swingTime.left);
agilityPredictors.swingTime.minRightSwingTime = min(agilityParams.swingTime.right);
agilityPredictors.swingTime.minBilateralSwingTime = min([agilityParams.swingTime.left agilityParams.swingTime.left]);
agilityPredictors.swingTime.maxLeftSwingTime = max(agilityParams.swingTime.left);
agilityPredictors.swingTime.maxRightSwingTime = max(agilityParams.swingTime.right);
agilityPredictors.swingTime.maxBilateralSwingTime = max([agilityParams.swingTime.left agilityParams.swingTime.right]);

% Double support time (i.e. time both feet touch the ground, i.e. time between heel strike of
% front foot and toe off of other foot)
agilityPredictors.dsTime.totalDoubleSupportTime = sum([agilityParams.dsTime.left_to_right(agilityParams.dsTime.left_to_right > 0) ...
    agilityParams.dsTime.right_to_left(agilityParams.dsTime.right_to_left > 0)]);
% rowIndex = find(completionTimes.trial == trialName);
% totalTime = table2array(completionTimes(rowIndex, {'total'}));
% totalTime = completion



agilityPredictors.dsTime.doubleSupportPercentage = agilityPredictors.dsTime.totalDoubleSupportTime / completionTime;

% Step length
agilityPredictors.stepLength.meanRightStepLength = mean(agilityParams.stepLength.right);
agilityPredictors.stepLength.meanLeftStepLength = mean(agilityParams.stepLength.left);
agilityPredictors.stepLength.meanBilateralStepLength = mean([agilityParams.stepLength.left agilityParams.stepLength.right]);
agilityPredictors.stepLength.minRightStepLength = min(agilityParams.stepLength.right);
agilityPredictors.stepLength.minLeftStepLength = min(agilityParams.stepLength.left);
agilityPredictors.stepLength.minBilateralStepLength = min([agilityParams.stepLength.left agilityParams.stepLength.right]);
agilityPredictors.stepLength.maxRightStepLength = max(agilityParams.stepLength.right);
agilityPredictors.stepLength.maxLeftStepLength = max(agilityParams.stepLength.left);
agilityPredictors.stepLength.maxBilateralStepLength = max([agilityParams.stepLength.left agilityParams.stepLength.right]);

% Gait speed
agilityPredictors.runningGaitSpeed = agilityPredictors.stepLength.meanBilateralStepLength/agilityPredictors.stepTime.meanBilateralStepTime;

clearvars rowIndex totalTime

% jointAngles.kneeAngle(:,1) %right
% jointAngles.kneeAngle(:,2) %left
edges = 0:60:180;
agilityParams.kneeAnglesLowMidHigh(:,1) = histcounts(jointAngles.kneeAngle(:,1), edges)'; %right knee
agilityParams.kneeAnglesLowMidHigh(:,2) = histcounts(jointAngles.kneeAngle(:,2), edges)'; %left knee

agilityParams.hipAnglesLowMidHigh(:,1) = histcounts(jointAngles.hipAngle(:,1), edges)'; %right knee
agilityParams.hipAnglesLowMidHigh(:,2) = histcounts(jointAngles.hipAngle(:,2), edges)'; %left knee

agilityParams.ankleAnglesLowMidHigh(:,1) = histcounts(jointAngles.ankleAngle(:,1), edges)'; %right knee
agilityParams.ankleAnglesLowMidHigh(:,2) = histcounts(jointAngles.ankleAngle(:,2), edges)'; %left knee

clearvars edges

%%

agilityPredictors.rangeOfMotion.knee.left = max(jointAngles.kneeAngle(:,2)) - min(jointAngles.kneeAngle(:,2));
agilityPredictors.rangeOfMotion.knee.right = max(jointAngles.kneeAngle(:,1)) - min(jointAngles.kneeAngle(:,1));
agilityPredictors.rangeOfMotion.hip.left = max(jointAngles.hipAngle(:,2)) - min(jointAngles.hipAngle(:,2));
agilityPredictors.rangeOfMotion.hip.right = max(jointAngles.hipAngle(:,1)) - min(jointAngles.hipAngle(:,1));
agilityPredictors.rangeOfMotion.ankle.left = max(jointAngles.ankleAngle(:,2)) - min(jointAngles.ankleAngle(:,2));
agilityPredictors.rangeOfMotion.ankle.right = max(jointAngles.ankleAngle(:,1)) - min(jointAngles.ankleAngle(:,1));

agilityPredictors.rangeOfMotion.knee.diff = abs(agilityPredictors.rangeOfMotion.knee.left - agilityPredictors.rangeOfMotion.knee.right);
agilityPredictors.rangeOfMotion.hip.diff = abs(agilityPredictors.rangeOfMotion.hip.left - agilityPredictors.rangeOfMotion.hip.right);
agilityPredictors.rangeOfMotion.ankle.diff = abs(agilityPredictors.rangeOfMotion.ankle.left - agilityPredictors.rangeOfMotion.ankle.right);

%%

%knee range of motion
agilityPredictors.disbalance.knee.LowRangeDisbalance = ...
    2*((abs(agilityParams.kneeAnglesLowMidHigh(1,1)-mean(agilityParams.kneeAnglesLowMidHigh(1,:))))/sum(agilityParams.kneeAnglesLowMidHigh(1,:)));
agilityPredictors.disbalance.knee.MidRangeDisbalance = ...
    2*((abs(agilityParams.kneeAnglesLowMidHigh(2,1)-mean(agilityParams.kneeAnglesLowMidHigh(2,:))))/sum(agilityParams.kneeAnglesLowMidHigh(2,:)));
agilityPredictors.disbalance.knee.HighRangeDisbalance = ...
    2*((abs(agilityParams.kneeAnglesLowMidHigh(3,1)-mean(agilityParams.kneeAnglesLowMidHigh(3,:))))/sum(agilityParams.kneeAnglesLowMidHigh(3,:)));

agilityPredictors.disbalance.knee.LowRangeDisbalance(isnan(agilityPredictors.disbalance.knee.LowRangeDisbalance))=0;
agilityPredictors.disbalance.knee.MidRangeDisbalance(isnan(agilityPredictors.disbalance.knee.MidRangeDisbalance))=0;
agilityPredictors.disbalance.knee.HighRangeDisbalance(isnan(agilityPredictors.disbalance.knee.HighRangeDisbalance))=0;

%hip range of motion
agilityPredictors.disbalance.hip.LowRangeDisbalance = ...
    2*((abs(agilityParams.hipAnglesLowMidHigh(1,1)-mean(agilityParams.hipAnglesLowMidHigh(1,:))))/sum(agilityParams.hipAnglesLowMidHigh(1,:)));
agilityPredictors.disbalance.hip.MidRangeDisbalance = ...
    2*((abs(agilityParams.hipAnglesLowMidHigh(2,1)-mean(agilityParams.hipAnglesLowMidHigh(2,:))))/sum(agilityParams.hipAnglesLowMidHigh(2,:)));
agilityPredictors.disbalance.hip.HighRangeDisbalance = ...
    2*((abs(agilityParams.hipAnglesLowMidHigh(3,1)-mean(agilityParams.hipAnglesLowMidHigh(3,:))))/sum(agilityParams.hipAnglesLowMidHigh(3,:)));

agilityPredictors.disbalance.hip.LowRangeDisbalance(isnan(agilityPredictors.disbalance.hip.LowRangeDisbalance))=0;
agilityPredictors.disbalance.hip.MidRangeDisbalance(isnan(agilityPredictors.disbalance.hip.MidRangeDisbalance))=0;
agilityPredictors.disbalance.hip.HighRangeDisbalance(isnan(agilityPredictors.disbalance.hip.HighRangeDisbalance))=0;

%ankle range of motion
agilityPredictors.disbalance.ankle.LowRangeDisbalance = ...
    2*((abs(agilityParams.ankleAnglesLowMidHigh(1,1)-mean(agilityParams.ankleAnglesLowMidHigh(1,:))))/sum(agilityParams.ankleAnglesLowMidHigh(1,:)));
agilityPredictors.disbalance.ankle.MidRangeDisbalance = ...
    2*((abs(agilityParams.ankleAnglesLowMidHigh(2,1)-mean(agilityParams.ankleAnglesLowMidHigh(2,:))))/sum(agilityParams.ankleAnglesLowMidHigh(2,:)));
agilityPredictors.disbalance.ankle.HighRangeDisbalance = ...
    2*((abs(agilityParams.ankleAnglesLowMidHigh(3,1)-mean(agilityParams.ankleAnglesLowMidHigh(3,:))))/sum(agilityParams.ankleAnglesLowMidHigh(3,:)));

agilityPredictors.disbalance.ankle.LowRangeDisbalance(isnan(agilityPredictors.disbalance.ankle.LowRangeDisbalance))=0;
agilityPredictors.disbalance.ankle.MidRangeDisbalance(isnan(agilityPredictors.disbalance.ankle.MidRangeDisbalance))=0;
agilityPredictors.disbalance.ankle.HighRangeDisbalance(isnan(agilityPredictors.disbalance.ankle.HighRangeDisbalance))=0;

agilityPredictors.speed.meanSpeed = mean(velocityCoM(idx_start:idx_end));
agilityPredictors.speed.maxSpeed = max(velocityCoM(idx_start:idx_end));
agilityPredictors.acceleration.meanAcceleration = mean(accelerationCoM(idx_start:idx_end));
agilityPredictors.acceleration.minAcceleration = min(accelerationCoM(idx_start:idx_end));
agilityPredictors.acceleration.maxAcceleration = max(accelerationCoM(idx_start:idx_end));

accelerationCoM_filt = accelerationCoM(idx_start:idx_end);
agilityPredictors.acceleration.meanPosAcceleration = mean(accelerationCoM_filt(accelerationCoM_filt > 0));
agilityPredictors.acceleration.meanNegAcceleration = mean(accelerationCoM_filt(accelerationCoM_filt < 0));
%%

% accelerationCoM_filt = accelerationCoM(idx_start:idx_end);
% mean(accelerationCoM_filt(accelerationCoM_filt < 0))
% % mean(accelerationCoM__)
% mean(accelerationCoM(accelerationCoM(idx_start:idx_end) < 0))
%%

% fprintf("Mean Speed: %.2f m/s | Mean Acceleration: %.2f m/s^2 \n" + ...
%     "Max Speed: %.2f m/s | Max Acceleration: %.2f m/s^2", ...
%     agilityPredictors.speed.meanSpeed,agilityPredictors.acceleration.meanAcceleration, ...
%     agilityPredictors.speed.maxSpeed,agilityPredictors.acceleration.maxAcceleration);

% clearvars -except agility_parameters subjectNum trialName trialDir agilityPredictors testResults anthropoData test completionTime

% clear variables
% load("agility_parameters.mat");

%% Write all KPIs to Excel
clc

if isfile("D:\GP\modelData.xlsx")
    existingData = readtable("D:\GP\modelData.xlsx");
else
    % Make N by 2 matrix of fieldname + value type
    variable_names_types = [["PID", "double"]; ...
                ["age","double"], ...
                ["weight", "double"]; ...
                ["height", "double"]; ...
                ["sittingHeight", "double"]; ...
                ["CODTest", "string"]; ...
                ["direction", "string"]; ...
                ["CODTestTime", "double"]; ...
                ["sprint10m1st", "double"]; ...
                ["sprint10m2nd", "double"]; ...
                ["sprint10mAvg", "double"]; ...
                ["sprint30m1st", "double"]; ...
                ["sprint30m2nd", "double"]; ...
                ["sprint30mAvg", "double"]; ...
                ["AFCT505Left", "double"]; ...
                ["AFCT505Right", "double"]; ...
                ["AFCT505Avg", "double"]; ...
                ["AFCTFullLeft", "double"]; ...
                ["AFCTFullRight", "double"]; ...
                ["AFCTFullAvg", "double"]; ...
			    ["Left505", "double"]; ...
                ["Right505", "double"]; ...
                ["Avg505", "double"]; ...
                ["arrowheadLeft", "double"]; ...
                ["arrowheadRight", "double"]; ...
                ["arrowheadAvg", "double"]; ...
                ["endurance1000m", "double"]; ...
                ["endurance400m", "double"]; ...
			    ["meanLeftStanceTime", "double"]; ...
			    ["meanRightStanceTime", "double"]; ...
			    ["meanBilateralStanceTime", "double"]; ...
			    ["minLeftStanceTime", "double"]; ...
			    ["minRightStanceTime", "double"]; ...
			    ["minBilateralStanceTime", "double"]; ...
                ["maxLeftStanceTime", "double"]; ...
			    ["maxRightStanceTime", "double"]; ...
			    ["maxBilateralStanceTime", "double"]; ...
                ["totalLeftStanceTime", "double"]; ...
                ["totalRightStanceTime", "double"]; ...
                ["meanL2RStepTime", "double"]; ...
                ["meanR2LStepTime", "double"]; ...
                ["meanBilateralStepTime", "double"]; ...
                ["minL2RStepTime", "double"]; ...
                ["minR2LStepTime", "double"]; ...
                ["minBilateralStepTime", "double"]; ...
                ["maxL2RStepTime", "double"]; ...
                ["maxR2LStepTime", "double"]; ...
                ["maxBilateralStepTime", "double"]; ...
                ["meanLeftSwingTime", "double"]; ...
                ["meanRightSwingTime", "double"]; ...
                ["meanBilateralSwingTime", "double"]; ...
                ["minLeftSwingTime", "double"]; ...
                ["minRightSwingTime", "double"]; ...
                ["minBilateralSwingTime", "double"]; ...
                ["maxLeftSwingTime", "double"]; ...
                ["maxRightSwingTime", "double"]; ...
                ["maxBilateralSwingTime", "double"]; ...
                ["totalDoubleSupportTime", "double"]; ...
                ["doubleSupportPercentage", "double"]; ...
                ["meanRightStepLength", "double"]; ...
                ["meanLeftStepLength", "double"]; ...
                ["meanBilateralStepLength", "double"]; ...
                ["minRightStepLength", "double"]; ...
                ["minLeftStepLength", "double"]; ...
                ["minBilateralStepLength", "double"]; ...
                ["maxRightStepLength", "double"]; ...
                ["maxLeftStepLength", "double"]; ...
                ["maxBilateralStepLength", "double"]; ...
                ["runningGaitSpeed", "double"]; ...
                ["kneeLowRangeDisbalance", "double"]; ...
                ["kneeMidRangeDisbalance", "double"]; ...
                ["kneeHighRangeDisbalance", "double"]; ...
                ["hipLowRangeDisbalance", "double"]; ...
                ["hipMidRangeDisbalance", "double"]; ...
                ["hipHighRangeDisbalance", "double"]; ...
                ["ankleLowRangeDisbalance", "double"]; ...
                ["ankleMidRangeDisbalance", "double"]; ...
                ["ankleHighRangeDisbalance", "double"]; ...
                ["kneeRangeOfMotionLeft", "double"]; ...
                ["kneeRangeOfMotionRight", "double"]; ...
                ["kneeRangeOfMotionDiff", "double"]; ...
                ["hipRangeOfMotionLeft", "double"]; ...
                ["hipRangeOfMotionRight", "double"]; ...
                ["hipRangeOfMotionDiff", "double"]; ...
                ["ankleRangeOfMotionLeft", "double"]; ...
                ["ankleRangeOfMotionRight", "double"]; ...
                ["ankleRangeOfMotionDiff", "double"]; ...
                ["meanLeftTakeoffDistance", "double"]; ...
                ["meanRightTakeoffDistance", "double"]; ...
                ["meanBilateralTakeoffDistance", "double"]; ...
                ["maxLeftTakeoffDistance", "double"]; ...
                ["maxRightTakeoffDistance", "double"]; ...
                ["maxBilateralTakeoffDistance", "double"]; ...
                ["meanSpeed", "double"]; ...
                ["maxSpeed", "double"]; ...
                ["meanAcceleration", "double"]; ...
                ["minAcceleration", "double"]; ...
                ["maxAcceleration", "double"]; ...
                ["meanPosAcceleration", "double"]; ...
                ["meanNegAcceleration", "double"]; ...
                ["vidDuration","double"]
                ];

    % Make table using fieldnames & value types from above
    existingData = table('Size',[0,size(variable_names_types,1)],... 
	    'VariableNames', variable_names_types(:,1),...
	    'VariableTypes', variable_names_types(:,2));

end

% Create a new row of data
newRow = { ...
    subjectNum, ...
    ...
    ... anthropometric data
    anthropoData.age, ...
    anthropoData.weight, ...
    anthropoData.height, ...
    anthropoData.sittingHeight, ...
    ...
    ... information regarding test that is analyzed
    trialName, ...
    trialDir, ...
    completionTime, ...
    ...
    ... test results for all tests (i.e. sprint, agility FCT, 505, arrowhead, endurance)
    testResults.sprintTest10m1st, ...
    testResults.sprintTest10m2nd, ...
    testResults.sprintTest10mAvg, ...
    testResults.sprintTest30m1st, ...
    testResults.sprintTest30m2nd, ...
    testResults.sprintTest30mAvg, ...
    testResults.AFCT505Left, ...
    testResults.AFCT505Right, ...
    testResults.AFCT505Avg, ...
    testResults.AFCTFullLeft, ...
    testResults.AFCTFullRight, ...
    testResults.AFCTFullAvg, ...
    testResults.Left505, ...
    testResults.Right505, ...
    testResults.Avg505, ...
    testResults.arrowheadLeft, ...
    testResults.arrowheadRight, ...
    testResults.arrowheadAvg, ...
    testResults.endurance1000m, ...
    testResults.endurance400m, ...
    ...
    ... stance time KPIs
    agilityPredictors.stanceTime.meanLeftStanceTime, ...
    agilityPredictors.stanceTime.meanRightStanceTime, ...
    agilityPredictors.stanceTime.meanBilateralStanceTime, ...
    agilityPredictors.stanceTime.minLeftStanceTime, ...
    agilityPredictors.stanceTime.minRightStanceTime, ...
    agilityPredictors.stanceTime.minBilateralStanceTime, ...
    agilityPredictors.stanceTime.maxLeftStanceTime, ...
    agilityPredictors.stanceTime.maxRightStanceTime, ...
    agilityPredictors.stanceTime.maxBilateralStanceTime, ...
    agilityPredictors.stanceTime.totalLeftStanceTime, ...
    agilityPredictors.stanceTime.totalRightStanceTime, ...
    ...
    ... step time KPIs
    agilityPredictors.stepTime.meanL2RStepTime, ...
    agilityPredictors.stepTime.meanR2LStepTime, ...
    agilityPredictors.stepTime.meanBilateralStepTime, ...
    agilityPredictors.stepTime.minL2RStepTime, ...
    agilityPredictors.stepTime.minR2LStepTime, ...
    agilityPredictors.stepTime.minBilateralStepTime, ...
    agilityPredictors.stepTime.maxL2RStepTime, ...
    agilityPredictors.stepTime.maxR2LStepTime, ...
    agilityPredictors.stepTime.maxBilateralStepTime, ...
    ...
    ... swing time KPIs
    agilityPredictors.swingTime.meanLeftSwingTime, ...
    agilityPredictors.swingTime.meanRightSwingTime, ...
    agilityPredictors.swingTime.meanBilateralSwingTime, ...
    agilityPredictors.swingTime.minLeftSwingTime, ...
    agilityPredictors.swingTime.minRightSwingTime, ...
    agilityPredictors.swingTime.minBilateralSwingTime, ...
    agilityPredictors.swingTime.maxLeftSwingTime, ...
    agilityPredictors.swingTime.maxRightSwingTime, ...
    agilityPredictors.swingTime.maxBilateralSwingTime, ...
    ...
    ... double support time KPIs
    agilityPredictors.dsTime.totalDoubleSupportTime, ...
    agilityPredictors.dsTime.doubleSupportPercentage, ...
    ...
    ... step length KPIs
    agilityPredictors.stepLength.meanRightStepLength, ...
    agilityPredictors.stepLength.meanLeftStepLength, ...
    agilityPredictors.stepLength.meanBilateralStepLength, ...
    agilityPredictors.stepLength.minRightStepLength, ...
    agilityPredictors.stepLength.minLeftStepLength, ...
    agilityPredictors.stepLength.minBilateralStepLength, ...
    agilityPredictors.stepLength.maxRightStepLength, ...
    agilityPredictors.stepLength.maxLeftStepLength, ...
    agilityPredictors.stepLength.maxBilateralStepLength, ...
    ...
    ... running gait speed KPI
    agilityPredictors.runningGaitSpeed, ...
    ...
    ... disbalance KPIs
    agilityPredictors.disbalance.knee.LowRangeDisbalance, ...
    agilityPredictors.disbalance.knee.MidRangeDisbalance, ...
    agilityPredictors.disbalance.knee.HighRangeDisbalance, ...
    agilityPredictors.disbalance.hip.LowRangeDisbalance, ...
    agilityPredictors.disbalance.hip.MidRangeDisbalance, ...
    agilityPredictors.disbalance.hip.HighRangeDisbalance, ...
    agilityPredictors.disbalance.ankle.LowRangeDisbalance, ...
    agilityPredictors.disbalance.ankle.MidRangeDisbalance, ...
    agilityPredictors.disbalance.ankle.HighRangeDisbalance, ...
    ...
    ... range of motion KPIs
    agilityPredictors.rangeOfMotion.knee.left, ...
    agilityPredictors.rangeOfMotion.knee.right, ...
    agilityPredictors.rangeOfMotion.knee.diff, ...
    agilityPredictors.rangeOfMotion.hip.left, ...
    agilityPredictors.rangeOfMotion.hip.right, ...
    agilityPredictors.rangeOfMotion.hip.diff, ...
    agilityPredictors.rangeOfMotion.ankle.left, ...
    agilityPredictors.rangeOfMotion.ankle.right, ...
    agilityPredictors.rangeOfMotion.ankle.diff, ...
    ...
    ... takeoff distance KPIs
    agilityPredictors.takeoffDistance.meanLeftTakeoffDistance, ...
    agilityPredictors.takeoffDistance.meanRightTakeoffDistance, ...
    agilityPredictors.takeoffDistance.meanBilateralTakeoffDistance, ...
    agilityPredictors.takeoffDistance.maxLeftTakeoffDistance, ...
    agilityPredictors.takeoffDistance.maxRightTakeoffDistance, ...
    agilityPredictors.takeoffDistance.maxBilateralTakeoffDistance, ...
    ...
    ... speed and acceleration KPIs
    agilityPredictors.speed.meanSpeed, ...
    agilityPredictors.speed.maxSpeed, ...
    agilityPredictors.acceleration.meanAcceleration, ...
    agilityPredictors.acceleration.minAcceleration, ...
    agilityPredictors.acceleration.maxAcceleration, ...
    agilityPredictors.acceleration.meanPosAcceleration, ...
    agilityPredictors.acceleration.meanNegAcceleration, ...
    ...
    vid.Duration
    };

% Append the new row to the existing data
updatedData = [existingData; newRow];

% Write the updated data to the Excel file
writetable(updatedData, "D:\GP\modelData.xlsx");

% f = figure;
% f('WindowState','maximized');
% hold on
% grid on
% 
% ylabel('Y')
% xlabel('X')
% xlim([-9000 9000]);
% ylim([-10000 15000]);
% % view([90 0]); %plot ZY plane (side view)
% view([180 90]); %plot XY plane (topview)
% title('Running pattern')
% plot([-10000 10000],[0 0],'k-');
% 
% plot(centerofmass(:,1),centerofmass(:,2),'y.','MarkerSize',2);
% 
% for i = idx_start:idx_end
%     % plot(centerofmass(i,1),centerofmass(i,2)-thresholdX,'Color',[i/size(pos,1) (size(pos,1)-i)/size(pos,1) 1 0.2],'Marker','.','MarkerSize',2);
%     plot(centerofmass(i,1),centerofmass(i,2),'Color',lineColor,'Marker','.','MarkerSize',2);
%     % while waitforbuttonpress ~= 1
%     % end
% end
% 
% axis square equal

% while waitforbuttonpress ~= 1
% end

% catch
%     continue
% end 
end
end

winopen("D:\GP\modelData.xlsx")


toc