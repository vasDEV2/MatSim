%Read paths
projPath = getenv('PROJECT');
simPath = getenv('SIM');
cityPath = getenv('MCITY');
savePath = getenv('PYTHON_SCRIPTS_DIR');
currPath = getenv('SCRIPT_DIR')



openProject(projPath);
%sceneName = 'VirtualMCity';
%[sceneImage, sceneRef] = helperGetSceneImage(sceneName);
%hFig = helperSelectSceneWaypoints(sceneImage, sceneRef);

sceneImage = imread(cityPath);
data = load('sim3d_SpatialReferences.mat');
sceneRef = data.spatialReference.VirtualMCity;
hFig = helperSelectSceneWaypoints(sceneImage, sceneRef);


while true
        if exist('wayPoints', 'var') % Check workspace
            break; % Waypoints found in workspace
        else
            % Pause briefly to avoid excessive CPU usage
            pause(0.1);
        end
end
%waypoints = readmatrix('/home/vasudevan/VDMGEN/SceneCameraRay53/SceneCameraRay/wp.csv');
fprintf('Waypoints are available. Proceeding...\n');

waypoints = cell2mat(wayPoints);
waypoints(:, 2) = waypoints(:, 2) * -1;
waypoints(:, 1) = waypoints(:, 1)-8;
waypoints(:, 2) = waypoints(:, 2)-3;

try
    sim(simPath); 
catch
    % If simulation is terminated manually, continue to next step
    warning('Simulation was terminated manually or failed. Proceeding to next step...');
    simOut = []; % Set simOut to an empty value, so the next part can run
end

% for i = 1:logsout.numElements
%     % Get the signal name and values
%     signalName = logsout{1}.Name; % Signal name
%     signalData = logsout{1}.Values; % Extract timeseries object
% 
%     % Extract time and data from the timeseries object
%     signalTime = signalData.Time; % Extract time
%     signalValues = signalData.Data; % Extract data
% 
%     % Save the signal data into a structured variable
%     try
%         assignin('base', [signalName '_time'], signalTime); % Save time
%         assignin('base', [signalName '_data'], signalValues); % Save data
%     catch
%         warning('Proceeding')
%     end
% 
%     % Optional: Display the extraction status
%     fprintf('Signal "%s" extracted as "%s_time" and "%s_data".\n', ...
%             signalName, signalName, signalName);
% end

signalTime = logsout{1}.Values.Time;
signalValues = logsout{1}.Values.Data;


INPUT_CODE = table(VX.time, VX.signals.values, VY.signals.values,R.signals.values,signalValues,omega.signals.values,lamda.signals.values,TAU.signals.values);
INPUT_CODE.Properties.VariableNames = {'TIME','VX','VY','R','Fz','omega','slip','TAU'};
writetable(INPUT_CODE, strcat(savePath,'/vdm_input.csv'));
S = table(STRA.time, STRA.signals.values);
S.Properties.VariableNames = {'TIME','STRA'};
writetable(S, strcat(savePath,'/steering.csv'));
OUTPUT_CODE = table(AX.time, AX.signals.values, AY.signals.values, AngACC.signals.values, X.signals.values, Y.signals.values, R.signals.values,PSI.signals.values);
OUTPUT_CODE.Properties.VariableNames = {'TIME','AX','AY','AngACC','X','Y','R','PSI'};
writetable(OUTPUT_CODE, strcat(savePath,'/output_vdm.csv'));
WAYPOINTSE = table(waypoints);
writetable(WAYPOINTSE, 'wp.csv')

% Define the name of the shell script file
fileName = strcat(currPath,'/flag.sh');

% Open the file for writing
fileID = fopen(fileName, 'w');

% Write the shell script content
fprintf(fileID, 'export MATLAB_F="true"');
% Close the file
fclose(fileID);

% Make the shell script executable
system(['chmod +x ', fileName]);
% setenv("MATLAB_F", "true");

disp('Simulation completed.');


quit;
