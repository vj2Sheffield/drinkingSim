% Housekeeping
clear, clc

% Declare modifiable variables
mapSize = 1000; % Size of the grid on which agents are placed
nAgents = 10;   % Number of agents used in simulation
nConnections = 12;  % Average number of connections needed
preferredDistribution = 'poisson';

wSex = 1;               % Weight parameters as named
wAge = 1;               % Age
wBlack = 1;             % Race - black
wHispanic = 1;          % Race - hispanic
wWhite = 1;             % Race - white
wOther = 1;             % Race - other
wEmpStatus = 1;         % Employment status 
wParentStatus = 1;      % Parenthood status
wMaritalStatus = 1;     % Marital status
wEdu = 1;               % Level of education
wIncome = 1;            % Income bracket
wDrinkStatus = 1    ;   % Drinking status
wDrinkFrequency = 1;    % Drink frequency
wHED = 1;               % Heavy episodic drinking
wGPD = 1;               % GPD
wDrinkPerMonth = 1;     % Drinks per month
xCoord = 1;             % Spatial x-coordinate
yCoord = 1;             % Spatial y-coordinate
   
w_p = [wSex, wAge, wBlack, wHispanic, wWhite, wOther, wEmpStatus, ...
       wParentStatus, wMaritalStatus, wEdu, wIncome, wDrinkStatus, ...
       wDrinkFrequency, wHED, wGPD, wDrinkPerMonth, xCoord, yCoord]; % 18 variables

w_p = w_p(:)/sum(w_p);  % Normalised weights of parameters

clear wSex wAge wBlack wHispanic wWhite wOther wEmpStatus wParentStatus...  % Clear the above (unneeded)
    wMaritalStatus wEdu wIncome wDrinkStatus wHED wGPD wDrinkPerMonth...
    xCoord yCoord wDrinkFrequency

%% Generate Random Distributions
queryPoint = nConnections/nAgents*100;
if queryPoint > 100; queryPoint = 100; end

%% Read in Data Files
socialReachDistribution = readmatrix('poiss_RNG.csv');    % Heterogeneous social reach distribution
x_y_coord = readmatrix('xy_coord_RNG.csv'); % Get x and y coordinates
unprocessedUSA_data = readtable('USA1000.csv');    % Read in data table

%% Process Non-numeric Data
unprocessedUSA_data = makeDataNumeric(x_y_coord, unprocessedUSA_data, mapSize);
% Below: check if the reshuffled data exists for the number of agents
% needed. If not, write it using function randomiseRows().
if exist('RandUSA_data.csv', 'file') ~= 2 || length(readmatrix('RandUSA_data.csv')) ~= nAgents
    randomiseRows(unprocessedUSA_data, nAgents);  
end

%% Final Calculations and Transfer Data to Workspace Variables
USA_data = readmatrix('RandUSA_data.csv'); % Read randomised data in
USdat_stdev = std(USA_data(1:nAgents, :)); % Get standard deviations for parameters

%% Calculate Adjacency Matrix from Social Distances and get Network Metrics
clc
[adjacencyMatrix, socialDistances] = getSocialDistances(nAgents, USA_data, USdat_stdev, socialReachDistribution, w_p);
metricsCell = networkMetrics(adjacencyMatrix);

averageConnectionsMade = sum(adjacencyMatrix, 'all')/(2*nAgents);
fprintf('Average connections made per node: %f\n', averageConnectionsMade)

%% Metrics (more)
temp = metricsCell{2,2};
temp = temp(1:end, 1);
[count, bins] = hist(temp);
count = count/nAgents*100;
% bins = bins/(2*nAgents)*100;
frequencyMatrix = [bins; count];
plot(bins, count)
title('DoC - Exponential');
xlabel('Degree of Connectivity (% of total possible connections)'); 
ylabel('Frequency %');
grid on
% csvwrite('powerFrequency.csv', frequencyMatrix)  % Write to file
% sum(adjacencyMatrix,'all')

%% Plotting it all together:
% This section plots the agents and their connections in a toroidal space.
% The connections are taken from the adjacency matrix (adjacencyMatrix)
% and plotted using setloop(). Connections that aren't made are also
% tracked in the count variable. The output of this section is a plot
% showing nodes and connections.

figure
count = 0;  % Counting connections not made
for k = 1:nAgents
    for m = 1:nAgents 
        x2 = x_y_coord(m, 1); y2 = x_y_coord(m, 2); % Set point 2 as used in plot
        if m == k || adjacencyMatrix(k, m) == 0 % Check a connection is made and nodes aren't connecting with themselves
            count = count + 1;  % Count connections not made
        elseif adjacencyMatrix(k, m) == 1   % If a connection is made
            x1 = x_y_coord(k, 1); y1 = x_y_coord(k, 2); % Set point 1 as used in plot
            setloop(x1, x2, y1, y2, mapSize);   % Plot nodes and connections
            hold on; 
        end
    end
end

clear x1 x2 y1 y2 k m
scatter(x_y_coord(1:nAgents, 1), x_y_coord(1:nAgents, 2), 'k', 'filled');   % Plot all agents
grid on
xlim([0 mapSize]); ylim([0 mapSize]);   % Set plot size
