% Housekeeping
clc, clear
tic

% Declare modifiable variables
mapSize = 1000; % Size of the grid on which agents are placed
nAgents = 1000;   % Number of agents used in simulation
nConnections = 12;  % Average number of connections needed
preferredDistribution = 'power';

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
wDrinkStatus = 1;       % Drinking status
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
% clc
% generateDistributions(preferredDistribution, lambda);
% queryPoint = nConnections/nAgents*100;
% 
% switch preferredDistribution
%     case 'poisson'
%         interTable = readmatrix('poissonInterpTable.csv');
%         interTableV = interTable(:, 1)';
%         interTableX = interTable(:, 2)';
%         if queryPoint > max(interTableX); queryPoint = max(interTableX); end
% 
%         lambda = interp1(interTableX, interTableV, queryPoint);
%         generateDistributions(preferredDistribution, lambda);
%         histDist = readmatrix('randomDistribution.csv');
%         histogram(histDist);
%         
%     case 'uniform'
%         interTable = readmatrix('uniformInterpTable.csv');
%         interTableV = interTable(:, 1)';
%         interTableX = interTable(:, 2)';
%         if queryPoint > max(interTableX); queryPoint = max(interTableX); end
% 
%         lambda = interp1(interTableX, interTableV, queryPoint);
%         generateDistributions(preferredDistribution, [lambda - 1, lambda + 1]);
%         histDist = readmatrix('randomDistribution.csv');
%         histogram(histDist);
%         
%     case 'power'
%         generateDistributions(preferredDistribution, [0.2, 18, 1.0986, 0.5]);
% end

%% Generate Random Distributions
clc
switch preferredDistribution
    case 'poisson'
        if nAgents == 10
            lambda = 9;
        elseif nAgents == 100
            lambda = 3.5;
        elseif nAgents == 1000
            lambda = 1.061;
        end
        generateDistributions(preferredDistribution, lambda);
                
    case 'uniform'
        if nAgents == 10
            upper = 20;
            lower = 18;
            range = [lower, upper];
        elseif nAgents == 100
            upper = 5;
            lower = 3.5;
            range = [lower, upper];
        elseif nAgents == 1000
            upper = 2.4;
            lower = 1;
            range = [lower, upper];
        end
        generateDistributions(preferredDistribution, range);
        
    case 'power'
        if nAgents == 10
            upper = 18;
            lower = 0.2;
            mu = 10*0.0592;
            sigma = 0.1;
            powerParameters = [lower, upper, mu, sigma];
        elseif nAgents == 100
            upper = 18;
            lower = 0.2;
            mu = 10*0.0592;
            sigma = 0.1;
            powerParameters = [lower, upper, mu, sigma];
        elseif nAgents == 1000
            upper = 18;
            lower = 0.2;
            mu = 10*0.0592;
            sigma = 0.1;
            powerParameters = [lower, upper, mu, sigma];
        end
        generateDistributions(preferredDistribution, powerParameters);
end

%% Read in Data Files
socialReachDistribution = readmatrix('randomDistribution.csv');    % Heterogeneous social reach distribution
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
degreeConnections = metricsCell{2,2};
degreeConnections = degreeConnections(:, 1);
histogram(degreeConnections);
xlabel('Number of Connections'); ylabel('Frequency'); title('Histogram of Connections');

G = digraph(adjacencyMatrix);
averageConnectionsMade = mean(centrality(G,'outdegree'));
fprintf('Average connections made per node: %f\n', averageConnectionsMade)

% %% Plotting it all together:
% % This section plots the agents and their connections in a toroidal space.
% % The connections are taken from the adjacency matrix (adjacencyMatrix)
% % and plotted using setloop(). The output of this section is a plot
% % showing nodes and connections.
% 
% figure
% for k = 1:nAgents
%     x1 = x_y_coord(k, 1); y1 = x_y_coord(k, 2); % Set point 1 as used in plot
%     for m = 1:nAgents 
%         if adjacencyMatrix(k, m) == 1 && m ~= k  % If a connection is made
%             x2 = x_y_coord(m, 1); y2 = x_y_coord(m, 2); % Set point 2 as used in plot
%             setloop(x1, x2, y1, y2, mapSize);   % Plot nodes and connections
%             hold on; 
%         end
%     end
% end
% 
% clear x1 x2 y1 y2 k m
% scatter(x_y_coord(1:nAgents, 1), x_y_coord(1:nAgents, 2), 'k', 'filled');   % Plot all agents
% grid on
% xlim([0 mapSize]); ylim([0 mapSize]);   % Set plot size
