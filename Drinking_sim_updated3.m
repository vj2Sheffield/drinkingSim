% Suggestions:
% 1. Use version control (git) and annotate m-files
% 2. Come to meetings with Powerpoint that sets out your progress
% and thinking (including on method and results).
% 3. Use random number generators in Matlab stats toolbox to
% generate reaches. Need to account for scaling on Hamill. Check
% output using histograms.
% 4. Use graph theory utilities in the Matlab graph and network
% algoriths toolbox to calculate all the network metrics
% 5. Work with equally weighted population. Generate results for
% 10, 100 and 1000 agents.
% 6. Allow the whole weight vector to be easily manipulated, e.g.
% wghts=[wSex, wAge, wEthnicity]; normalise it so it sums to 1
% wSex = 10; Make the script so it generates all the results from this.
% 7. Use easily interpretable variable names, e.g. numOfAgents


% 16-Dec-19 Took Robin through script. Changed some variable names. 
%           Changed weight vector for easy modification
% 17-Dec-19 Changed distance metric to Euclidean. Put mfiles in git.
%           Improved on annotation and added in row randomiser for data table.
%           Changed poisson distribution to pre-made matlab function and began
%           working on gamma distribution
% 18-Dec-19 Finished working on gamma distribution. Fixed git stuff.
%           Removed old 'metrics' section and began to construct improved one using
%           inbuilt Matlab functionality.

% Housekeeping
clear, clc
if ~exist('Social Spheres Plots', 'dir')    % Check if directory exists
       mkdir 'Social Spheres Plots'         % For use in plot storage
end

% Declare modifiable variables
mapSize = 1000; % Size of the grid on which agents are placed
nAgents = 10;   % Number of agents used in simulation

wSex = 1;  % Weight parameters as named
wAge = 1;   % Age
wBlack = 1; % Race - black
wHispanic = 1;  % Race - hispanic
wWhite = 1; % Race - white
wOther = 1; % Race - other
wEmpStatus = 1; % Employment status 
wParentStatus = 1;  % Parenthood status
wMaritalStatus = 1; % Marital status
wEdu = 1;   % Level of education
wIncome = 1;    % Income bracket
wDrinkStatus = 1;   % Drinking status
wDrinkFrequency = 1;    % Drink frequency
wHED = 1;   % Heavy episodic drinking
wGPD = 1;   % GPD
wDrinkPerMonth = 1; % Drinks per month
xCoord = 1; % Spatial x-coordinate
yCoord = 1; % Spatial y-coordinate
   
w_p = [wSex, wAge, wBlack, wHispanic, wWhite, wOther, wEmpStatus, ...
       wParentStatus, wMaritalStatus, wEdu, wIncome, wDrinkStatus, ...
       wDrinkFrequency, wHED, wGPD, wDrinkPerMonth, xCoord, yCoord]; % 18 variables

w_p = w_p(:)/sum(w_p);  % Normalised weights of parameters

clear wSex wAge wBlack wHispanic wWhite wOther wEmpStatus wParentStatus...  % Clear the above (unneeded)
    wMaritalStatus wEdu wIncome wDrinkStatus wHED wGPD wDrinkPerMonth...
    xCoord yCoord wDrinkFrequency

%% Preallocate variables 
socialReachDistribution = readmatrix('poiss_RNG.csv');    % Heterogeneous social reach distribution
x_y_coord = readmatrix('xy_coord_RNG.csv'); % Get x and y coordinates
USA_data = readtable('USA1000.csv');    % Read in data table

USA_data_temp = zeros(mapSize, 17); % Preallocate temporary variable for data
adjacencyMatrix = zeros(nAgents, nAgents);  % Preallocate variable for adjacency matrix
distance_ij = zeros(1, 18); % Preallocate variable for distance calculations
raceEducationIncome = table2array(USA_data(:, [4, 8, 9]));  % Separate non-numeric data and place in table for comparison
raceEducationIncome_temp = zeros(mapSize, 6);   % Set temporary variable for storage
socialDistances = zeros(nAgents);   % Preallocate temporary variable for social distance measurements

tempEdu = ["highschool", 1; "somecollege", 2; "collegeplus", 3];    % Set data conversion reference table
tempIncome = ["X0", 0; "X1", 1; "X2", 2; "X3", 3; "X4", 4; "X5", 5; "X6", 6; "X7", 7]; % Set data conversion reference table
tempRace = ["BLA", "SPA", "WHI", "OTH"]; % Set data conversion reference table

headers = {'Sex', 'Age', 'BLA', 'SPA', 'WHI',...    % Set headers for future use
    'OTH', 'Employment Status', 'Parenthood Status',... 
    'Marital Status', 'Education', 'Income', 'Drinking Status',...
    'Drink Frequency', 'HED', 'GPD', 'Drinks per Month', ...
    'X-coord', 'Y-coord'};

figure; 
histogram(socialReachDistribution); % Plot histogram of social reaches
title('Histogram of Social reaches');   % Title
xlabel('Size of Social Reach'); ylabel('Number of Agents'); % Plot labels
%% Non-numeric data processing
% Some of the data is non-numeric. To use these in this model, they need
% to be converted to numeric data. Some (such as income and education) can
% be ranked from numbers 1 to 7 and 1 to 3 respectively. Race can be 
% converted into a logical array. This section processes this data into
% numerical values and reconstructs the data using these values.

for n = 1:height(USA_data)
    mask_edu = ismember(tempEdu(:, 1), raceEducationIncome(n, 2));  % Find logical array matching reference table (as above)
    index_edu = find(mask_edu == 1);    % Find which values are true
    raceEducationIncome_temp(n, 5) = index_edu; % Store the value for later use
    
    % Repeat for income
    mask_inc = ismember(tempIncome(:, 1), raceEducationIncome(n, 3));
    index_inc = find(mask_inc == 1);
    raceEducationIncome_temp(n, 6) = index_inc - 1;
    
    % Repeat for race
    mask_rac = ismember(tempRace, raceEducationIncome(n, 1));
    raceEducationIncome_temp(n, 1:4) = mask_rac;
end

% Table reconstruction
USA_data_temp(:, [1:3, 8:10, 13:17]) = table2array(USA_data(:, [1:3, 5:7, 10:14])); % Insert numerical data into array
USA_data_temp(:, [4:7, 11, 12]) = [raceEducationIncome_temp(:, 1:4), raceEducationIncome_temp(:, 5), raceEducationIncome_temp(:, 6)]; % Insert processed data into array
USA_data_temp(:, 18:19) = x_y_coord;    % Insert x and y-coordinates
USA_data_temp = USA_data_temp(:, 2:19); % Remove the ID value (as it is arbitrary)

% Below: check if the reshuffled data exists for the number of agents
% needed. If not, write it using function randomiseRows().
if exist('RandUSA_data.csv', 'file') ~= 2 || length(readmatrix('RandUSA_data.csv')) ~= nAgents
    randomiseRows(USA_data_temp, nAgents);  
end

USA_data2 = readmatrix('RandUSA_data.csv'); % Read randomised data in
USdat_stdev = std(USA_data2(1:nAgents, :)); % Get standard deviations for parameters

% Cleaning Workspace
clear T_e T_i T_r n USA_data_temp REI REI_temp tempIncome tempEdu ...
      mask_edu mask_rac mask_inc index_edu index_inc USA_data tempRace ...
      raceEducationIncome raceEducationIncome_temp

%% Calculate social distances
% This section calculates the social distances between nodes using the 
% square of the sum of the Euclidean distances. If the social distance is
% less than the predetermined social reach, a connection is made and 
% stored in the adjacency matrix.

for k = 1:nAgents
    agent_i = USA_data2(k, :)'; % Get agent i details
    for m = 1:nAgents
        if k~= m
            dTemp_ij = 0;   % Set temp variable to 0
            agent_j = USA_data2(m, :)'; % Get agent j details
            for n = 1:size(USA_data2, 2)       % Find Euclidean distance
                distance_ij(n) = (w_p(n)*(agent_i(n) - agent_j(n))^2)/USdat_stdev(n); % Calculate social distance between i and j
                dTemp_ij = distance_ij(n) + dTemp_ij;   % Get sum of social distances
            end
            socialDistances(k, m) = sqrt(dTemp_ij); % Store social distance

            if socialDistances(k, m) < socialReachDistribution(m) && socialDistances(m, k) < socialReachDistribution(k) % If a connection is made (undirected)
                adjacencyMatrix(k, m) = 1;  % Set adjacency matrix element to 1
            end
        end
    end
end

%% Metrics
adjGraph = graph(adjacencyMatrix);
adjDegree = degree(adjGraph);
adjEdges = numedges(adjGraph);

adjDensity = adjEdges/(nAgents*(nAgents - 1))*100;  % Calculate network centrality

%% Test Metrics
clc
A = [0 1 0 0 1 0;...
     0 0 1 1 0 0;...
     0 1 0 1 0 0;...
     0 1 1 0 0 0;...
     1 0 1 0 0 0;...
     0 0 1 0 0 0];
G = digraph(A);
nAgents = 6;
nEdges = numedges(G);
nNodes = numnodes(G);

% -----------Degree Centrality Test-----------
inDegree = centrality(G, 'indegree');
outDegree = centrality(G, 'outdegree');
totalDegree = inDegree + outDegree;

% -----------Katz Centrality Test-----------
alpha = 0.9*(1/max(eig(A)));
beta = 1;

katzCentrality = beta*(inv(eye(nAgents) - alpha*A'))*ones(nAgents, 1);
katzCentrality = round(katzCentrality, 3, 'significant');

% -----------Closeness Centrality Test-----------
pathLength = zeros(nAgents, nAgents);
pathLengthForAPL = zeros(nAgents, nAgents);
for n = 1:nAgents
    for m = 1:nAgents
        [~, pathLength(n, m)] = shortestpath(G, m, n);
        pathLengthForAPL(n, m) = pathLength(n, m);
        pathLength(n, m) = 1/pathLength(n, m);
    end
end
pathLength(pathLength == Inf) = 0;
closenessCentrality = sum(pathLength)'*(1/(nAgents-1));
inCloseness = centrality(G, 'incloseness');
outCloseness = centrality(G, 'outcloseness');

% -----------Betweenness Centrality Test-----------
betweenness = centrality(G, 'betweenness');

% -----------Transitivity-----------


% -----------Reciprocity-----------
count = 0;
for m = 1:nAgents
    for n = 1:nAgents
        if A(m, n) == A(n, m) && n ~= m && A(m, n) == 1
            count = count + 1;
        end
    end
end
reciprocity = count/nEdges;

% -----------Average Path Length-----------
tempAPL = 0;
validNodes = 0;
for m = 1:nAgents
    for n = 1:nAgents
        if pathLengthForAPL(n, m) ~= Inf && n ~= m
            tempAPL = pathLengthForAPL(n, m) + tempAPL;
            validNodes = validNodes + 1;
        end
    end
end
averagePathLength = tempAPL/(2*nAgents^2);

%% Putting it all together:
% This section plots the agents and their connections in a toroidal space.
% The connections are taken from the adjacency matrix (adjacencyMatrix)
% and plotted using setloop(). Connections that aren't made are also
% tracked in the count variable. The output of this section is a plot
% showing nodes and connections.

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

scatter(x_y_coord(1:nAgents, 1), x_y_coord(1:nAgents, 2), 'k', 'filled');   % Plot all agents
grid on
xlim([0 mapSize]); ylim([0 mapSize]);   % Set plot size

