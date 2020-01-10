% Calculate social distances
% This section calculates the social distances between nodes using the 
% square of the sum of the Euclidean distances. If the social distance is
% less than the predetermined social reach, a connection is made and 
% stored in the adjacency matrix.

function [adjacencyMatrix, socialDistances] = getSocialDistances(nAgents, inData, USdat_stdev, socialReachDistribution, w_p)
    adjacencyMatrix = zeros(nAgents, nAgents);  % Preallocate variable for adjacency matrix
    distance_ij = zeros(1, 18); % Preallocate variable for distance calculations
    socialDistances = zeros(nAgents);   % Preallocate temporary variable for social distance measurements

    for k = 1:nAgents
        agent_i = inData(k, :)'; % Get agent i details
        for m = 1:nAgents
            if k~= m
                dTemp_ij = 0;   % Set temp variable to 0
                agent_j = inData(m, :)'; % Get agent j details
                for n = 1:size(inData, 2)       % Find Euclidean distance
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
end