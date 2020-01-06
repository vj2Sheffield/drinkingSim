function metricsCell = networkMetrics(A)
    % Test Metrics
    G = digraph(A);
    nAgents = numnodes(G);
    nEdges = numedges(G);

    % -----------Network Density-----------
    networkDensity = nEdges/(nAgents*(nAgents - 1));

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
    closenessCentrality = sum(pathLength)'*(1/(nAgents - 1));

    % -----------Betweenness Centrality Test-----------
    betweenness = centrality(G, 'betweenness');

    % -----------Transitivity-----------
    tempMatrix = zeros(100, 2);
    count = 0;
    for m = 1:nAgents
        for n = 1:nAgents
            [~, transPathLength] = shortestpath(G, m, n);
            if transPathLength == 1
                count = count + 1;
                tempMatrix(count, 1) = m;
                tempMatrix(count, 2) = n;
            end
        end
    end
    tempMatrix((count + 1):end, :) = [];

    count = 0;
    tempMatrix2 = zeros(100, 3);
    for m1 = 1:length(tempMatrix)
        for m2 = 1:length(tempMatrix)
            if tempMatrix(m1, 2) == tempMatrix(m2, 1)
                if tempMatrix(m1, 1) ~= tempMatrix(m2, 2)
                    count = count + 1;
                    tempMatrix2(count, 1:3) = [tempMatrix(m1, 1), tempMatrix(m2, 1:2)];
                end
            end
        end
    end
    tempMatrix2((count + 1):end, :) = [];

    nOpenTriplets = length(tempMatrix2);
    tempMatrix2 = sort(tempMatrix2, 2);

    [~, ~, ic] = unique(tempMatrix2, 'rows');
    a_counts = accumarray(ic, 1);
    nClosedTriplets = sum(a_counts(a_counts > 1));

    clusteringCoefficient = nClosedTriplets/nOpenTriplets;

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

    % -----------Return Cell-----------
    % Density, Centrality, Transitivity, Reciprocity, APL
    metricsCell = cell(5, 2);

    metricsCell{1, 1} = 'Density';
    metricsCell{1, 2} = networkDensity;

    metricsCell{2, 1} = 'Centrality';
    metricsCell{2, 2} = [totalDegree, katzCentrality, closenessCentrality,...
                         betweenness];

    metricsCell{3, 1} = 'Transitivity';
    metricsCell{3, 2} = clusteringCoefficient;

    metricsCell{4, 1} = 'Reciprocity';
    metricsCell{4, 2} = reciprocity;

    metricsCell{5, 1} = 'Average Path Length';
    metricsCell{5, 2} = averagePathLength;

    clearvars -except metricsCell
end