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
    nOpenTriplets = sum(sum(A'^2)) - trace(A'^2);
    nClosedTriplets = trace(A'^3);
    clusteringCoefficient = nClosedTriplets/nOpenTriplets;

    % -----------Reciprocity-----------
    reciprocity = (1/nEdges)*trace(A'^2);
    
    % -----------Average Path Length-----------
    pathLengths = distances(G);
    pathLengths(~isfinite(pathLengths)) = 0;
    pathLengths = sum(pathLengths, 'all');
    APL = pathLengths/(nAgents*(nAgents - 1));

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
    metricsCell{5, 2} = APL;

    clearvars -except metricsCell
end