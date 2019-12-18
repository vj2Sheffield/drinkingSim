% This function takes a table with a length nAgents and randomises the
% rows, ensuring that the sum of each column is larger than zero. This is
% important because the standard deviations of the table columns are used
% in distance calculations, and if this is zero, all social distances
% become NaN 
% Victoria Johnson 
% 17/12/19

function randomiseRows(USA_data, nAgents)
    random_USA_data = USA_data(randperm(size(USA_data, 1)), :); % Randomise the rows
    ifSumColumnsZero = sum(random_USA_data(1:nAgents, :));  % Find the sum of each column up to a length of nAgents

    while sum(ifSumColumnsZero == 0) ~= 0   % If the sum of any of the columns is 0 (repeats until not true)
        random_USA_data = USA_data(randperm(size(USA_data, 1)), :); % Randomise again
        ifSumColumnsZero = sum(random_USA_data(1:nAgents, :));  % Find the sum of the columns again
    end

    writematrix(random_USA_data(1:nAgents, :), 'RandUSA_data.csv'); % Write the acceptable matrix to a csv
    clear ifSumColumnsZero random_USA_data  % Clear irrelevant variables.
end