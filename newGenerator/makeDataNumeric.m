% Non-numeric data processing
% Some of the data are non-numeric. To use these in this model, they need
% to be converted to numeric data. Some (such as income and education) can
% be ranked from numbers 1 to 7 and 1 to 3 respectively. Race can be 
% converted into a logical array. This section processes this data into
% numerical values and reconstructs the data using these values.

function numericMatrix = makeDataNumeric(xyCoordinate, inData, mapSize)
    if size(inData, 2) ~= 14
        error('Some columns may be missing');
    end

    USA_data_temp = zeros(mapSize, 17); % Preallocate temporary variable for data
    raceEducationIncome = table2array(inData(:, [4, 8, 9]));  % Separate non-numeric data and place in table for comparison
    raceEducationIncome_temp = zeros(mapSize, 6);   % Set temporary variable for storage
    tempEdu = ["highschool", 1; "somecollege", 2; "collegeplus", 3];    % Set data conversion reference table
    tempIncome = ["X0", 0; "X1", 1; "X2", 2; "X3", 3; "X4", 4; "X5", 5; "X6", 6; "X7", 7]; % Set data conversion reference table
    tempRace = ["BLA", "SPA", "WHI", "OTH"]; % Set data conversion reference table
    
    for n = 1:height(inData)
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
    USA_data_temp(:, [1:3, 8:10, 13:17]) = table2array(inData(:, [1:3, 5:7, 10:14])); % Insert numerical data into array
    USA_data_temp(:, [4:7, 11, 12]) = [raceEducationIncome_temp(:, 1:4), raceEducationIncome_temp(:, 5), raceEducationIncome_temp(:, 6)]; % Insert processed data into array
    USA_data_temp(:, 18:19) = xyCoordinate;    % Insert x and y-coordinates
    USA_data_temp = USA_data_temp(:, 2:19); % Remove the ID value (as it is arbitrary)

    numericMatrix = USA_data_temp;
end