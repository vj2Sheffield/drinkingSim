% This function firstly calculates the distances between two points
% supplied in the function input. Then, if the distance is smaller than the
% radius of a circle with the length of half the map size, draw the
% connection between the points across the map. Otherwise, create a tile
% surface, find the shortest distance between points 1 and 2 on there, draw
% lines between all of them, then crop to the map size.
% Victoria Johnson
% December 2019

function setloop(x1, x2, y1, y2, limit)
    euclid_tile_distance_p1 = zeros(9, 1);  % Preallocating variable for tiles
    d1 = sqrt((x2 - x1)^2 + (y2 - y1)^2);   % Finding the Euclidean distance between points
        
    if d1 > sqrt((0.5*limit)^2 + (0.5*limit)^2) % If the distance between the points is larger than a circle with a radius of d1
        % Establish a 3x3 tile surface using the tiles() function
        temp_coord_p1 = tiles(x1, y1, limit);   % Tiles from point 1
        temp_coord_p2 = tiles(x2, y2, limit);   % Tiles from point 2

        for n = 1:9 % Find Pythagoras length between tiles  
            euclid_tile_distance_p1(n) = sqrt((x2 - temp_coord_p1(n, 1))^2 + (y2 - temp_coord_p1(n, 2))^2); % Find the distance between points on each tile
        end
        smallest_distance = min(euclid_tile_distance_p1);   % Get the smallest distance calculated above
        coord_idx = euclid_tile_distance_p1 == smallest_distance;   % Find the tile number with the shortest distance (binary logical array)
        smallest_coord_p1 = temp_coord_p1(coord_idx(:, 1) == 1, :); % As above (converstion to decimal)
        smallest_coord_p2 = temp_coord_p2(5, :);    % Set the other coordinate to the original point

        if length(smallest_coord_p1) > 1    % If there is more than one 'smallest point...
            smallest_coord_p1 = smallest_coord_p1(1, :);    % Select the first one
        end

        coord_line_1 = tiles(smallest_coord_p1(1, 1), smallest_coord_p1(1, 2), limit);  % Get the lines between points 1 and 2 translated to tiles
        coord_line_2 = tiles(smallest_coord_p2(1, 1), smallest_coord_p2(1, 2), limit);  % Get the lines between points 1 and 2 translated to tiles

        for n = 1:9 % Draw the lines between points on all tiles
            line([coord_line_1(n, 1), coord_line_2(n, 1)],[coord_line_1(n, 2), coord_line_2(n, 2)], 'Color', [0.5 0.5 0.5]);
        end     

    elseif d1 <= sqrt((0.5*limit)^2 + (0.5*limit)^2) % If the distance between the points is smaller than a circle with a radius of d1
        line([x1, x2], [y1, y2], 'Color', [0.5 0.5 0.5]);   % Draw line between points
    end
end