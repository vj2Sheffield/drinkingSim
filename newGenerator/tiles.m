% This function accepts two points as well as a map size, lim. It then
% produces a 3-by-3 tile surface with a centre on the two points, offset by
% lim, and will then return the tile surface f.
% Victoria Johnson
% December 2019

function f = tiles(x, y, lim)
    f = [x - lim, y + lim;
            x, y + lim;
            x + lim, y + lim;
            x - lim, y;
            x, y;
            x + lim, y;
            x - lim, y - lim;
            x, y - lim;
            x + lim, y - lim];
end