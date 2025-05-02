function meanScatter = meanDistNeighbor(x, y)

    % Ensure x_values and y_values are sorted by x_values
    [x, idx] = sort(x);
    y = y(idx);
    
    distance = [];
    pointAngles = [];
    for i = 1:length(x)-1

        dx = x(i+1) - x(i);
        dy = y(i+1) - y(i);

        pointAngles = [pointAngles, atan2(dy, dx)/pi]; % in multiples of pi because I'm physicist not an engineer
        distance = [distance, sqrt(dx^2 + dy^2)];
    end

    deltaAngles = [];
    for i = 1:length(pointAngles)-1
        deltaAngles = [deltaAngles, pointAngles(i+1) - pointAngles(i)]; % angles between the points
        deltaAngles = abs(deltaAngles); % normalize by (pi/2) since largest angle is pi/2; this will give us a value between 0 and 1
    end
    meanDeltaAngle = mean(deltaAngles)*2; % mean angle between points
    meanDist = mean(distance); 
    % meanScatter = meanDist * meanDeltaAngle; % mean distance between points
    meanScatter = meanDeltaAngle;
end