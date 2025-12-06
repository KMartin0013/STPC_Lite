function [in, on]=check_polygon_in(XY,lonlon,latlat)

    if sum(isnan(XY(:,1)))==0
        [in,on] = inpolygon(lonlon,latlat,XY(:,1),XY(:,2));
    elseif sum(isnan(XY(:,1)))==1
        idx = find(isnan(XY(:,1)));
        % inpolygon need the outer loop in a counterclockwise direction
        % while inner loop in a clockwise direction.
        XY_poly=[XY(idx-1:-1:1,:); NaN NaN; XY(end:-1:idx,:)];

        [in,on] = inpolygon(lonlon,latlat,XY_poly(:,1),XY_poly(:,2));
    else
        error('There should be NO or ONE nan value in XY. Please check your code defining the region.')
    end

end