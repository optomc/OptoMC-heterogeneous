function t = rayPlaneIntersects(plane, photonPos, photonDir)
    % compute ray/plane intersection test
    dotProduct = dot(photonDir, plane.Normal);
    if ((dotProduct < 1E-8) && (dotProduct > -1E-8))
        t = 0;
        return;    
    end
    t = dot(plane.Normal, (plane.Pos - photonPos)) / dotProduct;  
    if (t < -1E-8)
        t = 0;
        return;
    end
end