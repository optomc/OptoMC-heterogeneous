for i = 1:length(face)
    f = face(i,1:3);
    cord = node(f,1:3);

    area(i) = triangleArea3d(cord(1,:), cord(2,:), cord(3,:));
    length(i, 1) = pdist([cord(1,:); cord(2,:)], 'euclidean');
    length(i, 2) = pdist([cord(2,:); cord(3,:)], 'euclidean');
    length(i, 3) = pdist([cord(3,:); cord(1,:)], 'euclidean');    
end
lengthList = length(:);

sprintf('Area - mean: %f, max: %f, min: %f, std: %f', mean(area), max(area), min(area), std(area))
sprintf('Length - mean: %f, max: %f, min: %f, std: %f', mean(lengthList), max(lengthList), min(lengthList), std(lengthList))
