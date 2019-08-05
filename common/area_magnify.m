function Y = area_magnify(X, areaStartEnd, times, pos)
%  X: original imgae, matrix or 3-dimensional tensor
%  areaStartEnd: 2 x 2 matrix, the coordinate of the area needed to be magnified.
%                Left up corner is (0,0)
%                [x_A,y_A;
%                 x_B,y_B]
%  times: times of magnify, scalar
%  pos: the position to put the magnify area; 1:left up corner; 2:left down vcorner
%                                             3:right up corner; 4: right down corner
%

X = std_image(X);
sizeX = size(X);
indX  = areaStartEnd(1,1):areaStartEnd(2,1);
indY  = areaStartEnd(1,2):areaStartEnd(2,2);

if length(sizeX) == 3
    Y = X;
else
    Y = repmat(X, [1, 1, 3]);
end

YOri = Y;

Y(areaStartEnd(1,2):areaStartEnd(1,2)+1, indX, :) = cat(3, ones(2, length(indX)), zeros(2, length(indX)), zeros(2, length(indX)));
Y(areaStartEnd(2,2):areaStartEnd(2,2)+1, indX, :) = cat(3, ones(2, length(indX)), zeros(2, length(indX)), zeros(2, length(indX)));
Y(indY, areaStartEnd(1,1):areaStartEnd(1,1)+1, :) = cat(3, ones(length(indY), 2), zeros(length(indY), 2), zeros(length(indY), 2));
Y(indY, areaStartEnd(2,1):areaStartEnd(2,1)+1, :) = cat(3, ones(length(indY), 2), zeros(length(indY), 2), zeros(length(indY), 2));


areaNeedMagnify = YOri(indY, indX, :);
bigArea         = imresize(areaNeedMagnify, times);


[bigAreaH, bigAreaW, ~] = size(bigArea);
bigArea(1:2, :, :) = cat(3, ones(2, bigAreaW), zeros(2, bigAreaW), zeros(2, bigAreaW));
bigArea(end-1:end, :, :) = cat(3, ones(2, bigAreaW), zeros(2, bigAreaW), zeros(2, bigAreaW));
bigArea(:, 1:2, :) = cat(3, ones(bigAreaH, 2), zeros(bigAreaH, 2), zeros(bigAreaH, 2));
bigArea(:, end-1:end, :) = cat(3, ones(bigAreaH, 2), zeros(bigAreaH, 2), zeros(bigAreaH, 2));

switch pos
    case 1
        indBigH = 3;
        indBigW = 3;
    case 2
        indBigH = (sizeX(1)-3) - bigAreaH + 1;
        indBigW = 3;
    case 3
        indBigH = 3;
        indBigW = (sizeX(2)-3) - bigAreaW + 1;
    case 4
        indBigH = (sizeX(1)-3) - bigAreaH + 1;
        indBigW = (sizeX(2)-3) - bigAreaW + 1;
end

Y(indBigH:(indBigH+bigAreaH-1), indBigW:(indBigW+bigAreaW-1), :) = bigArea;

end
