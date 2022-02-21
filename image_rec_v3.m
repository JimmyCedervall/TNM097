function [image] = image_rec_v3(smallCellSize,img,database)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    % rgb till lab
    imgLAB = rgb2lab(double(img));

    % create cells from the original img
    imgTiles = mat2tiles(imgLAB, [smallCellSize,smallCellSize]);



    % convert matrix to 1 dimension
    imgTiles = imgTiles(:);

    % create a image
    image = zeros(size(img,1),size(img,1),3);
    % turn into 1 dimensional image
    %image = image(:);
    siz = smallCellSize;
    siz = int16(siz(1));
    % create ref image
    imgREF = zeros(smallCellSize,smallCellSize,3);
    for i = 1:size(imgTiles,1)
        % create and reset the value of ref
        ref = 0;
        imageMatrix = cell2mat(imgTiles(i));
        for k = 1:size(database,2)
            % Temporary image created from database
            imgTEMP = cell2mat(database(k));
            imgTEMP = imresize(imgTEMP,[smallCellSize,smallCellSize]);
            imgTEMP = rgb2lab(imgTEMP);

            % Compare TEMP image to "cell" from original image
            [value, ~] = ssim(imageMatrix, imgTEMP);

            if(value > ref)
                ref = value;
                imgREF = imgTEMP;
            end
        end

        image(i*siz - siz+1:i*siz,i*siz - siz+1:i*siz,:) = imgREF;
    end
    image = lab2rgb(image);
end

