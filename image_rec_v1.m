function [image] = image_rec_v1(smallCellSize,img)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    %New code that uses database TRIPLE NESTED
    
% rgb till lab
imgLAB = rgb2lab(double(img));
    
% create cells from the original img
imgTiles = mat2tiles(imgLAB, [smallCellSize,smallCellSize]);

image = zeros(size(img,1),size(img,1),3);
siz = smallCellSize;
siz = int16(siz(1));
imgREF = zeros(smallCellSize,smallCellSize,3);

% for-loop for each "small cell"
for i = 1:size(img,1)/smallCellSize
    for j = 1:size(img,1)/smallCellSize
        %Create cell image from original image
        imageMatrix = cell2mat(imgTiles(i,j));
        ref = 0;
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
            % Add small image to final image
            image(i*siz - siz+1 : i*siz, j*siz - siz+1 : j*siz,:) = imgREF;
    end
end

end

