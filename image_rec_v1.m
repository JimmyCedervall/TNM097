function [imgOUT, imgSSIM, imgSNR, imgDE] = image_rec_v1(fileName, database, smallCellSize, n)
%
% This function was created by:
% Jimmy Cedervall Lamin (jimla401)
% Edvin Nordin (edvno177)
% Carl Melin (carme007)
% 
%
% * The function takes in an image (must be an image with 3 RGB channels)
% * The size of the image does not matter

% * Three outputs are generated:
% * imgOUT is the recreated image using small images from a database
% * imgSSIM is the SSIM value when comparing imgOUT with img
% * imgDE is the delta E value between imgOUT and img

% Make image double
img = im2double(imread(fileName));

% the size of the sectioning of the image, prefferable of size 8 but is
% interchangable
%smallCellSize = 128;

% number of optimal images from deltaE
%n = 2;

boolean = 0;

% Check dimensions of matrix and resize if needed to make it divideable by
% '8' or smallCellSize
if mod(size(img,1), smallCellSize) ~= 0
    img = imresize(img, [smallCellSize * ceil(size(img,1) / smallCellSize), size(img,2)]);
    boolean = 1;
end

if mod(size(img,2), smallCellSize) ~= 0
    img = imresize(img, [size(img,1), smallCellSize * ceil(size(img,2) / smallCellSize)]);
    boolean = 1;
end

if boolean == 1
    disp(strcat(('Image dimensions are not divisible by:'), int2str(smallCellSize)));
    disp(strcat('The size of the image has been changed to:', strcat((strcat(int2str(size(img,1)),'x'))),int2str(size(img,2))));
    disp('Some information from the original image has therefore been lost');
end

% rgb till lab
img = rgb2lab(double(img));

% create cells from the original img
imgTiles = mat2tiles(img, [smallCellSize,smallCellSize]);

% imgOUT scale factor
imgOUTscale = 3;

% final image matrix
imgOUT = zeros(size(img,1)*imgOUTscale,size(img,2)*imgOUTscale,3);

% går att göra dessa två rader finare
siz = int16(smallCellSize*imgOUTscale);
%siz = int16(siz(1));

% create image reference to later be used
imgREF = zeros(smallCellSize,smallCellSize,3);

tic;

progress = 1;
% for-loop for each "small cell" to create the imgOUT of small images
for i = 1:size(img,1)/smallCellSize
    for j = 1:size(img,2)/smallCellSize
        
        %Create cell image from original image
        imageMatrix = cell2mat(imgTiles(i,j));
        
        % Matrix for holding all deltaE values
        dE = zeros(1,size(database,2));
        for k = 1:size(database,2)
            % Temporary image created from database
            imgTEMP = cell2mat(database(k));
            imgTEMP = imresize(imgTEMP,[smallCellSize,smallCellSize]);
            imgTEMP = rgb2lab(imgTEMP);
            
            % calculating delta E values for all images
            dE(1,k) = mean(mean(sqrt( (imgTEMP(:,:,1) - imageMatrix(:,:,1)).^2 + (imgTEMP(:,:,2) - imageMatrix(:,:,2)).^2 + (imgTEMP(:,:,3) - imageMatrix(:,:,3)).^2)));
        end
                      
        % array for the 'n' best values
        val = zeros(n,1);
        dEtemp = dE;
        for k=1:n
          [val(k),idx] = min(dEtemp);
          % remove for the next iteration the last smallest value: (maybe not needed)
          dEtemp(idx) = [];
        end
        
        bestMatches = cell(1,n);
        % find the positions of the images containing the best values
        for k=1:n
            [~, col] = find(abs(dE-val(k))<1e-3);

            for h = 1:size(col,2)
                bestMatches(k+h-1) = database(col(1,h));    
            end
        end
                
        % go through the best matches and find the ONE best match
        ref = 0;
        for k = 1:size(bestMatches,2)
            tinyRef = imresize(cell2mat(bestMatches(k)), [smallCellSize*imgOUTscale, smallCellSize*imgOUTscale]);
            value = ssim(imresize(lab2rgb(imageMatrix), imgOUTscale), tinyRef);
                         
            % compare each match to find the best one
            if(value > ref)
                ref = value;
                imgREF = tinyRef;
            end
        end

        % Add to image
        imgOUT(i*siz - siz+1 : i*siz, j*siz - siz+1 : j*siz,:) = imgREF;
        
    end
    % Give quick feedback for how long time is left
    if (size(img,1)/smallCellSize) * progress / 10  == i
        progress = progress + 1;
        disp(strcat('Image is processing...', int2str((progress - 1) * 10), '%'));
    end
end

toc;

% Save final image
if n == 1
    imwrite(imgOUT, strcat('recreated_noSSIM_', fileName));
else
    imwrite(imgOUT, strcat('recreated_', fileName));
end

% Objektiva kvalitetsmått

% SSIM
imgOUT = imresize(imgOUT,[(size(img,1)),(size(img,2))]);
imgSSIM = ssim(img, imgOUT);

% SNR
imgSNR = mysnr(img, img - imgOUT);

% DELTA E
img = rgb2lab(img);
imgOUT = rgb2lab(imgOUT);
imgDE = mean(mean(sqrt( (imgOUT(:,:,1) - img(:,:,1)).^2 + (imgOUT(:,:,2) - img(:,:,2)).^2 + (imgOUT(:,:,3) - img(:,:,3)).^2)));

% Elapsed time
elapsedTime = toc;
elapsedTime = elapsedTime / 60;

% print information
disp(strcat('SSIM:', sprintf('%.6f',imgSSIM)));
disp(strcat('SNR:', sprintf('%.6f',imgSNR)));
disp(strcat('DELTA E:', sprintf('%.6f',imgDE)));
disp(strcat('Elapsed time', elapsedTime));

end

