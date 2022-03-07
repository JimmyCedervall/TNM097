function [imgOUT, imgSSIM, imgSNR, imgDE] = image_rec_v1(fileName, database, meanDatabase, smallCellSize, n)
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
    disp(strcat(('Image dimensions are not divisible by: '), int2str(smallCellSize)));
    disp(strcat('The size of the image has been changed to:', strcat((strcat(int2str(size(img,1)),'x'))),int2str(size(img,2))));
    disp('Some information from the original image has therefore been lost');
end

% rgb till lab
img = rgb2lab(double(img));

% create cells from the original img
imgTiles = mat2tiles(img, [smallCellSize,smallCellSize]);

% imgOUT scale factor.
% Factor of 3 results in an output image of 3 times the size
imgOUTscale = 3;

% final image matrix
imgOUT = zeros(size(img,1)*imgOUTscale,size(img,2)*imgOUTscale,3);

% går att göra dessa två rader finare
siz = int16(smallCellSize*imgOUTscale);
%siz = int16(siz(1));

% create image reference to later be used
imgREF = zeros(smallCellSize*imgOUTscale,smallCellSize*imgOUTscale,3);

timerVal = tic;

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
            %imgTEMP = cell2mat(database(k));
            %imgTEMP = imresize(imgTEMP,[smallCellSize,smallCellSize]);
            %imgTEMP = rgb2lab(imgTEMP);
                     
            % calculating delta E values for all images
            %dE(1,k) = mean(mean(sqrt( (imgTEMP(:,:,1) - imageMatrix(:,:,1)).^2 + (imgTEMP(:,:,2) - imageMatrix(:,:,2)).^2 + (imgTEMP(:,:,3) - imageMatrix(:,:,3)).^2)));

            L = (meanDatabase(1,k) - mean(mean(imageMatrix(:,:,1))))^2;
            A = (meanDatabase(2,k) - mean(mean(imageMatrix(:,:,2))))^2;
            B = (meanDatabase(3,k) - mean(mean(imageMatrix(:,:,3))))^2;
            dE(1,k) = (sqrt(L + A + B));
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
                
        % go through the best matches and find THE best match
        ref = 0;
        for k = 1:size(bestMatches,2)
            %tinyRef = imresize(cell2mat(bestMatches(k)), [smallCellSize*imgOUTscale, smallCellSize*imgOUTscale]);
            tinyRef = cell2mat(bestMatches(k));
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
    % Give quick feedback for how long time is left aproximetly
    check = size(img,1) * size(img,2);
    if check > 5000000
        if round((size(img,1)/smallCellSize) * progress / 100)  == i
            progress = progress + 1;
            disp(strcat('Image is processing...', int2str((progress - 1) * 1), '%'));
        end
    elseif check > 2000000
        if round((size(img,1)/smallCellSize) * progress / 20)  == i
            progress = progress + 1;
            disp(strcat('Image is processing...', int2str((progress - 1) * 5), '%'));
        end
    else
        if round((size(img,1)/smallCellSize) * progress / 10)  == i
            progress = progress + 1;
            disp(strcat('Image is processing...', int2str((progress - 1) * 10), '%'));
        end
    end
    
end

elapsedTime = toc(timerVal);

% Save final image
% if n == 1
%     imwrite(imgOUT, strcat('final_', n, '_', fileName));
% else
%     imwrite(imgOUT, strcat('final_', n, '_', fileName));
% end
imwrite(imgOUT, strcat('final_N', int2str(n), '_A', int2str(size(database,2)), fileName));

%%%%%%%%% Objektiva kvalitetsmått

% SSIM
imgOUT = imresize(imgOUT,[(size(img,1)) (size(img,2))]);
imgSSIM = ssim(double(img), double(imgOUT));

% SNR
imgSNR = mysnr(double(img), double(img) - double(imgOUT));

% DELTA E spacial cielab
%img = rgb2lab(img);
%imgOUT = rgb2lab(imgOUT);
%imgDE = mean(mean(sqrt( (img(:,:,1) - imgOUT(:,:,1)).^2 + (img(:,:,2) - imgOUT(:,:,2)).^2 + (img(:,:,3) - imgOUT(:,:,3)).^2)));

% SCIELAB computation
img = rgb2xyz(img);
imgOUT = rgb2xyz(imgOUT);
whitePoint = [95.05 100 108.9];
sampDegree = 72 * 30 * 0.0175;
imgDE = scielab(sampDegree, img, imgOUT, whitePoint, 'xyz');

% print information
disp(strcat('SSIM:', sprintf('%.6f',imgSSIM)));
disp(strcat('SNR:', sprintf('%.6f',imgSNR)));
disp(strcat('DELTA E:', sprintf('%.6f',imgDE)));
disp(strcat(strcat('Elapsed time:', int2str(elapsedTime)/60), '-min'));

end

