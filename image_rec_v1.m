function [imgOUT, imgSSIM, imgSNR, imgDE, imgSCIELAB] = image_rec_v1(fileName, database, meanDatabase, smallCellSize, n, tinyImgSize)
%
% This function was created by:
% Jimmy Cedervall Lamin (jimla401)
% Edvin Nordin (edvno177)
% Carl Melin (carme007)
% 
%
% * The function takes in an image (must be an image with 3 RGB channels)
% * The size of the image does not matter

% * Four outputs are generated:
% * imgOUT is the recreated image using small images from a database
% * imgSSIM is the SSIM value when comparing imgOUT with img
% * imgSNR is the SNR value of the image compared to input img
% * imgDE is the delta E value between imgOUT and img
% * imgSCIELAB is the S-CIELAB delta E value

% Make image double
img = im2double(imread(fileName));

boolean = 0;

% Check dimensions of matrix and resize if needed to make it divideable by smallCellSize
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
imgOUTscale = tinyImgSize/smallCellSize;

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
            % Compute delta E with Lab values
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
            value = ssim(imresize(rgb2gray(lab2rgb(imageMatrix)), imgOUTscale), rgb2gray(tinyRef));
                         
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
imwrite(imgOUT, strcat('final_N', int2str(n), '_A', int2str(size(database,2)), '_C_', int2str(smallCellSize), fileName));

%%%%%%%%% Objektiva kvalitetsmått
imgOUTtemp = imgOUT;
% SSIM
imgOUTtemp = imresize(imgOUTtemp,[(size(img,1)) (size(img,2))]);
imgSSIM = ssim(rgb2gray(lab2rgb(img)), rgb2gray(imgOUTtemp));

% SNR
imgSNR = mysnr(rgb2gray(lab2rgb(img)), rgb2gray(lab2rgb(img)) - rgb2gray(imgOUTtemp));

% DELTA E
imgOUTde = rgb2lab(imgOUTtemp);
imgDE = sqrt((img(:,:,1) - imgOUTde(:,:,1)).^2 + (img(:,:,2) - imgOUTde(:,:,2)).^2 + (img(:,:,3) - imgOUTde(:,:,3)).^2);
imgDE = mean(mean(imgDE));

% SCIELAB computation
imgXYZ = rgb2xyz(lab2rgb(img));
imgOUTxyz = rgb2xyz(imgOUTtemp);
whitePoint = [95.05 100 108.9];
sampDegree = 109 * 50 * 0.0175;
% computed value for our setup = 38
SCIELAB = scielab(38, imgXYZ, imgOUTxyz, whitePoint, 'xyz');
imgSCIELAB = mean(mean(SCIELAB));

% print information
disp(strcat('SSIM:', sprintf('%.4f',imgSSIM)));
disp(strcat('SNR:', sprintf('%.4f',imgSNR)));
disp(strcat('DELTA E:', sprintf('%.4f',imgDE)));
disp(strcat('DELTA E (scielab):', sprintf('%.4f',imgSCIELAB)));

if elapsedTime < 60
    disp(strcat(strcat('Elapsed time:', int2str(elapsedTime)), '-sec'));
else
    disp(strcat(strcat('Elapsed time:', int2str(elapsedTime/60)), '-min'));
end

end

