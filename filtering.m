function [smallDatabase,meanSmallDB] = filtering(numberOfCluster,clusterSize,tinyImgSize)
<<<<<<< Updated upstream
%link to database: https://www.kaggle.com/theblackmamba31/landscape-image-colorization

% This function was created by:
% Jimmy Cedervall Lamin (jimla401)
% Edvin Nordin (edvno177)
% Carl Melin (carme007)
% 
%
% * The function takes in three inputs
% * numberOfCluster is the amount of clusters
% * clusterSize is the amount of images taken from each cluster
% * tinyImgSize is how big the resize images should be

% * Three outputs are generated:
% * smallDatabase is the filtered database
% * meanSmallDB is all the mean values of images in smallDatabase
=======
%https://www.mathworks.com/content/dam/mathworks/tag-team/Objects/c/88360_93001v00_Color-Based_Seg_K-Means_Clustering_2016.pdf
>>>>>>> Stashed changes

addpath('databas')
allSize = 7128;
database = cell(1,allSize);

%add to database
for i = 1:size(database,2)
    fileName = strcat(int2str(i),'.jpg');
    img = im2double(imread(fileName));
    database{i} = img;
end
%idx is which cluster each number belongs to
allChannels = zeros(3,allSize);
for i = 1:allSize
    img = cell2mat(database(1,i));
    allChannels(1,i) = mean(mean(img(:, :, 1)));
    allChannels(2,i) = mean(mean(img(:, :, 2)));
    allChannels(3,i) = mean(mean(img(:, :, 3)));
end
X = allChannels';
[idx,C] = kmeans(X,numberOfCluster);

% plot 3D clusters
scatter3(X(:,1),X(:,2),X(:,3),15,idx,'filled')

%hitta hur stor smalldatabase ska vara
most = mode(idx');
amountMost=0;
for i=1:size(idx,1)
    if idx(i,1) == most
        amountMost = amountMost+1;
    end
end
% go through all clusters in 'idx' 
cellOfCluster = cell(numberOfCluster, amountMost);

for row = 1:amountMost
    col = 1;
    % Inre loopen ska gå igenom alla idx värden för att hitta match
    for i = 1:allSize
        % idx(i,1) borde vara "gruppvärdet" som dess 'id' tillhör
        if idx(i,1) == row
            cellOfCluster(row,col) = database(1,i);
            col = col + 1;
        end
    end
end

%from cellOfCluster to smallDatabase
for i = 1:numberOfCluster
    smallDatabase(1,(i-1)*clusterSize+1:(i-1)*clusterSize+clusterSize) = cellOfCluster(i,1:clusterSize);
end
%rezise images in smallDatabase
for i=1:size(smallDatabase,2)
    if size(cell2mat(smallDatabase(1,i)),3) > 2  
        resized = imresize(cell2mat(smallDatabase(1,i)),[tinyImgSize tinyImgSize]);
        smallDatabase(1,i) = mat2cell(resized,tinyImgSize); 
    end
end

%remove the empty cells
smallDatabase = smallDatabase(~cellfun('isempty',smallDatabase));

%create the mean of small database
meanSmallDB = zeros(3,size(smallDatabase,2));
for i=1:size(smallDatabase,2)
    matSmallImg = rgb2lab(cell2mat(smallDatabase(1,i)));   
     
    meanSmallDB(1,i) = mean(mean(matSmallImg(:,:,1)));
    meanSmallDB(2,i) = mean(mean(matSmallImg(:,:,2)));
    meanSmallDB(3,i) = mean(mean(matSmallImg(:,:,3)));
end

%removes empty values
disp(strcat(strcat('Database has been created containing:', int2str(size(smallDatabase,2))),' images'));


%uncomment if you want to se all images in the smalldatabase
% for i=1:size(smallDatabase,2)
%     imgTEMP = cell2mat(smallDatabase(1,i));
%     nameOfFile=strcat(string(i),".jpg");
%     path = strcat("folderName/",nameOfFile);
%     imwrite(imgTEMP,path);    
% end
end

