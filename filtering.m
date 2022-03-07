function [smallDatabase,meanSmallDB] = filtering(numberOfCluster,clusterSize,tinyImgSize)
%https://www.mathworks.com/content/dam/mathworks/tag-team/Objects/c/88360_93001v00_Color-Based_Seg_K-Means_Clustering_2016.pdf

addpath('databas')
allSize = 7128;
database = cell(1,allSize);
% numberOfCluster = 25;
% clusterSize = 4;

for i = 1:size(database,2)
    fileName = strcat(int2str(i),'.jpg');
    img = im2double(imread(fileName));
    database{i} = img;
end

allChannels = zeros(3,allSize);

for i = 1:allSize
    fileName = strcat(int2str(i),'.jpg');
    img = im2double(imread(fileName));
    allChannels(1,i) = mean(mean(img(:, :, 1)));
    allChannels(2,i) = mean(mean(img(:, :, 2)));
    allChannels(3,i) = mean(mean(img(:, :, 3)));
end
X = allChannels';
[idx,C] = kmeans(X,numberOfCluster);

% plot 3D clusters
scatter3(X(:,1),X(:,2),X(:,3),15,idx,'filled')

% go through all clusters in 'idx' 
cellOfCluster = cell(numberOfCluster, allSize);

% Yttre loopen går igenom alla siffror 1 till allSize
for row = 1:allSize
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
tempDB = cell(1,clusterSize);
for i = 1:numberOfCluster
    tempDB = cellOfCluster(i,(1:clusterSize));
    for j = 1:clusterSize
        if (isempty(tempDB(1,j)))
            cellOfCluster(i:j) = database(1, randi([1 allSize]));
        end
    end
end

for i = 1:numberOfCluster
    smallDatabase(1,(i-1)*clusterSize+1:(i-1)*clusterSize+clusterSize) = cellOfCluster(i,1:clusterSize);
end

meanSmallDB = zeros(3,size(smallDatabase,2));
for i=1:size(smallDatabase,2)
    matSmallImg = rgb2lab(cell2mat(smallDatabase(1,i)));    
    
    meanSmallDB(1,i) = mean(mean(matSmallImg(:,:,1)));
    meanSmallDB(2,i) = mean(mean(matSmallImg(:,:,2)));
    meanSmallDB(3,i) = mean(mean(matSmallImg(:,:,3)));
    
    
    resized = imresize(cell2mat(smallDatabase(1,i)),[tinyImgSize tinyImgSize]);
    smallDatabase(1,i) = mat2cell(resized,tinyImgSize);
end

disp(strcat(strcat('Database has been created containing:', int2str(size(smallDatabase,2))),' images'));

% for i=1:size(smallDatabase,2)
%     imgTEMP = cell2mat(smallDatabase(1,i));
%     nameOfFile=strcat(string(i),".jpg");
%     path = strcat("databaseNew/",nameOfFile);
%     imwrite(imgTEMP,path);    
% end
end

