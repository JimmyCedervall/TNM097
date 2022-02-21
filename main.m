addpath('databas')

% Create database matrix
database = cell(1,40);

%imgTEMP = imresize(imgTEMP,[25,25]);
%imgTEMP = rgb2lab(imgTEMP);
%figure
%imshow(cell2mat(database(i)))

for i = 1:size(database,2)
    fileName = strcat(int2str(i),'.jpg');
    img = im2double(imread(fileName));
    database{i} = img;
end

fileName = 'godgaren.jfif';
img = image_rec_v2(8, im2double(imread(fileName)), database);

imshow(img)