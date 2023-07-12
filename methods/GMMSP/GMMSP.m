function GMMSP( imgPath, k, outPath, dtime_file, time_file )

e_c    = 8.000000;
T      = 10;
lambda = 8.000000;
e_s    = 2.000000;

addpath ./src
img_names = dir(imgPath);

t = nan(1, length(img_names));

for i=1:length(img_names)

if img_names(i).isdir == 0

image_name = img_names(i).name;
[~, image_stem, ext] = fileparts(image_name);

% image and v_x,v_y
img = imread([img_names(i).folder, '/', image_stem, ext]);

img_size = length(size(img));
if img_size == 3
    [rows,cols, aa] = size(img);
else
    [rows,cols] = size(img);
    img = repmat(reshape(img, [rows, cols, 1]), [1, 1, 3, 1]);
end

v_x = sqrt((k*cols)/rows);
v_y = (rows/cols)*v_x;

v_x = floor(cols/ceil(v_x));
v_y = floor(rows/ceil(v_y));

if(v_x < 8)
    v_x = 8;
end

if(v_y < 8)
    v_y = 8;
end

% call GMMSP
tic;
label = mx_GMMSP(img, v_x, v_y, e_c, T, lambda, e_s);
time = toc;
t(i) = time;

if strcmp(outPath,'-1') == 0
	imwrite(uint16(label), [outPath, '/', image_stem, '.pgm']);
end

if strcmp(dtime_file,'-1') == 0
	fileID = fopen(dtime_file,'a+');
	fprintf(fileID,'%s %.5f\n',image_stem, time);
	fclose(fileID);
end

end
end

if strcmp(time_file,'-1') == 0
	sum_time = t(~isnan(t));
	fileID = fopen(time_file,'a+');
	fprintf(fileID,'%d %.5f\n',k, mean(sum_time));
	fclose(fileID);
end

end
