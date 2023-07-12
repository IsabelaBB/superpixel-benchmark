function DRW_demo( imgPath, k, outPath, dtime_file, time_file, beta, sigma )

if nargin < 6
	beta = 1e-4;
end

if nargin < 7
	sigma = 50;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_names = dir(imgPath);

t = nan(1, length(img_names));

	for i=1:length(img_names)

		if img_names(i).isdir == 0

			image_name = img_names(i).name;
			[~, image_stem, ext] = fileparts(image_name);

			% image and v_x,v_y
			img = imread([img_names(i).folder, '/', image_stem, ext]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			img_size = length(size(img));
			if img_size == 3
			    [rows,cols, aa] = size(img);
			else
			    [rows,cols] = size(img);
			    img = repmat(reshape(img, [rows, cols, 1]), [1, 1, 3, 1]);
			end

			tic;
			[label_drw, seed_x, seed_y] = mexDRW(img, k, beta, sigma);
			time = toc;
			t(i) = time;

			if strcmp(outPath,'-1') == 0
				label_map = zeros(k);
				max=0;
				
				% count superpixels size
				for i=1:rows
				    for j=1:cols
    					if label_drw(i,j) > max
    						max = label_drw(i,j);
    	    			end
    	    			label_map(label_drw(i,j)) = label_map(label_drw(i,j)) + 1;
				    end
				end

				% create a mapping to superpixels > 0
				i = 1;
				j = 1;
				while i <= max && j <= max
				    while j <= max && j <= k && label_map(j) == 0
				        j = j+1;
				    end
				    label_map(j) = i;
				    i = i+1;
    				j = j+1;
				end

				% relabel superpixels
				for i=1:rows
				    for j=1:cols
				        label_drw(i,j) = label_map(label_drw(i,j));
    				end
				end
	
				imwrite(uint16(label_drw), [outPath, '/', image_stem, '.pgm']);
			end

			if strcmp(dtime_file,'-1') == 0
				fileID = fopen(dtime_file,'a+');
				fprintf(fileID,'%s %.6f\n',image_stem, time);
				fclose(fileID);
			end
		end
	end

	if strcmp(time_file,'-1') == 0
		sum_time = t(~isnan(t));
		fileID = fopen(time_file,'a+');
		fprintf(fileID,'%d %.6f\n',k, mean(sum_time));
		fclose(fileID);
	end
end

