function ERS(input_path, output_path, num_superpixels, dtime_file, time_file)
  pkg load image;
  
lambda_prime = 0.5;
sigma = 5.0;
conn8 = 1;
  
  img_names = dir(input_path);

	t = nan(1, length(img_names));

	for i=1:length(img_names)

		if img_names(i).isdir == 0

			image_name = img_names(i).name;
			[~, image_stem, ext] = fileparts(image_name);

		  img = imread([img_names(i).folder, '/', image_stem, ext]);
  		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		tic;
  		[labels] = mex_ers(double(img),num_superpixels);
  		time = toc;
		t(i) = time;
   
   		if strcmp(output_path,'-1') == 0
	  		imwrite(uint16(labels + 1), [output_path, '/', image_stem, '.pgm']);
  		endif
  		
  		if strcmp(dtime_file,'-1') == 0
				fileID = fopen(dtime_file,'a+');
				fprintf(fileID,'%s %.6f\n',image_stem, time);
				fclose(fileID);
			endif

		endif
	endfor

	if strcmp(time_file,'-1') == 0
		sum_time = t(~isnan(t));
		fileID = fopen(time_file,'a+');
		fprintf(fileID,'%d %.6f\n',num_superpixels, mean(sum_time));
		fclose(fileID);
	endif
	
endfunction
