% You are free to use, change or redistribute this code for any
% non-commrecial purposes.If you use this software,please cite the
% following in any resulting publication and email us:
% [1] Zhengqin Li, Jiansheng Chen, Superpixel Segmentation using Linear 
%     Spectral Clustering, IEEE Conference on Computer Vision and Pattern 
%     Recognition (CVPR), Jun. 2015 
% (C) Zhengqin Li, Jiansheng Chen, 2014
% li-zq12@mails.tsinghua.edu.cn
% jschenthu@mail.tsinghua.edu.cn
% Tsinghua University

function LSC(input_path, output_path, num_superpix, dtime_file, time_file)
  pkg load image;
  
  img_names = dir(input_path);

	t = nan(1, length(img_names));

	for i=1:length(img_names)

		if img_names(i).isdir == 0

			image_name = img_names(i).name;
			[~, image_stem, ext] = fileparts(image_name);
	
			% image and v_x,v_y
			readimg = imread([img_names(i).folder, '/', image_stem, ext]);
  
		  if(isgray(readimg) == true)
  		  img = uint8(zeros(size(readimg,1), size(readimg,2), 3));
  		  img(:,:,1) = img(:,:,2) = img(:,:,3) = readimg;
  		else
  		  img = readimg;
  		endif
  
	  	gaus=fspecial('gaussian',3);
  	
		  I=imfilter(img,gaus);
  
	  	ratio=0.075;

			tic;
			label=LSC_mex(I,num_superpix,ratio);
			time = toc;
			t(i) = time;
		  
			if strcmp(output_path,'-1') == 0
				imwrite(uint16(label), [output_path, '/', image_stem, '.pgm']);
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
		fprintf(fileID,'%d %.6f\n',num_superpix, mean(sum_time));
		fclose(fileID);
	end

endfunction




