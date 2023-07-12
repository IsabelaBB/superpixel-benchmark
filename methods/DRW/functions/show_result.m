function [ result ] = show_result( img, label )
%SHOW_RESULT 此处显示有关此函数的摘要
%   此处显示详细说明
%显示
color = [255,255,0];
result = img;
[rows,cols, ~] = size(img); 
for i=1:rows
    for j=1:cols
        b=0;
        if ((i-1)>=1)&&((j-1)>=1)
            if label(i-1,j-1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
            if label(i-1,j)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
            if label(i,j-1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end   
        elseif ((i-1)>=1)&&((j+1)<=cols)
            if label(i-1,j+1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
            if label(i,j+1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
        elseif ((i+1)<=rows)&&((j-1)>=1)
            if label(i+1,j-1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
            if label(i+1,j)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
        elseif ((i+1)<=rows)&&((j+1)<=cols)
            if label(i+1,j+1)~=label(i,j)
                result(i,j,1)=color(1);
                result(i,j,2)=color(2);
                result(i,j,3)=color(3);
                grid(i, j)=1;
                b=1;
            end
        end
    end
end

end

