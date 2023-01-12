function imgResult = Filter(img, kernel)
    imgResult = img;
    sizeImg = size(img);
    sizeKernel = size(kernel);
    tmpImg = zeros(sizeImg+fix(sizeKernel/2)*2);
    tmpImg(fix(sizeKernel(1)/2)+1:fix(sizeKernel(1)/2)+sizeImg(1), fix(sizeKernel(2)/2)+1:fix(sizeKernel(2)/2)+sizeImg(2)) = img;
    for row = 1:sizeImg(1)
        for col = 1:sizeImg(2)
            imgResult(row, col) = sum(tmpImg(row:row+sizeKernel(1)-1, col:col+sizeKernel(2)-1) .* kernel, "all");
        end
    end
end