function out = Canny(in, threshold1, threshold2, L2gradient)
   %Gaussian Blur
   %Processing====================================================
   StandardDeviation_GaussianBlur = 0.707;
   KernelSize_GaussianBlur = 3;
   [gridX, gridY] = meshgrid(-(KernelSize_GaussianBlur-1)/2:(KernelSize_GaussianBlur-1)/2,-(KernelSize_GaussianBlur-1)/2:(KernelSize_GaussianBlur-1)/2);
   Kernel_GaussianBlur = (1/(2*pi*StandardDeviation_GaussianBlur^2))*exp(1).^(-((gridX.^2+gridY.^2)/(2*StandardDeviation_GaussianBlur^2)));
   Kernel_GaussianBlur = Kernel_GaussianBlur / sum(Kernel_GaussianBlur,"all");
   in_After_GaussianBlur = Filter(in, Kernel_GaussianBlur);
   
   %Sobel Gradient Direction====================================================
   KernelSize_Sobel = 3;
   coefficients = [2 20 780 132600];
   [gridX, gridY] = meshgrid(-(KernelSize_Sobel-1)/2:(KernelSize_Sobel-1)/2,-(KernelSize_Sobel-1)/2:(KernelSize_Sobel-1)/2);
   
   %kernel for Vertical
   Kernel_Sobel_Vertical = gridY ./ (gridX .* gridX + gridY .* gridY); 
   Kernel_Sobel_Vertical((KernelSize_Sobel+1)/2, (KernelSize_Sobel+1)/2) = 0;
   Kernel_Sobel_Vertical = Kernel_Sobel_Vertical * coefficients((KernelSize_Sobel-1)/2);
   
   %kernel for Horizon
   Kernel_Sobel_Horizon = gridX ./ (gridX .* gridX + gridY .* gridY);
   Kernel_Sobel_Horizon((KernelSize_Sobel+1)/2, (KernelSize_Sobel+1)/2) = 0;
   Kernel_Sobel_Horizon = Kernel_Sobel_Horizon * coefficients((KernelSize_Sobel-1)/2);

   Gx = double(Filter(in_After_GaussianBlur, Kernel_Sobel_Horizon));
   Gy = double(Filter(in_After_GaussianBlur, Kernel_Sobel_Vertical));
   
   TAN22_5 = 0.414;
   TAN67_5 = 2.414;
   VERTICAL = 0;
   HORIZON = 1;
   SLASH = 2;
   BACK_SLASH = 3;
   ZERO = 4;
   
   eta = 10e-6;
   theta = ((Gy./(Gx + eta)));
   thetaDirection = zeros(size(theta));
   thetaDirection(theta < TAN67_5 & theta >= TAN22_5) = SLASH;
   thetaDirection(theta < TAN22_5 & theta >= -TAN22_5) = HORIZON;
   thetaDirection(theta < -TAN22_5 & theta >= -TAN67_5) = BACK_SLASH;

   G = (Gx.*Gx + Gy.* Gy).^0.5;
   
   %NMS====================================================
   [sDy, sDx] = size(thetaDirection);
   thetaDirection(1,:) = ZERO;
   thetaDirection(sDy,:) = ZERO;
   thetaDirection(:, 1) = ZERO;
   thetaDirection(:, sDx) = ZERO;
   
   %VERTICAL
   direction = thetaDirection==VERTICAL;
   direction = direction(2:sDy-1, 2:sDx-1);
   G_NM = true(size(G));
   pixel = G(2:sDy-1, 2:sDx-1);
   v1 = G(1:sDy-2, 2:sDx-1);
   v2 = G(3:sDy, 2:sDx-1);
   G_NM(2:sDy-1, 2:sDx-1) = direction & (pixel <= v1 | pixel <= v2);
   G(G_NM) = 0;
   
   %HORIZON
   direction = thetaDirection==HORIZON;
   direction = direction(2:sDy-1, 2:sDx-1);
   G_NM = true(size(G));
   pixel = G(2:sDy-1, 2:sDx-1);
   v1 = G(2:sDy-1, 1:sDx-2);
   v2 = G(2:sDy-1, 3:sDx);
   G_NM(2:sDy-1, 2:sDx-1) = direction & (pixel <= v1 | pixel <= v2);
   G(G_NM) = 0;

   %SLASH
   direction = thetaDirection==SLASH;
   direction = direction(2:sDy-1, 2:sDx-1);
   G_NM = true(size(G));
   pixel = G(2:sDy-1, 2:sDx-1);
   v1 = G(3:sDy, 3:sDx);
   v2 = G(1:sDy-2, 1:sDx-2);
   G_NM(2:sDy-1, 2:sDx-1) = direction & (pixel <= v1 | pixel <= v2);
   G(G_NM) = 0;

   %BACK_SLASH
   direction = thetaDirection==BACK_SLASH;
   direction = direction(2:sDy-1, 2:sDx-1);
   G_NM = true(size(G));
   pixel = G(2:sDy-1, 2:sDx-1);
   v1 = G(1:sDy-2, 3:sDx);
   v2 = G(3:sDy, 1:sDx-2);
   G_NM(2:sDy-1, 2:sDx-1) = direction & (pixel <= v1 | pixel <= v2);
   G(G_NM) = 0;

   %ZERO
   direction = thetaDirection==ZERO;
   direction = direction(2:sDy-1, 2:sDx-1);
   G_NM = true(size(G));
   G_NM(2:sDy-1, 2:sDx-1) = direction;
   G(G_NM) = 0;

   %Double Thresholding====================================================
   StrongEdge = threshold2 <= G;

   %Hysteresis
   while true
       WeakEdge = threshold1 <= G & threshold2 > G & (StrongEdge == 0);
       Cnt_Hysteresis = 0;
    
       %VERTICAL
       direction = thetaDirection==VERTICAL;
       direction = direction(2:sDy-1, 2:sDx-1);
       StrongEdge_Hysteresis = false(size(G));
       p1 = WeakEdge(3:sDy, 2:sDx-1);
       p2 = WeakEdge(1:sDy-2, 2:sDx-1);
       StrongEdge_Hysteresis(3:sDy, 2:sDx-1) = StrongEdge_Hysteresis(3:sDy, 2:sDx-1) | (direction & p1);
       StrongEdge_Hysteresis(1:sDy-2, 2:sDx-1) = StrongEdge_Hysteresis(1:sDy-2, 2:sDx-1) | (direction & p2);
       StrongEdge(StrongEdge_Hysteresis) = 1;
       Cnt_Hysteresis = Cnt_Hysteresis + sum(StrongEdge_Hysteresis,"all");
%        figure,imshow(StrongEdge);
       
       %HORIZON
       direction = thetaDirection==HORIZON;
       direction = direction(2:sDy-1, 2:sDx-1);
       StrongEdge_Hysteresis = false(size(G));
       p1 = WeakEdge(2:sDy-1, 3:sDx);
       p2 = WeakEdge(2:sDy-1, 1:sDx-2);
       StrongEdge_Hysteresis(2:sDy-1, 3:sDx) = StrongEdge_Hysteresis(2:sDy-1, 3:sDx) | (direction & p1);
       StrongEdge_Hysteresis(2:sDy-1, 1:sDx-2) = StrongEdge_Hysteresis(2:sDy-1, 1:sDx-2) | (direction & p2);
       StrongEdge(StrongEdge_Hysteresis) = 1;
       Cnt_Hysteresis = Cnt_Hysteresis + sum(StrongEdge_Hysteresis,"all");
%        figure,imshow(StrongEdge);
       
       %SLASH
       direction = thetaDirection==SLASH;
       direction = direction(2:sDy-1, 2:sDx-1);
       StrongEdge_Hysteresis = false(size(G));
       p1 = WeakEdge(1:sDy-2, 3:sDx);
       p2 = WeakEdge(3:sDy, 1:sDx-2);
       StrongEdge_Hysteresis(1:sDy-2, 3:sDx) = StrongEdge_Hysteresis(1:sDy-2, 3:sDx) | (direction & p1);
       StrongEdge_Hysteresis(3:sDy, 1:sDx-2) = StrongEdge_Hysteresis(3:sDy, 1:sDx-2) | (direction & p2);
       StrongEdge(StrongEdge_Hysteresis) = 1;
       Cnt_Hysteresis = Cnt_Hysteresis + sum(StrongEdge_Hysteresis,"all");
%        figure,imshow(StrongEdge);
    
       %BACK_SLASH
       direction = thetaDirection==BACK_SLASH;
       direction = direction(2:sDy-1, 2:sDx-1);
       StrongEdge_Hysteresis = false(size(G));
       p1 = WeakEdge(3:sDy, 3:sDx);
       p2 = WeakEdge(1:sDy-2, 1:sDx-2);
       StrongEdge_Hysteresis(3:sDy, 3:sDx) = StrongEdge_Hysteresis(3:sDy, 3:sDx) | (direction & p1);
       StrongEdge_Hysteresis(1:sDy-2, 1:sDx-2) = StrongEdge_Hysteresis(1:sDy-2, 1:sDx-2) | (direction & p2);
       StrongEdge(StrongEdge_Hysteresis) = 1;
       Cnt_Hysteresis = Cnt_Hysteresis + sum(StrongEdge_Hysteresis,"all");
%        figure,imshow(StrongEdge);
       
       if Cnt_Hysteresis == 0
           break;
       end
   end

   out = StrongEdge;

end