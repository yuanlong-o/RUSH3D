function [WDF] = wdfload_module(frame,realign_datafilename,wdf_size)
%% WDF Load Module Load WDF
% This program is used to load specified frame of WDF
% Last update: 05/16/2021. MW

% Input:
% frame
% realign_datafilename
% wdf_size                    size of WDF (4)

% Output:
% WDF

WDF = zeros(wdf_size);
for u = 1 : Nnum
    for v = 1 : Nnum
        temp = double(imread(strcat(realign_datafilename,num2str(frame),'.tif'),(u-1)*Nnum+v));
        WDF(:,:,u,v) = temp;
        disp(strcat('u = ',num2str(u),', v = ',num2str(v),' has been load...'));
    end
end

end