
clc

fname = 'may12_WORKS.mi';

f = fopen(fname, 'r');

l1 = fgets(f);

T = 1;
n_buffers = 0;
buffers = {};
while(T)
   l = fgets(f);
   if startsWith(l, 'bufferLabel')
%        disp(l) 
       n_buffers = n_buffers + 1;
       buffer_i = strsplit(l);
       buffer_i = buffer_i{2};
       buffers = {buffers, buffer_i};
   end
   
   if startsWith(l, 'xPixels')
       xpix = strsplit(l);
       xpix = str2num(xpix{2});
   end
   if startsWith(l, 'yPixels')
       ypix = strsplit(l);
       ypix = str2num(ypix{2});
   end
   
   
   
   if startsWith(l, 'data')
       T = 0;
   end
    
end


% Get buffer 1 data

data = struct();

datmat = zeros(xpix, ypix);
for jj = 1:ypix
    for kk = 1:xpix
       datmat(jj, kk) = str2num(fgets(f));

    end

end
%%
csvwrite('grating01.csv', datmat)
cmp = contrast(datmat);
HeatMap(datmat, 'Colormap', cmp)


