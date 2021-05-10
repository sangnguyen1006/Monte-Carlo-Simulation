%noi suy cong suat cho den LED
function led_interp()

   global lamda power l p lamda_interp
   %doc file.txt
   fileID = fopen('led_data.txt','r');
   value = fscanf(fileID,'%f %f',[2 Inf]);
   fclose(fileID);

   %lay lamda va phan tram nang luong
   l=value(1,1:size(value,2));
   p=value(2,1:size(value,2));

   %thiet lap
   lamda_interp=[l(1) l(length(l))];
   delta_lamda= 1;
   interp_type='spline';
   lamda= lamda_interp(1):delta_lamda:lamda_interp(2);

   %noi suy nang luong cho tung lamda
   [~, ind] = unique(l); %xu ly su trung lap trong matrix power
   power=abs(interp1(l(ind), p(ind), lamda, interp_type));
   leddata_interp=[lamda',power'];

   %luu vao file.mat
   save('led_interp.mat','leddata_interp');

end