function look()
  global dz dr 
  global H1 F
  global lamda power l p lamda_interp
  
  figure (1)
   plot(l, p ,'o', lamda, power, 'r');
   xlabel('wavelength(nm)');
   ylabel('Relative Spectral Powver');
   title('Typical relative spectral power vs wavelength');
   set(gca,'fontsize',12)%kich thuoc phong chu
   xlim ([lamda_interp(1) lamda_interp(2)]);

  figure(2)
   imagesc(log10(H1(1:150,1:150)));
   title('Su phan bo mat do nang luong, J/mm3');
   xlabel(['do sau x ' num2str(dz) 'mm' ]);
   ylabel(['ban kinh x ' num2str(dr) 'mm']);
   set(gca,'fontsize',12)%kich thuoc phong chu
   colorbar
   colormap(makec2f)% <---------------makec2f.m
   set(colorbar,'fontsize',12)%phong chu colorbar

 figure(3)
   imagesc(log10(F(1:150,1:150)));
   title('Su phan bo fluence rate, W/mm2');
   xlabel(['do sau x ' num2str(dz) 'mm' ]);
   ylabel(['ban kinh x ' num2str(dr) 'mm']);
   set(gca,'fontsize',12)
   colorbar
   colormap(makec2f) % <---------------makec2f.m
   set(colorbar,'fontsize',12)
end