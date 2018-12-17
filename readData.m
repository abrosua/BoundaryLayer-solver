function [xb,yb,m,mp1,Uinf,alpha] = readData(filename)
  %% this function is to read data from given airfoil coordinate data and also read a user-based input for freestream velocity and angle of attack
  %% 
  % specifying a coordinates (xb,yb) of boundary points on airfoil surface. the
  % last point coincides with the first
  fileID = fopen(filename, 'r');
  formatSpec = '%f';
  sizeA = [2 inf];
  A = fscanf(fileID, formatSpec, sizeA);
  fclose(fileID);
  
  xy = A';

  % making coordinate to be column-wise
  xb0=xy(:,1)'; yb0=xy(:,2)';

  % calculating total panel point 
  m   = length(xb0)-1; % total point with TE & LE is treated as 1 point 
  mp1 = m+1; % total point including TE and LE
  % reverse point index, so the indexing is began from TE and numbered from the lower surface of airfoil
  for i=1:mp1
      xb(i)=xb0(mp1+1-i);
      yb(i)=yb0(mp1+1-i);
  end

  fprintf('input parameter:\n');
  % setting freestream velocity through user input by a control condition 
  Uinf = input('freestream velocity (1-100m/s) : ');
  while (Uinf<=0) || (Uinf>100)
    Uinf = input('incorrect input, freestream velocity (1-100m/s) : ');
  end
  % setting up the angle of attack(AoA) in [radian] through user input by a control condition 
  alphaDeg = input('angle of attack (-5 < AoA < 10deg) : '); % AoA in [degree]
  while (alphaDeg>10) || (alphaDeg<-4)
    alphaDeg = input('incorrect input, angle of attack (-5 < AoA < 10deg) : '); % AoA in [degree]
  end
  alpha = alphaDeg*pi/180; % AoA in [radian]

end
