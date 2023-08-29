function define_HA_vectors(psi_range,phi_range,theta_range)
%==========================================================================
% This function models configurations of the second Fab arm for a given
% range of transformation angles (illustrated in Figure 6 of the
% manuscript). This saves a file which can be analyzed using
% 'analyze_HA_vectors.m'
%==========================================================================

load('aligned_Fabs.mat');

for k = 1:length(pdb);
  disp(pdb{k});

  n = 1;
  while n <= 1000;

  pHA = xyz_HA{k}; 
  pHC = xyz_HC{k}; 
  pLC = xyz_LC{k};

  % Sample a range of positions for the second Fab arm.
  % First, rotate about the x-axis, displaced by 10nm:
  D = 100;
  qHA = pHA; qHA(:,2) = -pHA(:,2); qHA(:,1) = -pHA(:,1)+D; 
  qHC = pHC; qHC(:,2) = -pHC(:,2); qHC(:,1) = -pHC(:,1)+D; 
  qLC = pLC; qLC(:,2) = -pLC(:,2); qLC(:,1) = -pLC(:,1)+D;
  C0 = mean([pHC;pLC;qHC;qLC]);

  % Introduce twist and symmetry-breaking wobble:
  d = 10*(rand(1)-1/2);                     % Introduce random translation of +/- ~0.5nm
  qHA(:,1) = qHA(:,1)+d; qHC(:,1) = qHC(:,1)+d; qLC(:,1) = qLC(:,1)+d;
  C = mean([qHC;qLC]);                      % Recenter at Fab2
  [pHA,pHC,pLC] = shift_center(pHA,pHC,pLC,C);
  [qHA,qHC,qLC] = shift_center(qHA,qHC,qLC,C);
  t = [psi_range*randn(1)]*pi/180;          % Twist angle, in degrees
  Rx = [1 0 0 ; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)];
  qHA = rotate_points(qHA,Rx); qHC = rotate_points(qHC,Rx); qLC = rotate_points(qLC,Rx);
  w = [10*randn(1,2)]*pi/180;               % Introduce rotational wobble of +/- ~10 degrees;
  Rz = [cos(w(1)) -sin(w(1)) 0 ; sin(w(1)) cos(w(1)) 0 ; 0 0 1];
  Ry = [cos(w(2)) 0 sin(w(2)) ; 0 1 0 ; -sin(w(2)) 0 cos(w(2))];
  qHA = rotate_points(qHA,Rz); qHC = rotate_points(qHC,Rz); qLC = rotate_points(qLC,Rz);
  qHA = rotate_points(qHA,Ry); qHC = rotate_points(qHC,Ry); qLC = rotate_points(qLC,Ry);

  % Recenter at the original midpoint of the Fabs:
  [pHA,pHC,pLC] = shift_center(pHA,pHC,pLC,C0-C);
  [qHA,qHC,qLC] = shift_center(qHA,qHC,qLC,C0-C);
  
  % Next, rotate about the z-axis:
  a = [phi_range*randn(1)]*pi/180;
  Rz = [cos(a) -sin(a) 0 ; sin(a) cos(a) 0 ; 0 0 1];
  qHA = rotate_points(qHA,Rz);
  qHC = rotate_points(qHC,Rz);
  qLC = rotate_points(qLC,Rz);
  
  % Finally, rotate one of the Fab arms about the y-axis:
  a = [theta_range*randn(1)-60]*pi/180; % 'Resting position' ~ 60 degrees.
  Ry = [cos(a) 0 sin(a) ; 0 1 0 ; -sin(a) 0 cos(a)];
  qHA = rotate_points(qHA,Ry);
  qHC = rotate_points(qHC,Ry);
  qLC = rotate_points(qLC,Ry);

  % Re-center things at the midpoint of the first Fab:
  C = mean([pHC;pLC]);
  [pHA,pHC,pLC] = shift_center(pHA,pHC,pLC,C);
  [qHA,qHC,qLC] = shift_center(qHA,qHC,qLC,C);

  % Check the distance between the two Fabs:
  p = [pHC;pLC]; q = [qHC;qLC];
  dx = p(:,1)*ones(1,length(p)) - ones(length(q),1)*q(:,1)';
  dy = p(:,2)*ones(1,length(p)) - ones(length(q),1)*q(:,2)';
  dz = p(:,3)*ones(1,length(p)) - ones(length(q),1)*q(:,3)';
  dr = sqrt((dx.^2 + dy.^2 + dz.^2));

  % Define a vector for the axis of the two HAs (if Fabs are far enough apart):
  if(min(min(dr)) > 5)
    h1{n} = [pHA(177,:) ; pHA(441,:)]; % For central stalk, use 387 over 177
    h2{n} = [qHA(177,:) ; qHA(441,:)]; % For central stalk, use 387 over 177
    n = n + 1;
  end

  end
  H1{k} = h1;
  H2{k} = h2;

end
savename = ['HA_vectors_psi',num2str(psi_range),'_phi',num2str(phi_range),'_theta',num2str(theta_range),'.mat'];
save(savename,'H1','H2','pdb')
end

%==========================================================================
function y = rotate_points(x,R);
for k = 1:size(x,1);
  y(k,:) = [R*x(k,:)']';
end
end
%==========================================================================
function [x,y,z] = shift_center(X,Y,Z,C);
x(:,1) = X(:,1) - C(1);
x(:,2) = X(:,2) - C(2);
x(:,3) = X(:,3) - C(3);
y(:,1) = Y(:,1) - C(1);
y(:,2) = Y(:,2) - C(2);
y(:,3) = Y(:,3) - C(3);
z(:,1) = Z(:,1) - C(1);
z(:,2) = Z(:,2) - C(2);
z(:,3) = Z(:,3) - C(3);
end
%==========================================================================