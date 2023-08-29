function analyze_HA_vectors()
%==========================================================================
% This function analyzes the output from 'define_HA_vectors.m' and assigns
% cis and trans crosslinking scores for each antibody. For each antibody,
% the function displays the two HA vectors and the histogram of vector 
% products between them. 
%==========================================================================

loadname = 'HA_vectors_psi60_phi30_theta30_manuscript.mat';
load(loadname);
p = 4;  % The exponent 'p' is used to calculate crosslinking scores. 
        % Higher values place greater weight on HA orientations that 
        % are parallel or antiparallel.

for k = 1:length(pdb);
  PDB(k,:) = pdb{k};
end

N = length(H1{1}); % Number of sampled configurations for each antibody

for k = 1:length(H1);
  az = 0; ay = 0;
  for n = 1:N
    X = H1{k}{n}; dX = X(2,:)-X(1,:);
    Y = H2{k}{n}; dY = Y(2,:)-Y(1,:);

    % Rotate the HA vectors so that H1 is aligned with the z-axis:
    C = X(2,:); 
    X(1,:) = X(1,:) - C; X(2,:) = X(2,:) - C;
    Y(1,:) = Y(1,:) - C; Y(2,:) = Y(2,:) - C;
    
    if(n==1);
    a = 0.1*pi/180;
    while abs(dX(2))>0.1
      Rz = [cos(a) -sin(a) 0 ; sin(a) cos(a) 0 ; 0 0 1];
      X = rotate_points(X,Rz); dX = X(2,:)-X(1,:);
      Y = rotate_points(Y,Rz); dY = Y(2,:)-Y(1,:);
      az = az + a;
    end

    while abs(dX(1))>0.1
      Ry = [cos(a) 0 sin(a) ; 0 1 0 ; -sin(a) 0 cos(a)];
      X = rotate_points(X,Ry); dX = X(2,:)-X(1,:);
      Y = rotate_points(Y,Ry); dY = Y(2,:)-Y(1,:);
      ay = ay + a;
    end
    else
      Rz = [cos(az) -sin(az) 0 ; sin(az) cos(az) 0 ; 0 0 1];
      X = rotate_points(X,Rz); dX = X(2,:)-X(1,:);
      Y = rotate_points(Y,Rz); dY = Y(2,:)-Y(1,:);
      Ry = [cos(ay) 0 sin(ay) ; 0 1 0 ; -sin(ay) 0 cos(ay)];
      X = rotate_points(X,Ry); dX = X(2,:)-X(1,:);
      Y = rotate_points(Y,Ry); dY = Y(2,:)-Y(1,:);
    end

    if(X(2,3)>X(1,3))
      a = pi;
      Ry = [cos(a) 0 sin(a) ; 0 1 0 ; -sin(a) 0 cos(a)];
      X = rotate_points(X,Ry); dX = X(2,:)-X(1,:);
      Y = rotate_points(Y,Ry); dY = Y(2,:)-Y(1,:);
    end

    % Remove HAs that are 'out of bounds':
    d1 = sqrt(sum((X(2,:)-Y(1,:)).^2)); %The base of one HA should not be too close to the head of another
    d2 = sqrt(sum((X(1,:)-Y(2,:)).^2));
    d3 = sqrt(sum((X(1,:)-Y(1,:)).^2)); %The heads of two HAs should not be too close together
    d4 = sqrt(sum((X(2,:)-Y(2,:)).^2)); %The base of two HAs should not be too close together
    trans_flag(n) = (Y(1,3)-X(1,3))>0;  %The head of the second HA must be above the head of the first for trans-crosslinking
    cis_flag(n) = (Y(1,3)-X(2,3))>0;    %The head of the second HA must be above the bottom of the first for cis-crosslinking
    if( (Y(1,3)<-100) + (d1<100) + (d2<100) + (d3<50) + (d4<50) );
      s(n) = nan;
      if(mod(n,10)==0)
        subplot(1,2,1);
        plot3(X(:,1),X(:,2),X(:,3),'-k','LineWidth',3); hold on;
        plot3(X(2,1),X(2,2),X(2,3),'sg','LineWidth',3);
        plot3(Y(:,1),Y(:,2),Y(:,3),'-','Color',[0.5 0.5 0.5],'LineWidth',3);
        plot3(Y(2,1),Y(2,2),Y(2,3),'s','Color',[0.7 0.7 0.7],'LineWidth',3);
      end
    else
      s(n) = sum(dX.*dY)./(sqrt(sum(dX.^2))*sqrt(sum(dY.^2)));
      if(mod(n,10)==0)
        subplot(1,2,1);
        plot3(X(:,1),X(:,2),X(:,3),'-k','LineWidth',3); hold on;
        plot3(X(2,1),X(2,2),X(2,3),'sg','LineWidth',3);
        plot3(Y(:,1),Y(:,2),Y(:,3),'-r','LineWidth',3);
        plot3(Y(2,1),Y(2,2),Y(2,3),'sg','LineWidth',3);
      end
    end
  end
  axis equal; title(pdb{k}); hold off;
  subplot(1,2,2)
  [temp,x] = hist(s,[-1:0.1:1]); 
  hist(s,[-1:0.1:1]); xlim([-1 1]); title(pdb{k});
  drawnow;

  i = find(s>0); cis(k) = sum(s(i).^p.*cis_flag(i))./N;        % Assign a weighted score 
  i = find(s<0); trans(k) = sum(s(i).^p.*trans_flag(i))./N;    % Assign a weighted score
end

savename = ['xlink_scores_psi',num2str(psi_range),'_phi',num2str(phi_range),'_theta',num2str(theta_range),'.mat'];
save(savename,'x','cis','trans','pdb');
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