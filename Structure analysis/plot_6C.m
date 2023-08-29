function plot_6C;

load('xlink_scores_psi60_phi30_theta30_manuscript.mat')

% Filter to remove PDB IDs of duplicate mAbs:
load duplicate_ids.mat
q = 1;
for k = 1:length(pdb)
  current = pdb{k};
  s = sum(dup == current,2);
  if (sum(s==4)==0);
    PDB{q}=pdb{k};
    CIS(q) = cis(k);
    TRANS(q) = trans(k);
    q = q + 1;
  end
end
pdb = PDB; cis = CIS; trans = TRANS;

cis(cis<1e-6) = 1e-6;
trans(trans<1e-6) = 1e-6;
c0 = 0.8;

loglog(cis,trans,'s',MarkerFaceColor=[c0 c0 c0], MarkerEdgeColor='none',MarkerSize=8);
hold on; axis square

test = '4fqr';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor='r',MarkerSize=8); hold on;
    disp(test)
  end
end

test = '4gms';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0 0.5 0.5],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '4fqi';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0 1 0.5],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '3ztj';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[.5 0 .5],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '6oc3';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[1 0.5 0],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '6hjp';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0 0.2 .6],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '3sdy';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor='m',MarkerSize=8); hold on;
    disp(test)
  end
end

test = '5k9o';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0.3 0.3 0.3],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '5w6g';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0.3 0.3 0.3],MarkerSize=8); hold on;
    disp(test)
  end
end

test = '4o58';
for k = 1:length(cis)
  current = pdb{k};
  if(test==current);
    loglog(cis(k),trans(k),'sk',MarkerFaceColor=[0.3 0.3 0.3],MarkerSize=8); hold on;
    disp(test)
  end
end

xlim([1e-4 1])
ylim([0.6e-6 1])