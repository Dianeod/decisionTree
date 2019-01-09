np:.p.import`numpy
shape:{-1_count each first scan x};
readsam:.p.import[`pyclustering.utils]`:read_sample;
SIMPLE_SAMPLES:.p.import[`pyclustering.samples.definitions]`:SIMPLE_SAMPLES;
FCPS_SAMPLES:.p.import[`pyclustering.samples.definitions]`:FCPS_SAMPLES;
imax:{x?max x};
imin:{x?min x};
plt3D:.p.import[`mpl_toolkits.mplot3d]`:Axes3D;

plt:.p.import[`matplotlib]`:pyplot;

cure:{[qu;com;numR;sam;pointDist]
  old:qu j:j,qu[;`closestIdx]j@:imin qu[j:where qu`valid]`closestDist;
  qu[j;`valid]:0b;
  mean:avg pts:sam idx:raze old`idx;
  maxFromMean:idx imax sum each{x*x}mean-/:pts;
  rep:sam numR{[idx;pointDist;i]union[i]idx imax -0w^min pointDist[i;idx]}[idx;pointDist]/enlist maxFromMean;
  rep:(rep*1-com)+\:com*mean;
  new:`closestIdx`closestDist!(k;d)@\:imin d:{min sum each x*x:raze x-/:\:y}[rep]each qu[k:where qu`valid]`rep;
  qu[l;`closestIdx`closestDist]:(j0:first j;d k?l:exec i from qu where valid,closestDist>d);
  qu[j0]:update rep,idx,valid:1b from new;
  j:exec i from qu where valid,closestIdx in j;
  qu[j;`closestIdx`closestDist]:flip(kk@'l;d@'l:imin each d:{min sum each x*x:raze x-/:\:y}/:'[qu[j]`rep]qu[`rep]kk:(k,j0)except/:j);
  qu};

cureClust:{[sample;numRep;comp;numClust]
  pointDist:{[fx;x;i]@[;i;:;0n]sum fx*fx-:x}[flip sample]'[sample;j:til count sample];
  queue:update idx:enlist each i,valid:1b from([]rep:enlist each sample);
  queue:update clustDist:pointDist from queue,'flip`closestIdx`closestDist!
   flip{[pointDist;x;y]y[i],d i:imin d:pointDist[x]y}[pointDist]\:[j;j];

  res:{[n;x]n<sum x`valid}[numClust]cure[;comp;numRep;sample;pointDist]/queue;
  delete valid from update pts:sample idx from select from res where valid}


