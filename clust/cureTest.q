imax:{x?max x};
imin:{x?min x};


plt:.p.import`matplotlib.pyplot;

cure:{[newqu;com;numR;sam;pointDist] 
 
 qu:newqu[0];
 newq:newqu[1];
 
 old:qu j:j,qu[;`closestIdx]j@:imin qu[j:where qu`valid]`closestDist;
 qu[j;`valid]:0b;
 
 mean:avg pts:sam idx:raze old`idx;
 maxFromMean:idx imax sum each{x*x}mean-/:pts;
 rep:sam numR{[idx;pointDist;i]union[i]idx imax -0w^min pointDist[i;idx]}[idx;pointDist]/enlist maxFromMean;
 rep:(rep*1-com)+\:com*mean;
 new:`closestIdx`closestDist!(k;d)@\:imin d:{min sum each x*x:raze x-/:\:y}[rep]each qu[k:where qu`valid]`rep;
 qu[l;`closestIdx`closestDist]:(j0:first j;d k?l:exec i from qu where valid,closestDist>d);

 qu[j0]:update rep,idx,valid:1b from new;
 
 newq[j0;k]:d;

 j1:exec i from qu where valid,closestIdx in j;

 newq[j1;j0]:d exec i from qu k where closestIdx in j;
 
 n:where qu`valid;


 qu[j1;`closestIdx`closestDist]:{[newq;n;x] (newq[x]?min newq[x] n;min newq[x] n)}[newq;n]each j1;

 (qu;newq)}
 
cureClust:{[sample;numRep;comp;numClust] 
  pointDist:{[fx;x;i]@[;i;:;0n]sum fx*fx-:x}[flip sample]'[sample;j:til count sample];
  queue:update idx:enlist each i,valid:1b from([]rep:enlist each sample);
  queue:queue,'flip`closestIdx`closestDist!
   flip{[pointDist;x;y]y[i],d i:imin d:pointDist[x]y}[pointDist]\:[j;j];

  res:{[n;x]n<sum x[0]`valid}[numClust]cure[;comp;numRep;sample;pointDist]/(queue;pointDist);
  delete valid from update pts:sample idx from select from res[0] where valid}

