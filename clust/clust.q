
imed:{x?med x};
imax:{x?max x};
imin:{x?min x};




cure:{[newqu;com;numR;sam;pointDist] 
 
 qu:newqu[0];
 newq:newqu[1];
 
 
 old:qu j:j,qu[;`closestIdx]j@:imin qu[j:where qu`valid]`closestDist;
 qu[j;`valid]:0b;
 mean:avg pts:sam idx:raze old`idx;
 maxFromMean:idx imax sum each{x*x}mean-/:pts;
 
 rep:sam numR{[idx;pointDist;i]union[i]idx imax -0w^min pointDist[i;idx]}[idx;pointDist]/enlist maxFromMean;
 j0:first j;
 rep:(rep*1-com)+\:com*mean;
 
 
 k:(where qu`valid) except raze raze (old`furthestIdx);
 
 d:{min sum each x*x:raze x-/:\:y}
   [rep]each qu[k]`rep;
 
 new:`closestIdx`closestDist!raze (k;d)@\:imin d;
 
 qu[l;`closestIdx`closestDist]:(j0;d k?l:exec i from (qu k) where valid,closestDist>d);

 newq[j0;k]:d;
 
 newq[j1;j0]:d j2:k?j1:exec i from qu where valid,closestIdx in j; 
 n:where qu`valid;
 
 qu[j0]:update furthestIdx:enlist where (where qu`valid) in where (newq[j0])>med newq[j0], rep,idx,valid:1b from new;
 
 qu[j1;`closestIdx`closestDist]:{[newq;n;x] (newq[x]?min newq[x] n;min newq[x] n)}[newq;n]each j1;
 
 (qu;newq)}
 
cureClust2:{[sample;numRep;comp;numClust] 
 pointDist:{[fx;x;i]@[;i;:;0n]sum fx*fx-:x}[flip sample]'[sample;j:til count sample];
 queue:update idx:enlist each i,valid:1b from([]rep:enlist each sample);
 queue:queue,'flip `closestIdx`closestDist`furthestIdx!
    flip {[pointDist;x;y] (y[i];d i:imin d:pointDist[x]y;
     (enlist y[count sample]))}[pointDist]\:[j;j];
  newq: pointDist; `num xcols flip((`$/:string til count pointDist),`num)!(pointDist upsert til count pointDist);
  
  res:{[n;x]n<sum x[0]`valid}[numClust]cure[;comp;numRep;sample;pointDist]/(queue;pointDist);
   delete valid from update pts:sample idx from select from res[0] where valid}
