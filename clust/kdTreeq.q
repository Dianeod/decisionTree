\l p.q
np:.p.import`numpy
readsam:.p.import[`pyclustering.utils]`:read_sample;
SIMPLE_SAMPLES:.p.import[`pyclustering.samples.definitions]`:SIMPLE_SAMPLES;
FCPS_SAMPLES:.p.import[`pyclustering.samples.definitions]`:FCPS_SAMPLES;
imax:{x?max x};
imin:{x?min x};
plt3D:.p.import[`mpl_toolkits.mplot3d]`:Axes3D;

plt:.p.import[`matplotlib]`:pyplot;



/insert cluster into tree

insertKd:{[kd;sample2;L;clust] 
 
 
 dir:{(first x[0]`idx)>=0}{[sample2;kd;x] a:x[0];
  $[first (sample2[first x[0]`dim]>(first x[0]`rDim));
  i:update left:enlist 0,right:enlist 1 from exec idx,dim,rDim from kd where right=1,parent=first a`idx;
  i:update left:enlist 1,right:enlist 0 from exec idx,dim,rDim from kd where left=1,parent=first a`idx];
  (i;a)}[sample2;kd]/(i;i:exec from kd where parent=0,left=0,right=0); /insert cluster into the tree by looking at splitting dimension of each node &going left or right

 
 i:dir[0]; /indicates whether to the left or right 
 a:dir[1]; /parent node
 dim:((first a`dim)+1)mod 2; /get its new splitting dimension

 
 root:kd upsert flip update idx:(max kd`idx)+1,initi:L,clust:clust,rep:enlist sample2,
    rDim:enlist sample2[dim],dim:enlist dim,
    parent:enlist first a`idx from exec left,right from i; /update the info of new node into tree
 root}


/Calculating distances between clusters using kd tree. 


distC:{[kd;pt] 
   
 distCalc:{[kd;query;bestD]
   
 X:bestD[3]; / nodes that were already searched, not to be searched again
    
 cl:bestD[5]; / what cluster the node belongs to, as not to search any points in the same cluster
 n:bestD[1]; /nodes to search
 a:n where {[cl;kd;x](first exec clust from kd where idx=x)<>cl}[cl;kd]each n; /nodes to search that arent in the same cluster
  
 newD:imins,a i?imins:min i:{[kd;query;x] 
    sum m*m:(raze exec rep from kd where idx=x)-query}[kd;query] each a; /get minimum dist of all searched nodes

  $[(newD[0]<bestD[0])&(count a)<>0;(bestD[0]:newD[0];bestD[2]:newD[1])
     ;(bestD:bestD)]; /if new dist is less than current best dist, then that becomes new best dist
  
  axisD:(raze {[kd;bestD;query;x] $[(m*m:
    (first exec rDim from kd where idx=x)-
    query(first exec dim from kd where idx=x))<=bestD[0];
    (exec idx from kd where parent=x),exec parent from kd where idx=x;
    $[(query(first exec dim from kd where idx=x))<
    first exec rDim from kd where idx=x;
    (exec idx from kd where parent=x,left=1),exec parent from kd where idx=x;
    (exec idx from kd where parent=x,right=1),exec parent from kd where idx=x]]
    }[kd;bestD;query]each n)except bestD[3]:X,n; /get dists between node and search pts based on splitting dimension.
    /if =< than best Dist, than search the children of that node &parents.
    /Go up the tree and to the left/right based on whether that pt is <or> search pt
  
 (bestD[0];distinct axisD;bestD[2];bestD[3];bestD[4];cl)};
 dist:{[distCalc;kd;X] ({(count x[1])<>0}distCalc[kd;
    raze exec rep from kd where idx=X]/(0W;
    (raze (exec idx from kd where parent=X),exec parent from kd where idx=X)except X;X;X;X
    ;first exec clust from kd where idx=X))
    }[distCalc;kd;pt];
    
    kdC:update closDist:dist[0],closIdx:dist[2] from kd where idx=pt; /update new closDist and idx in tree
    kdC}



/ delete from kd tree

deleteN:{[kd;X] 
 n:first exec idx from kd where initi=X;
 delN:{
  kd:x[0]; /kdtree
  X:x[1]; /point to be deleted
 
  delNode:select from kd where idx=X; /details of deleted pt
  axis:delNode`dim; /splitting dim

  mindim:$[(count exec idx from kd where right=1,parent=X)=0;
    raze {[kd;x](count exec idx from kd where
       parent=first x)<>0}[kd]{[kd;x]raze exec idx from kd where parent=last x}[kd]\
    first exec idx from kd where parent=X,left=1; / if has no right child then left child replaces 
    raze {[kd;x] (count exec idx from kd where
       parent=first x)<>0}[kd]{[kd;x]raze exec idx from kd where parent=last x}[kd]\
    first exec idx from kd where parent=X,right=1]; / get all the right children if there

  newP:mindim repPt?min repPt:(raze{[kd;x] exec rep from kd where idx=x}[kd]each mindim)axis; /the min value of rep pts from mindim based on splitting dimension
  newNode:select from kd where idx=newP; /get info from kdtree of newP
  tree:update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,
    clust:newNode`clust,clustIdx:newNode`clustIdx,rDim:(first newNode`rep)axis,
    closIdx:newNode`closIdx from kd where idx=X; /newNode replaces the deleted pt info in the tree, but dim and idx stay the same as deleted pt
  

  tree:update closIdx:X from tree where closIdx=first newNode`idx; /any pt that had the replacing pt as its closest idx has to update closIdx to its new value in the tree
  

    (tree;newP)};

 delCl:{(count select from (first x) where parent=last x)<>0}delN/(kd;n); /repeat this until reach a node with no children
 
 
 delete from first delCl where idx=last delCl} /delete the node with no children



/Create tree, search init nearest neighbours
createTree:{[sample]

 root:flip `idx`initi`rep`left`right`dim`parent`rDim`clust`clustIdx!
   (0;0;enlist sample[0];enlist 0;enlist 0;enlist 0;enlist 0;enlist sample[0;0];0;enlist (0 0)); /insert first cluster

 kds:insertKd/[root;(1_sample);1_til count sample;1_til count sample]; /insert the rest of the clusters

 kds:update clustIdx:enlist each til count sample from kds; /insert the cluster indices

 kds:distC/[kds;kds`initi]; /get closest cluster to each cluster
 kds
 }


clust:{[sample;numR;com;kd]
 j:first select from kd where closDist=first min closDist; /get min dist from tree
 
 j2:select clust from kd where idx=first j`closIdx; /get cluster of the closestPt
 
 old:select from kd where clust in ((j2`clust),j`clust); /get info of two clostest pts

 j0:(exec initi from kd where closIdx in old`idx); /get the initial indexes of those pts
 
 mean:avg pts:sample idxs:(distinct raze old`clustIdx);
 
 deleteClust:deleteN/[kd;idxs]; /delete merged clusters from tree
 
 maxFromMean:idxs imax sum each{x*x}mean-/:pts; /get pt with max dist from mean of cluster
 
 rep:distinct sample sami:numR{[sample;idxs;x] x,maxIdx imax {[sample;maxMean;x] min{[sample;x;maxMean]sum x*x:sample[maxMean]-sample[x]
     }[sample;x] each maxMean}[sample;x]each maxIdx:idxs except x}[sample;idxs]/maxFromMean;
  /get rep pts based on the most spread out pts in cluster

 

 rep:(rep*1-com)+\:com*mean; /move rep pts based on compressed
 
 newCl:{[rep;x]min (sum each x*x:rep-\:x)}[rep]each deleteClust`rep; /dist of new cluster to all other clusters

 j3:(deleteClust`idx)[n:where newCl<deleteClust`closDist]; /any clusters where a<closDist

 closCl: (deleteClust`idx)where newCl=closD:min newCl; /what cluster is closest to the new clusters

 insertClust:insertKd/[deleteClust;rep;sami;first idxs]; /insert new rep pts into kdtree

 insertClust2:{[n;j3;newCl;insertClust;x] update closDist:newCl(n x),closIdx:last (insertClust`idx) from 
    insertClust where idx=j3(x)}[n;j3;newCl]/[insertClust;til count j3]; /update the tree with j3

 clustDist:update closIdx:first closCl,closDist:closD from insertClust2 where clust=first idxs; /update closDist of new clusters

 j5:exec idx from clustDist where clust=first idxs; /indexes of new clusters

  j4:(exec idx from clustDist where initi in j0) except j5; /any clusters that had the merged clusters as closDist

  recalc:distC/[clustDist;j4]; /recalc those distances
  insertIdx:enlist idxs; /initial indexes of merged cluster
 
  kd:{[insertIdx;kd;x] update clustIdx:(insertIdx) from kd where initi=x}[insertIdx]/[recalc;sami]; /update into kdtree
 
  kd}



cure:{[sample;numR;com;numCl] 
 cureTab:{[numCl;kd] (count distinct kd`clust)>numCl}[numCl]clust[sample;numR;com]/createTree sample;
 distinct cureTab`clustIdx}


