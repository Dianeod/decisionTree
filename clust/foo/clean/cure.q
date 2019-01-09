/ for tree nearest neighbors searches
\l kdtree.q

\d .cure
i.:(::)
clib:2:[`:./cure]
impf:{i[u]:clib(u:`$"_"sv string(x;y)),z}
impf'[`cure_cluster_dists;`F`E;3];
impf'[`cure_cluster_reps;`F`E;3];
impf'[`cure_nn;`F`E;4];

/ distances between x - cluster1, y - list of other clusters, z - list of representative points, x and y are indexes into z
ccdists:{$["f"=.Q.t type z 0;i.cure_cluster_dists_F;"e"=t;i.cure_cluster_dists_E;'`type][x;y;z]}
/ shrink points towards their mean
shrink:{y+x*y-\:avg y}
/ x - representative points shrunk by y - compression towards mean from z - a list of points
crep:{t$shrink[y]$["f"=t;i.cure_cluster_reps_F;"e"=t;i.cure_cluster_reps_E;'`type][z;x;(t:.Q.t type z 0)$avg z]}
/ from pi - a list of indices into reps - a list of representative points, find the closest (other) cluster and distance to it
ccnn:{[pi;tree;clusters;reps]@[;0;clusters]{mi:u?mn:min u:x[;1];(x[;0]mi;mn)}cnn[;tree;clusters;reps]each pi}
/ closest (other) cluster and distance to it from x - a point index (into w - a list of reps) given a tree y and mapping from 
/ w - a list 
cnn:{[x;y;z;w]$["f"=.Q.t type w 0;i.cure_nn_F;"e"=t;i.cure_nn_E;'`type][x;y;z;w]}
imin:{x?min x}


cure:{[r;c;n;d]
 reps:flip d;
 tree:.kdtree.create[r;d];
 r2l:((pc:count d 0)#0N){[tree;x;y]@[x;tree[3]y;:;y]}[tree]/where tree 2; / rep to leaf node it's in
 r2c:til pc; / rep to cluster it's in
 gone:pc#0b; / cluster is no longer valid
 c2r:enlist each r2c; / cluster to its representatives
 c2p:enlist each r2c; / cluster to its points
 ndists:flip cnn[;tree;r2c;reps]each r2c; / nearest neighbors and distances to each for each initial representative
 i:0;N:pc-n;   / loop counter, and number of iterations required
 ST1:.z.P;
 while[i<N;
  mci:u,ndists[0;u:imin ndists 1];                               / clusters to merge
  ori:raze c2r mci;                                              / old representatives
  orl:r2l ori:raze c2r mci;                                      / old representatives and the leaf nodes they belong to
  npts:flip d[;npi:raze c2p mci];                                / points in new cluster
  nreps:crep[r;c;npts];                                          / representatives of new cluster
  reps[nri:count[nreps]#ori]:nreps;                              / overwrite any old representatives with the new ones
  r2l[nri]:nrl:.kdtree.searchfrom[tree;;0]each nreps;            / leaf nodes for new representatives
  r2c[nri]:r2c ori 0;                                            / new cluster is first of old clusters
  c2p[mci]:(npi;0#0);c2r[mci]:(nri;0#0);                         / update cluster -> points and representatives
  gone[mci 1]:1b;                                                / mark second of merged clusters as removed
  / update the tree with the new representatives, probably a faster way
  tree:.[tree;(3;distinct orl);{y except x}ori];                 / delete old representatives
  tree:tree{.[x;(3;y 0);{y,x}y 1]}/flip(nrl;nri);                / add new representatives 
  / now update cluster distances and closest cluster for everything
  cnc:ccnn[nri;tree;r2c;reps];                                   / nearest cluster to the newly merged cluster and the distance to it
  u:where not[gone]and(ndists[1]>cnc 1)and not ndists[0]in mci;  / clusters which could be nearer to the new cluster than their current closest
  cdists:@[(pc#0w);u;:;ccdists[nri;;reps]c2r u];                 / calculate cluster distance for each cluster which could be closer
  tab:([]cdist:cdists;ndist:ndists 1;nc:ndists 0);        
  nearer:exec i from tab where ndist>cdist;
  invalid:exec i from tab where nc in mci,not i in nearer,not gone i;
  ndists[0 1;nearer]:(count[nearer]#mci 0;cdists nearer);
  ndists[0 1;invalid]:$[count invalid;flip ccnn[;tree;r2c;reps]each c2r invalid;(0#0;0#0f)];  
  ndists[;mci 0]:cnc;
  ndists[;mci 1]:(0N;0w);
  i+:1;
 ];
 :c2p where not gone;
 }

cure1:{[r;c;n;d]
 reps:flip d;
 tree:.kdtree.create[r;d];
 }
 
\
