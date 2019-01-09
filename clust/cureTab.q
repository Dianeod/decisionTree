\l p.q
\l cure.q

pyTimeAll:{[dataPt;dimCl;nCl;com;rep]
 pycure:.p.import[`pyclustering.cluster.cure];
 dim:dataPt,dimCl;
 d:(`int$(dataPt%nCl));
 clusters:nCl;
 init:d#dim#(*/[dim])?(`float$(first dim)*(last dim));
 merge:{[dim;init;x]flip @[flip init;x;+;rand((*/[dim]))*first dim]}[dim;init]each (clusters-1)#(til count flip init);

 clust:{[merge;x;y] x,merge[y]}[merge]/[init;til count merge];
 sample:clust;


 cureinst:pycure[`:cure][sample;nCl;rep;com;`ccore pykw 1b];
 a:.z.t;
 cureinst[`:process][];

 clusts:cureinst[`:get_clusters][]`;
 timePy:.z.t-a; 
 b:.z.t;
 cureQ:.cure.cure[rep;com;nCl;flip sample];
 timeQ:.z.t-b;
 compare:(asc asc each clusts)~asc asc each (cureQ);
 flip `dataPts`dim`Numclust`com`timePy`timeQ`rep`clust!(dataPt;dimCl;nCl;com;timePy;timeQ;rep;enlist compare)}


data:{[dataPt;dimCl;nCl]
 dim:dataPt,dimCl;
 d:(`int$(dataPt%nCl));
 clusters:nCl;
 init:d#dim#(*/[dim])?(`float$(first dim)*(last dim));
 merge:{[dim;init;x]flip @[flip init;x;+;rand((*/[dim]))*first dim]}[dim;init]each (clusters-1)?(til count flip init);

 clust:{[merge;x;y] x,merge[y]}[merge]/[init;til count merge];
 clust}


pyTime:{[dataPt;dimCl;nCl;com;rep;sample]
 pycure:.p.import[`pyclustering.cluster.cure];

 cureinst:pycure[`:cure][sample;nCl;rep;com;`ccore pykw 1b];
 a:.z.t;
 cureinst[`:process][];

 clusts:cureinst[`:get_clusters][]`;
 timePy:"i"$(.z.t-a) mod 1000;
 b:.z.t;
 cureQ:.cure.cure[rep;com;nCl;flip sample];
 timeQ:"i"$(.z.t-b) mod 1000;
 compare:(asc asc each clusts)~asc asc each cureQ;
 flip `dataPts`dim`Numclust`com`timePy`timeQ`rep`clust!(dataPt;dimCl;nCl;com;timePy;timeQ;rep;enlist compare)}

pyCl:{[nCl;com;rep;sample]
 pycure:.p.import[`pyclustering.cluster.cure];
 cureinst:pycure[`:cure][sample;nCl;rep;com;`ccore pykw 1b];
 cureinst[`:process][];
 clusts:cureinst[`:get_clusters][]`;
 clusts}
 
clustsQ:{[rep;com;nCl;sample] asc asc each .cure.cure[rep;com;nCl;flip sample]}

clustsPy:{[rep;com;nCl;sample] 
	pycure:.p.import[`pyclustering.cluster.cure];
	cureinst:pycure[`:cure][sample;nCl;rep;com;`ccore pykw 1b];
	cureinst[`:process][];
	clusts:cureinst[`:get_clusters][]`;
	asc asc each clusts}



