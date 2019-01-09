dsc:{[n;r;s](r+n?s-r)*/:(cos;sin)@\:4*n?acos 0};
data :flip(-5 1)+(1 1.5)*dsc[1000; 0;1.8],'dsc[2000;3.1;4.2],'dsc[4000;5.2;6.5];
data,:flip(4 -1)+(1 8)*dsc[4000; 0;1.];
data@:neg[4000&count data]?count data;

sample:data;


startPt:sample?Pt:(asc sample)(`int$(count sample)%2);


root:flip `idx`rep`left`right`dim`parent`rDim!
   (enlist startPt;enlist Pt;enlist 0;enlist 0;enlist 0;enlist startPt;enlist Pt[0]);


kdtree:{[root;sample;newSam]
 
 X:(count root)-1;
 
 $[0N!newSam[X;0]>Pt[0];(i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep from root where 
   right=1,parent=startPt);
   i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep from root where left=1,parent=startPt];
 
 $[(first i`idx)<>0N;
  (while[(first i`idx)<>0N,;a:i;i:$[first (newSam[X;first i`dim]>(raze i`rep)[first i`dim]);
    i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep from root where right=1,parent=first first a;
    i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep from root where left=1,parent=first first a]]);
    a:update idx:startPt,dim:0 from i]; 
  newa:update left:i`left,right:i`right  from a;
  root:root upsert flip update idx:(sample?newSam[X]),rep:enlist newSam[X],
   rDim:enlist newSam[X;((first a`dim)+1)mod 2],dim:enlist ((first a`dim)+1)mod 2,
   parent:enlist first a`idx from exec left,right from newa;
    show root;
   root}
