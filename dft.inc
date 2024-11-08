      PARAMETER (maxel=5001,maxmon=1601)
      COMMON/VECT/dVex2A(0:maxel),dVex2B(0:maxel),
     *fdmon(0:maxel),VexA(0:maxel),
     *ebelam(0:maxel),VexB(0:maxel),diffVex(0:maxel),
     *P(0:maxel),bP(0:maxel),fdA(0:maxel),fdB(0:maxel),
     *PhiAA(0:maxel),convs(0:maxel),cdmon(0:maxel),eblam(0:maxel),
     *ehbclam(0:maxel),Phi(0:maxel)
      COMMON/VAR/dz,closew,pis,pif,pit,vk,rrT,vkrrT,hvk,scalem,scales,
     *dzpie,AA1,AA2,BB1,BB2,c1,c2,Y,emscale,eblemb,ehblcmb,eblsmb,bcdt,
     *rrjkdiff,threqz,rtwelve,pie,rthree,rdz,btrams,
     *sclosew,q1,q2,q3,p1,p2,p3,r2,r1,r0,s2,s1,s0,b2,b1,b0,r2sq,r1sq,
     *r0sq,Yfact,veq,rnmon,rrnmon,rrcmon,rq3,cdmbulk,cdsbulk,
     *cdmlbulk,cdslbulk,elblemb,elhblcmb,elblsmb,distp1,dmitt,eplag,
     *seplag,rrnarm,rnarm,tscalem,rp3,h,aveq,aYfact
      COMMON/HELTAL/istart,istp1,islut,ism,inw,nfack,imitt,nmon,
     *ist,ifin,istp1s,isluts,isms,inws,kst,kfin,ismms





