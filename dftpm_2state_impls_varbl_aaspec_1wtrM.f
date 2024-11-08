      program platem
      implicit double precision (a-h,o-z)
      include 'dft.inc'
      dimension c(0:maxel,maxmon),chA(0:maxel),chB(0:maxel),
     *ttP(0:maxel),VexS(0:maxel),VAlrc(0:maxel),vBlrc(0:maxel),
     *VSlrc(0:maxel)
      write(*,*) 'SIMPLE L-J 1-component 2-state system!'    
      write(*,*) 'implicit solvent model'       
      write(*,*) 'eAA specifically included!'       
      ifc = 38
      ins = 49
      iPd = 50
      iamu = 43
      ibdm = 45
      pi = acos(-1.d0)
      bk = 1.38066D-23
      avno = 6.02214D23
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0
      es22 = -32.d0*pi/9.d0 
      ddtol = 0.000002d0
      open (ifc,file='fcdfil',form='formatted')    
      open (iPd,file='Pdistr',form='formatted')    
      open (ins,file='input.dftp',form='formatted')
      open (iamu,file='amufil',form='formatted')
      open (ibdm,file='bdmfil',form='formatted')      
      rewind ifc
      rewind ins
      rewind iPd
      rewind iamu
      rewind ibdm
      read(ins,*) bdm
c      read(ibdm,*) bdm
      
      read(ins,*) gA,gB
      read(ins,*) T
      read(ins,*) eBB,eAB
      read(ins,*) eAS,eBS
      read(ins,*) h
      bdist = h

      read(ins,*) refdz
      read(iamu,*) amu
      read(iamu,*) dz

      read(iamu,*) bl            
      
      es22 = es22+amu

      read(ins,*) dmm
      read(ins,*) dms
      read(ins,*) kread
      read(ins,*) hwdisttrams
      read(ins,*) Af,Bf
      read(ins,*) zwc,zpmin
      read(ins,*) maxniter
      read(ins,*) nmon
      read(ins,*) Rctrams
      read(ins,*) eAA,Sf
      bdtot = 1.d0
      bds = 1.d0-bdm
      write(*,*) 'bdm = ',bdm
      write(*,*) 'gA,gB = ',gA,gB
      write(*,*) 'T = ',T
      write(*,*) 'eBB,eAB = ',eBB,eAB
      write(*,*) 'eBS,eAS = ',eBS,eAS
      write(*,*) 'h = ',h
      write(*,*) 'dz = ',dz
      write(*,*) 'dmm,dms = ',dmm,dms
      write(*,*) 'kread = ',kread
      write(*,*) 'Af,Bf = ',Af,Bf
      write(*,*) 'zwc,zpmin = ',zwc,zpmin
      write(*,*) 'maxniter = ',maxniter
      write(*,*) 'nmon = ',nmon
      write(*,*) 'Rctrams = ',Rctrams
      write(*,*) 'eAA = ',eAA
      write(*,*) 'Sf (FROM INPUT!) = ',Sf
      write(*,*) 'bl = ',bl      
      rrT = 1.d0/T 
      closew = 0.d0
      ufact = 1.d0
      q1 = 1.d0
      q2 = 1.d0
      q3 = 1.d0
      rq3 = 1.d0
      sclosew = closew
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
      rdz = 1.d0/dz

      rbl = 1.d0/bl
c      twopidz = twopi*dz
      twopidz = 0.5d0*dz*rbl
      
      irdz = int(rdz+0.001d0)
      nfack = int(h/dz+0.01d0)
      istart = int(closew/dz+0.01d0)
      istp1 = istart+1 
      islut = int((h-closew)/dz+0.01d0)  
      ists = int(sclosew/dz+0.01d0)
      istp1s = ists+1 
      isluts = int((h-sclosew)/dz+0.01d0)
      
c      ism = int(1.d0/dz+0.01d0)
c     isms = int(q1/dz+0.01d0)
      ism = int(bl/dz+0.01d0)
      isms = int(bl/dz+0.01d0)
      
      inw = ism+int(closew/dz+0.01d0)
      inws = isms+int(sclosew/dz+0.01d0)      
      pie = pi/8.d0
      dzpie = pie*dz

      
c     imitt = nfack/2
      imitt = nfack

      
c      es22AA = es22*rrT
      es22 = es22*rrT
      aphi = aphi*rrT
      rphi = rphi*rrT
      es22AA = es22*eAA
      es22BB = es22*eBB
      es22AB = es22*eAB
      es22BS = es22*eBS
      es22AS = es22*eAS
c      phzAA = -(6.d0*pi/5.d0)*rrT
      phz = -(6.d0*pi/5.d0)*rrT
      phzAA = phz*eAA
      phzBB = phz*eBB
      phzAB = phz*eAB
      phzBS = phz*eBS
      phzAS = phz*eAS
      write(*,*) 'es22AA = ',es22AA 
      write(*,*) 'es22BB = ',es22BB 
      write(*,*) 'es22AB = ',es22AB 
      write(*,*) 'es22AS = ',es22AS
      write(*,*) 'es22BS = ',es22BS 
      rkaj = (es22AB-0.5d0*es22AA-0.5d0*es22BB)
      write(*,*) 'rkaj = ',rkaj
c      twokaj = 2.d0*eAB-1.d0-eBB
      twokaj = 2.d0*eAB-eAA-eBB
      write(*,*) 'twokaj = ',twokaj
      Pb = 0.5d0
      nPb = 100000
      do i = 1,100000
      tPb = gA*(1.d0-Pb)*
     *dexp(bdm*(Pb*2.d0*rkaj+es22BB-es22AB)+(1.d0-bdm)*(es22BS-es22AS))/
     *gB
      Pb = 0.999d0*Pb+0.001d0*tPb
      if (mod(i,(nPb/10)).eq.0) write(*,*) i,Pb,tPb
      if (i.gt.(nPb-2)) write(*,*) i,Pb,tPb
      enddo
      write(*,*) 'dabs(Pb-tPb) = ',dabs(Pb-tPb)
      write(*,*) 'Pb = P(state A, in the bulk solution) = ',Pb
      bPb = 1.d0-Pb
      write(*,*) 'bPb = P(state B, in the bulk solution) = ',bPb

      rnmon = dfloat(nmon)
      rrnmon = 1.d0/rnmon
      bdpol = bdm/rnmon

      chempp=rnmon*bdm*(2.d0*Pb*bPb*es22AB+Pb*Pb*es22AA+bPb*bPb*es22BB)+
     *rnmon*(Pb*dlog(Pb/gA)+bPb*dlog(bPb/gB))+dlog(bdpol)+
     *rnmon*(Pb*es22AS-2.d0*bdm*Pb*es22AS+
     *bPb*es22BS-2.d0*bdm*bPb*es22BS)-rnmon*dlog(1.d0-bdm)
      emscale = chempp*rrnmon
      chemps = -rrnmon*chempp
c      pressure = 0.5d0*bdm*bdm*(2.d0*Pb*bPb*es22AB+Pb*Pb*es22AA+
c     *bPb*bPb*es22BB)+bdm*(1.d0-bdm)*(Pb*es22AS+tPb*es22BS)+bdpol+bds
      write(*,*) 'chempp  = ',chempp
c      write(*,*) 'pressure = ',pressure
      write(*,*) 'h,dz = ',h,dz
      write(*,*) 'chemps (=-chempp/rnmon)  = ',chemps
      write(*,*) 'closew,sclosew = ',closew,sclosew
      write(*,*) 'nfack = ',nfack
      write(*,*) 'istp1,islut,inw = ',istp1,islut,inw
      write(*,*) 'istp1s,isluts,inws = ',istp1s,isluts,inws
      write(*,*) 'ism,isms = ',ism,isms
      write(*,*) 'dmm,dms (density mixing param. mon.,solv.) = ',dmm,dms

      write(*,*) 'TRUNC. + SH. STANDARD (NOT STEEP) MORSE WALL POT.!'
      write(*,*)'w(z)=Af*((1.d0-exp(-(z-zpmin)))^2-1.)-w(zwc),(z<zwc)'
      write(*,*) 'Af,Bf = ',Af,Bf
      write(*,*) 'zpmin = ',zpmin
      write(*,*) 'zwc = ',zwc
      Af = Af*rrT
      Bf = Bf*rrT
      Sf = Sf*rrT
      write(*,*) 'Af,Bf (in kT units)= ',Af,Bf
      write(*,*) 'Sf (in kT units)= ',Sf
      refVex = Af*((1.d0-dexp(-(zwc-zpmin)))**2-1.d0)
      BrefVex = Bf*((1.d0-dexp(-(zwc-zpmin)))**2-1.d0)
      SrefVex = Sf*((1.d0-dexp(-(zwc-zpmin)))**2-1.d0)
      write(*,*) 'refVex = ',refVex
      write(*,*) 'BrefVex = ',BrefVex
      write(*,*) 'SrefVex = ',SrefVex
      hwdist = zpmin
      write(*,*) 'hwdist (=zpmin) = '
      write(*,'(1e25.14)') hwdist
      z = closew-0.5d0*dz
      rewind 78
      bdA = bdm*Pb
      bdB = bdm*bPb
      do iz = istp1s,isluts
      z = z+dz
      zw = z+hwdist
      zh = h-z
      zhw = h+hwdist-z
      VexA(iz) = 0.d0
      VexB(iz) = 0.d0
      VexS(iz) = 0.d0
c      dVex2A(iz) = 0.d0 
      if (zw.lt.zwc) then
      VexA(iz) = Af*((1.d0-dexp(-(zw-zpmin)))**2-1.d0)-refVex
      VexB(iz) = Bf*((1.d0-dexp(-(zw-zpmin)))**2-1.d0)-BrefVex
      VexS(iz) = Sf*((1.d0-dexp(-(zw-zpmin)))**2-1.d0)-SrefVex
      endif
      write(78,'(4e18.7)') z,VexA(iz),VexB(iz),VexS(iz)      
      enddo


      write(*,*) 'VexA(istp1) = ',VexA(istp1) 
      write(*,*) 'VexB(istp1) = ',VexB(istp1) 
      write(*,*) 'VexS(istp1) = ',VexS(istp1) 
c      write(*,*) 'diffVex(istp1) = ',diffVex(istp1) 
c      write(*,*) 'diffVex(imitt) = ',diffVex(imitt) 
      write(*,*) 'inw = ',inw,ism,imitt
      if (kread.eq.0) then
      do 12 iz = istp1,imitt
      P(iz) = Pb
      bP(iz) = 1.d0-P(iz)
      fdA(iz) = P(iz)*bdm
      fdB(iz) = bP(iz)*bdm
 12   fdmon(iz) = bdm
      else
      do iz = istp1,imitt
      read(ipD,*) trams,P(iz)
      read(ifc,*) trams,fdmon(iz),tramsA,tramsB
      bP(iz) = 1.d0-P(iz)
      fdA(iz) = P(iz)*fdmon(iz)
      fdB(iz) = bP(iz)*fdmon(iz)      
      enddo
      endif

      ifin = imitt+ism+1
      ist = imitt-ism
      rrjkdiff = 1.d0/dfloat(2*ism+1)      
      
      webdA = fdA(imitt-ism)
      webdB = fdB(imitt-ism)
      webdm = fdmon(imitt-ism)
      weP = P(imitt-ism)
      webP = bP(imitt-ism)
      do iz = imitt-ism+1,imitt+ism
      rkk = dfloat(iz-ist)
      rjj = dfloat(ifin-iz)
      fdA(iz) = (rjj*webdA+rkk*bdA)*rrjkdiff
      fdB(iz) = (rjj*webdB+rkk*bdB)*rrjkdiff
      fdmon(iz) = (rjj*webdm+rkk*bdm)*rrjkdiff
      P(iz) = (rjj*weP+rkk*Pb)*rrjkdiff
      bP(iz) = 1.d0-P(iz)
      enddo

      qqq = bdm*(2.d0*Pb*bPb*es22AB+Pb*Pb*es22AA+bPb*bPb*es22BB)+   
     *(Pb*es22AS-2.d0*bdm*Pb*es22AS+bPb*es22BS-2.d0*bdm*bPb*es22BS)
      yyy = Pb*dlog(Pb/gA)+bPb*dlog(bPb/gB)-dlog(1.d0-bdm)
      eblemb = dexp(emscale-qqq-yyy)
      ehblemb = dsqrt(eblemb)

      do iz = imitt+ism+1,maxel
      fdA(iz) = bdA
      fdB(iz) = bdB
      fdmon(iz) = bdm
      P(iz) = Pb
      bP(iz) = 1.d0-Pb
      ebelam(iz) = eblemb
      ehbclam(iz) = ehblemb
      enddo
      
      z = 0.d0
      iz = 0
      zp = z-dz
      do 212 jz = 0,maxel
      zp = zp+dz
      iii = iabs(jz-iz)
      diffz = abs(z-zp)
      if (diffz.gt.1.d0) then
      diffzsq = diffz**2
      rrdiffz2 = 1.d0/diffzsq
      rrdiffz4 = rrdiffz2**2
      rrdiffz10 = rrdiffz2*rrdiffz4**2
      phizmm = rphi*rrdiffz10-aphi*rrdiffz4
      else
c      phizmm = phzAA
      phizmm = phz
      endif
c 212  PhiAA(iii) = phizmm*dz
 212  Phi(iii) = phizmm*dz
c      twopidz = 0.5d0*dz
      rbdm = 1.d0/bdm

      rewind 76
      rewind 74      
      z = closew-0.5d0*dz      
      do i = istp1s,imitt
      z = z+dz

      uA = 0.d0
      um = 0.d0
      uB = 0.d0
      uu = 0.d0
      uP = 0.d0
      do j = imitt+10*ism+1,maxel
      iii = iabs(j-i)
      ttt = Phi(iii)
      uA = bdA*ttt+uA
      um = bdm*ttt+um
      uB = bdB*ttt+uB
      uu = ttt+uu
      uP = Pb*ttt+uP
      enddo
         
      zhw = dfloat(maxel)*dz-z      
      szhw3 = (1.d0/zhw)**3
      szhw9 = szhw3*szhw3*szhw3
      refV = (2.d0*szhw9/45.d0-szhw3/3.d0)
c      VAlrc(iz) = (2.d0*szhw9/45.d0-szhw3/3.d0)*((bdA*eAA+bdB*eAB)+
c     *eAS*(1.d0-bdm))
c      VBlrc(iz) = (2.d0*szhw9/45.d0-szhw3/3.d0)*((bdB*eBB+bdA*eAB)+
c     *eBS*(1.d0-bdm))
cc     alpha = -emscale+
cc     *VexA(i)*P(i)+VexB(i)*bP(i)-VexS(i)+
cc     *P(i)*(uA*eAA+uB*eAB)+bP(i)*(uB*eBB+uA*eAB)+
cc     *eAS*(P(i)*uS-uA)+eBS*(bP(i)*uS-uB)+
cc     *P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB)
c      write(76,'(4f21.12)')z,VAlrc(i),VBlrc(i),VSlrc(i)
      write(74,'(4e21.12)')z,refV*bdA,refV*bdB,refV*(1.d0-bdm)
      VAlrc(i) = uA+refV*bdA
      VBlrc(i) = uB+refV*bdB
      VSlrc(i) = uu-um+refV*(1.d0-bdm)
      write(76,'(4f21.12)')z,VAlrc(i),VBlrc(i),VSlrc(i)      
      enddo
      write(*,*) zhw,szhw3,refV
      write(*,*) refV*bdA,refV*bdB,refV*(1.d0-bdm)     
c      stop

      
      ddmax = 10000.
      dpmax = 10000.
      niter = 0
 100  continue
      niter = niter+1
      if (mod(niter,100).eq.0) write(*,*) niter,ddmax,dpmax
      if (niter.gt.maxniter) goto 200
      if (ddmax.lt.ddtol.and.dpmax.lt.ddtol) goto 200
      ddmax = 0.d0
      dpmax = 0.d0
      do i = istp1s,imitt-ism
      fdm = fdmon(i)
      uA = 0.d0
      um = 0.d0
      uB = 0.d0
      uu = 0.d0
      uP = 0.d0
      do j = istp1,islut+ism
      iii = iabs(j-i)
c      ttt = PhiAA(iii)
      ttt = Phi(iii)
      uA = fdA(j)*ttt+uA
      um = fdmon(j)*ttt+um
      uB = fdB(j)*ttt+uB
      uu = ttt+uu
c      uP = P(j)*ttt+uP
      enddo

      do j = islut+ism+1,islut+10*ism
      iii = iabs(j-i)
c      ttt = PhiAA(iii)
      ttt = Phi(iii)
      uA = bdA*ttt+uA
      um = bdm*ttt+um
      uB = bdB*ttt+uB
      uu = ttt+uu
c      uP = Pb*ttt+uP
      enddo
      
c      ubP = uu-uP
      uS = uu-um
      uA = uA+VAlrc(i)
      uB = uB+VBlrc(i)
      uS = uS+VSlrc(i)
      alpha = -emscale+
     *VexA(i)*P(i)+VexB(i)*bP(i)-VexS(i)+
     *P(i)*(uA*eAA+uB*eAB)+bP(i)*(uB*eBB+uA*eAB)+
     *eAS*(P(i)*uS-uA)+eBS*(bP(i)*uS-uB)+
     *P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB)
c     *P(i)*VAlrc(i)+bP(i)*VBlrc(i)  
      ebelam(i) = (1.d0-fdm)*dexp(-alpha)
      ehbclam(i) = dsqrt(ebelam(i))
      salpha = alpha

      alpha = uA*eAA+uB*eAB-uA*eAB-uB*eBB+
     *eAS*uS-eBS*uS+
     *VexA(i)-VexB(i)
c     *P(i)*VAlrc(i)+bP(i)*VBlrc(i) 
      tP = (1.d0-P(i))*gA*dexp(-alpha)/gB
      if (tP.gt.1.d0) tP = 1.d0-ddtol
      if (tP.lt.0.d0) tP = ddtol
      ddiff = abs(tP-P(i))/tP
      if (ddiff.gt.dpmax) dpmax = ddiff
      P(i) = P(i)*dms+tP*(1.d0-dms)
      bP(i) = 1.d0-P(i)

      enddo

c      write(*,*) i,imitt-ism
c      write(*,*) salpha,dexp(salpha),fdm
c      write(*,*) ebelam(imitt-ism),eblemb
c      stop

      webe = ebelam(imitt-ism)
c      write(*,*) 'webe,eblemb = ',webe,eblemb      
      webeh = ehbclam(imitt-ism)
      weP = P(imitt-ism)
      webP = bP(imitt-ism)
      do iz = imitt-ism+1,imitt+ism
      rkk = dfloat(iz-ist)
      rjj = dfloat(ifin-iz)
      ebelam(iz) = (rjj*webe+rkk*eblemb)*rrjkdiff
      ehbclam(iz) = (rjj*webeh+rkk*ehblemb)*rrjkdiff
      P(iz) = (rjj*weP+rkk*Pb)*rrjkdiff
      bP(iz) = 1.d0-P(iz)      
      enddo
      
c      rewind 67
c      z = closew-0.5d0*dz
c      do i = 1,imitt+ism
c      z = z+dz
c      write(67,'(3f18.9)') z,ebelam(i),P(i)   
c      enddo
c      stop

      
      nAB = 1
      do  245 iz = istp1,inw
      sume = 0.d0
      do 345 jz = istp1,iz+ism-1
 345  sume = ebelam(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*ebelam(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 245  c(iz,nmon-1) = tuu
  
      do 445 iz = inw+1,imitt
      sume = 0.5d0*ebelam(iz-ism)
      do 545 jz = iz-ism+1,iz+ism-1
 545  sume = ebelam(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*ebelam(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 445  c(iz,nmon-1) = tuu
      jz = imitt+1
      do 4445 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
 4445 c(iz,nmon-1) = tuu

      k = nmon-1
      do 745 mmm = 2,nmon-2
      k = k-1
      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do  845 iz = istp1,inw
      sume = 0.d0
      do 945 jz = istp1,iz+ism-1
 945  sume = chA(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
 845  c(iz,k) = tuu
      do 1045 iz = inw+1,imitt
      sume = 0.5d0*chA(iz-ism)
      do 1145 jz = iz-ism+1,iz+ism-1
 1145 sume = chA(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
 1045 c(iz,k) = tuu 
      jz = imitt+1
      do 7045 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chB(jz)
      chB(iz) = tuu
 7045 c(iz,k) = tuu
      else
      do  895 iz = istp1,inw
      sume = 0.d0
      do 995 jz = istp1,iz+ism-1
 995  sume = chB(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 895  c(iz,k) = tuu
      do 9045 iz = inw+1,imitt
      sume = 0.5d0*chB(iz-ism)
      do 9145 jz = iz-ism+1,iz+ism-1
 9145 sume = chB(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 9045 c(iz,k) = tuu 
      jz = imitt+1
      do 9345 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
 9345 c(iz,k) = tuu
      endif

 745  continue

      nAB = nAB+1
      if (mod(nAB,2).eq.0) then
      do  1245 iz = istp1,inw
      sume = 0.d0
      do 1345 jz = istp1,iz+ism-1
 1345 sume = chA(jz)+sume
 1245 c(iz,1) = ebelam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      do 1445 iz = inw+1,imitt
      sume = 0.5d0*chA(iz-ism)
      do 1545 jz = iz-ism+1,iz+ism-1
 1545 sume = chA(jz)+sume
 1445 c(iz,1) = ebelam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      else
      do  6245 iz = istp1,inw
      sume = 0.d0
      do 6345 jz = istp1,iz+ism-1
 6345 sume = chB(jz)+sume
 6245 c(iz,1) = ebelam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      do 6445 iz = inw+1,imitt
      sume = 0.5d0*chB(iz-ism)
      do 6545 jz = iz-ism+1,iz+ism-1
 6545 sume = chB(jz)+sume
 6445 c(iz,1) = ebelam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      endif

      z = -0.5d0*dz
      do 9 i = istp1s,imitt-ism
      z = z+dz
      dumsum = 0.d0 
      do 10 k = 2,nmon-1
 10   dumsum = c(i,k)*c(i,nmon+1-k)+dumsum
      tfem = 2.d0*c(i,1)
      tfdm = dumsum+tfem

      if (tfdm.lt.0.d0) tfdm = ddtol
      if (tfdm.gt.1.d0) tfdm = 1.d0-ddtol

      fdm = fdmon(i)
      ddiff = abs(tfdm-fdm)/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdmon(i) = fdm*dmm+tdmm*tfdm
      fdA(i) = P(i)*fdmon(i)
      fdB(i) = bP(i)*fdmon(i)      
 9    continue

      webdA = fdA(imitt-ism)
      webdB = fdB(imitt-ism)
      webdm = fdmon(imitt-ism)
      weP = P(imitt-ism)
      webP = bP(imitt-ism)
      do iz = imitt-ism+1,imitt+ism
      rkk = dfloat(iz-ist)
      rjj = dfloat(ifin-iz)
      fdA(iz) = (rjj*webdA+rkk*bdA)*rrjkdiff
      fdB(iz) = (rjj*webdB+rkk*bdB)*rrjkdiff
      fdmon(iz) = (rjj*webdm+rkk*bdm)*rrjkdiff
      P(iz) = (rjj*weP+rkk*Pb)*rrjkdiff
      bP(iz) = 1.d0-P(iz)
      enddo
      
      goto 100
 200  continue

      aaaFreen = 0.d0
      a4Freen = 0.d0
      a4bFreen = 0.d0
      Freen = 0.d0
      do i = istp1s,imitt
      fdm = fdmon(i)
      uA = 0.d0
      uB = 0.d0
c      uP = 0.d0
c      ubP = 0.d0
c      um = 0.d0

cc     do j = istp1,islut
c     do j = istp1,maxel
      do j = istp1,imitt+10*ism      
      iii = iabs(j-i)
      uA = fdA(j)*Phi(iii)+uA
      uB = fdB(j)*Phi(iii)+uB
c      uP = P(j)*Phi(iii)+uP
c      ubP = bP(j)*Phi(iii)+ubP
c      um = fdmon(j)*Phi(iii)+um
      enddo
      uA = uA+VAlrc(i)
      uB = uB+VBlrc(i)
c      uS = uS+VSlrc(i)
      
      fds = 1.d0-fdm
   
c      aaaFreen=aaaFreen+(fdA(i)*uB)*eAB+0.5d0*fdA(i)*uA*eAA+
c     *0.5d0*uB*eBB*fdB(i)+
c     *fdm*(dlog(ebelam(i))-rrnmon)+
c     *uA*(1.d0-fdm)*eAS+uB*(1.d0-fdm)*eBS+ 
c     *VexA(i)*fdA(i)+VexB(i)*fdB(i)+VexS(i)*fds+
c     *(1.d0-fdm)*dlog(1.d0-fdm)-1.d0+fdm+
c     *fdm*(P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB))

c      a4Freen = a4Freen+fdA(i)*uB*eAB+0.5d0*fdA(i)*uA*eAA+
c     *0.5d0*uB*eBB*fdB(i)-fdm*rrnmon+
c     *fdm*(P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB))+
c     *VexA(i)*fdA(i)+VexB(i)*fdB(i)+VexS(i)*fds+
c     *uA*(1.d0-fdm)*eAS+uB*(1.d0-fdm)*eBS+
c     *(1.d0-fdm)*dlog(1.d0-fdm)-1.d0+fdm+
c     *fdm*dlog(ebelam(i))

      aaaFreen = aaaFreen+fdm*(dlog(ebelam(i))-emscale)-fdm*rrnmon+fdm+
     *(bdtot-fdm)*dlog(bdtot-fdm)
     *+fdA(i)*uB*eAB+0.5d0*fdA(i)*uA*eAA+
     *0.5d0*uB*eBB*fdB(i)+
     *uA*(1.d0-fdm)*eAS+uB*(1.d0-fdm)*eBS+
     *fdm*(P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB))+
     *VexA(i)*fdA(i)+VexB(i)*fdB(i)+VexS(i)*fds

      a4Freen = a4Freen+fdA(i)*uB*eAB+0.5d0*fdA(i)*uA*eAA+
     *0.5d0*uB*eBB*fdB(i)-fdm*rrnmon+
     *uA*(1.d0-fdm)*eAS+uB*(1.d0-fdm)*eBS+
     *fdm*(dlog(ebelam(i))-emscale)+fdm+
     *(bdtot-fdm)*dlog(bdtot-fdm)+
     *fdm*(P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB))+
     *VexA(i)*fdA(i)+VexB(i)*fdB(i)+VexS(i)*fds
      
      a4bFreen = a4bFreen+fdA(i)*uB*eAB+0.5d0*fdA(i)*uA*eAA+
     *0.5d0*uB*eBB*fdB(i)-fdm*rrnmon+
     *fdm*(P(i)*dlog(P(i)/gA)+bP(i)*dlog(bP(i)/gB))+
     *VexA(i)*fdA(i)+VexB(i)*fdB(i)+VexS(i)*fds+
     *uA*(1.d0-fdm)*eAS+uB*(1.d0-fdm)*eBS+
     *(1.d0-fdm)*dlog(1.d0-fdm)-1.d0+fdm+
     *fdm*dlog(ebelam(i))-(1.d0-fdm)*chemps
      enddo
      aaaFreen = aaaFreen*dz
      a4Freen = a4Freen*dz
      a4bFreen = a4bFreen*dz
c      write(*,*) 'aaaFreen = ',aaaFreen
c      write(*,*) 'a4Freen = ',a4Freen
c      write(*,*) 'a4bFreen = ',a4bFreen
      Freen = a4bFreen
      write(*,*) 
      write(*,*) 'grand potential:'   
      write(*,*) Freen
      write(*,*) a4Freen
      write(*,*) aaaFreen            
      write(*,*) 
c      write(*,*) 'total free energy:' 
c      write(*,*) Freen+pressure*h  
      rewind ifc
      rewind ipd
c      Fs2m = 0.d0
      avfdm = 0.d0
      adsfdm = 0.d0
      adsfdA = 0.d0
      adsfdB = 0.d0
      adsfds = 0.d0
      niz = 0
      Pmacross = 0.d0
      z = sclosew-0.5d0*dz
      do 61 iz = istp1s,isluts
      niz = niz+1
      z = z+dz
c      Fs2m = fdA(iz)*dVex2A(iz)+fdB(iz)*dVex2B(iz)+Fs2m
      if (iz.lt.(imitt+1)) then
      write(ifc,'(1f14.7,3e21.9)') z,fdmon(iz),fdA(iz),fdB(iz)
      write(ipd,'(3e25.14)') z,P(iz),bP(iz)
      adsfdm = (fdmon(iz)-bdm)*dz+adsfdm
      adsfdA = (fdA(iz)-bdm*Pb)*dz+adsfdA
      adsfdB = (fdB(iz)-bdm*bPb)*dz+adsfdB
      adsfds = ((1.d0-fdmon(iz))-bds)*dz+adsfds
      endif
 61   avfdm = fdmon(iz)+avfdm
c      fp1S = fdmon(istp1)
c      fn1S = fdmon(istp1+2)
c      c0Skv = fdmon(istp1+1)
c      c1Skv = (fp1S-fn1S)*0.5d0
c      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
cC      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
c      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
c      if (fwcmq.lt.0.d0) fwcmq = 0.d0
c      write(*,*) 'monomer contact density - quad. extr.:'   
c      write(*,*) fwcmq
c      Fs2m = -Fs2m*dz
c      write(*,*) 'Fs2m = ',Fs2m
c      write(*,*) 'Fs2m+fwcmq: '
c      write(*,*) Fs2m+fwcmq
c      Pnet = Fs2m+fwcmq-pressure
c      write(*,*) 'Pnet = ',Pnet
 9998 continue
      write(*,*) 'ddmax,niter = ',ddmax,niter
      write(*,*) 
      fp1S = fdmon(imitt+1)
      fn1S = fdmon(imitt+3)
      c0Skv = fdmon(imitt+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'monomer mid plane density - quad. extr. = '
      write(*,'(1e25.14)')  fmm
      write(*,*) 
      write(*,*) 'adsfdm = ',adsfdm
      write(*,*) 'adsfdA = ',adsfdA
      write(*,*) 'adsfdB = ',adsfdB
      write(*,*) 'adsfds = ',adsfds
 9999 continue
      STOP
      END
