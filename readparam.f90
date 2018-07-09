
c ===================================================================== 
      subroutine readlist(INFILE,OUTFIL,PLTFIL,DRVFIL,DPLTFL,PTBFIL,
     $                    LSTFIL,LISTF,ITER,NLF)
      include 'REG6ST.FOR'
      CHARACTER*64 INFILE,OUTFIL,PLTFIL,DRVFIL,DPLTFL,PTBFIL
      CHARACTER*64 PRE1,PRE2,POST  
      CHARACTER*64 LSTFIL
      CHARACTER*64 CONV 
      INTEGER NLF,NITER,ITER,LISTF

      NAMELIST/SETN/INFILE,OUTFIL,PLTFIL,DRVFIL,DPLTFL,PTBFIL,
     $
     $        DERIVF,GAMMAF,ARF,PTAIF,TTAIF,VINF,ALPINF,ALPOTF,
     $        MONEF,MENDF,MRSPF,ZMSPF,RMSPF,BESPF,MRWF,OMRWF,ETARWF,
     $        DPTRWF,BLKRWF,VRWF,ALPRWF,SVCRCF,MREFF,PHREFF,MLEF,PHILEF,
     $        RAF,RBF,TAUF,MTEF,PHITEF,RTEF,PLEF,PUPSF,FLEF,FTEF,BTUPSF,
     $        BTUPEF,MUPF,PHIUPF,BTLWSF,BTLWEF,MLWF,PHILWF,BLDFCF,MBLDF,
     $        RPHIF,HEPSCF,SECTRF,
     $
     $        VINI,ALPINI,ALPOTI,MRSPI,ZMSPI,RMSPI,
     $        BESPI,MRWI,OMRWI,ETARWI,DPTRWI,BLKRWI,VRWI,ALPRWI,ETARWJ,
     $        DPTRWJ,BLKRWJ,VRWJ,ALPRWJ,SVCRCK,MREFI,MREFK,PHREFI,PHREFK
     $        ,MLEK,PHILEK,RAK,RBK,TAUK,MTEK,PHITEK,RTEK,PLEK,PUPSK,FLEK
     $        ,FTEK,BTUPSK,BTUPEK,MUPI,MUPK,PHIUPI,PHIUPK,BTLWSK,BTLWEK,
     $        MLWI,MLWK,PHILWI,PHILWK,BLDFCK,MBLDI,MBLDK,RPHII,RPHIK
      data PRE1,PRE2,POST/'./dvlists/dvlist','dvlist','.in'/
      call filename1(PRE1,iter,POST,CONV)
      LSTFIL=CONV
      OPEN(UNIT=NLF,FILE=LSTFIL,STATUS='OLD',ERR=2)
      READ(NLF,NML=SETN)
      LISTF=1
      print*,''
      print*,'FOUND LIST FILE ',LSTFIL
      close (unit=NLF)
      return
    2 call filename1(PRE2,iter,POST,CONV)
      LSTFIL=CONV
      OPEN(UNIT=NLF,FILE=LSTFIL,STATUS='OLD',ERR=3)
      READ(NLF,NML=SETN)
      LISTF=1
      print*,''
      print*,'FOUND LIST FILE ',LSTFIL
      close (unit=NLF)
      return
    3 print*,''
      print*,'CANNOT FIND LIST FILE ',LSTFIL
      INFILE=''
      LISTF=0
      NITER=1
      return 
      end
c ===================================================================== 
