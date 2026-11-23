!> Export D(g,k) = <psi(Rk)|g|psi(k)> and the target orbital representation
!> in the order consumed by Wannier90 sitesym_read.
 SUBROUTINE MLWF_WANNIER90_DMN(MLWF,W,P,CQIJ,T_INFO,LATT_CUR,SYMM,IO)
   USE base, ONLY : symmetry, in_struct
   USE constant, ONLY : TPI
   USE pseudo, ONLY : potcar
   USE poscar, ONLY : type_info
   USE lattice, ONLY : latt
   USE wave_high, ONLY : wavespin, wavedes1, wavefuna, NI_GLOBAL, &
                         SETWDES, ELEMENTS, NEWWAVA, DELWAVA
   USE full_kpoints, ONLY : KPOINTS_FULL, KPOINT_IN_FULL_GRID
   USE pead, ONLY : WDES_FULL_PEAD, PEAD_WA_ROTATE
   USE wave_rotate, ONLY : wave_rotator, wave_rotator_init, wave_rotator_free, &
                           wave_rotator_rotate_pw, wave_rotator_rotate_proj, rotate_vkpt_irot
   USE choleski, ONLY : CALC_PAW_OVERLAP
   USE tutor, ONLY : vtutor
   USE string, ONLY : str

   TYPE(MLWF_type),INTENT(IN) :: MLWF
   TYPE(wavespin),INTENT(IN) :: W
   TYPE(potcar),INTENT(IN) :: P(:)
   OVERLAP,INTENT(IN) :: CQIJ(:,:,:,:)
   TYPE(type_info),INTENT(IN) :: T_INFO
   TYPE(latt),INTENT(IN) :: LATT_CUR
   TYPE(symmetry),INTENT(IN) :: SYMM
   TYPE(in_struct),INTENT(IN) :: IO

#ifndef gammareal
   INTEGER,PARAMETER :: IDENTITY(3,3)=RESHAPE([1,0,0,0,1,0,0,0,1],[3,3])
   REAL(q),PARAMETER :: GRID_TOL=5E-6_q
   TYPE(wave_rotator) :: ROTATOR
   TYPE(wavespin) :: WP
   TYPE(wavedes1),TARGET :: WDES_SRC, WDES_DST
   TYPE(wavefuna) :: BRA,KET,ROTATED
   INTEGER,ALLOCATABLE :: K_TO_IR(:),IR_TO_K(:),K_IMAGES(:,:),ISYMOP(:),ITRANS(:),BANDS(:)
   REAL(q),ALLOCATABLE :: RSHIFT(:,:,:,:),TRANS(:,:)
   COMPLEX(q),ALLOCATABLE :: ORBITAL_ROT(:,:,:),D_WANN(:,:),S(:,:),D_BAND(:,:)
   INTEGER :: NB,NK,NW,NSYM,MAXOPS,NIR,IS,IT,IG,IR,IK,JK,ISP,IU,NI,NP,NPRO,SPINOR,I,IOS
   REAL(q) :: GK(3),DK(3),ANG,MAX_LEAKAGE
   COMPLEX(q) :: PHASE
   CHARACTER(LEN=:),ALLOCATABLE :: FILENAME
   LOGICAL :: IO_NODE

   IF (SYMM%ISYM<=0) CALL vtutor%error('LWRITE_DMN requires active VASP symmetry (ISYM > 0).')
   IF (W%WDES%LSPIRAL) CALL vtutor%error('LWRITE_DMN does not support spin spirals.')
#if defined(ACC_OFFLOAD) || defined(OMP_OFFLOAD)
   CALL vtutor%error('LWRITE_DMN currently requires a CPU build.')
#endif
#ifdef oldsym
   IF (MLWF%NRSPINORS==2) CALL vtutor%error('Spinor DMN needs the current VASP symmetry implementation.')
#endif
   NK=MLWF%NKPTS
   NW=MLWF%NUM_WANN
   NB=MLWF%NB_TOT-MLWF%NEXCLB
   IF (ANY(MLWF%mp_grid<=0).OR.PRODUCT(MLWF%mp_grid)/=NK) &
      CALL vtutor%error('LWRITE_DMN requires a full uniform Wannier90 k mesh.')

   MAXOPS=SYMM%NROT*SYMM%NPCELL
   ALLOCATE(ISYMOP(MAXOPS),ITRANS(MAXOPS),TRANS(3,MAXOPS))
   NSYM=0
   ! Keep the unitary operations within each collinear spin channel.
   DO IT=1,SYMM%NPCELL
      DO IS=1,SYMM%NROT
         IF (MLWF%ISPIN==2) THEN
            IF (SYMM%MAGROT(IS,IT)<0.5_q) CYCLE
         ENDIF
         NSYM=NSYM+1
         ISYMOP(NSYM)=IS
         ITRANS(NSYM)=IT
         TRANS(:,NSYM)=SYMM%GTRANS(:,IS)+SYMM%PTRANS(IT,:)
      ENDDO
   ENDDO
   IF (NSYM==0) CALL vtutor%error('LWRITE_DMN found no unitary operations.')
   IF (ANY(SYMM%ISYMOP(:,:,ISYMOP(1))/=IDENTITY).OR.ANY(ABS(TRANS(:,1))>GRID_TOL)) &
      CALL vtutor%error('LWRITE_DMN requires identity as the first symmetry operation.')

   ALLOCATE(K_TO_IR(NK),IR_TO_K(NK),K_IMAGES(NSYM,NK))
   K_TO_IR=0
   NIR=0
   DO IK=1,NK
      IF (K_TO_IR(IK)/=0) CYCLE
      NIR=NIR+1
      IR_TO_K(NIR)=IK
      DO IG=1,NSYM
         GK=rotate_vkpt_irot(MLWF%kpt_latt(:,IK),SYMM%IGRPOP,ISYMOP(IG))
         DK=(GK-MLWF%kpt_latt(:,1))*MLWF%mp_grid
         IF (MAXVAL(ABS(DK-ANINT(DK)))>GRID_TOL) &
            CALL vtutor%error('LWRITE_DMN k mesh is not closed under the VASP symmetries.')
         JK=KPOINT_IN_FULL_GRID(GK,KPOINTS_FULL)
         K_IMAGES(IG,NIR)=JK
         K_TO_IR(JK)=NIR
      ENDDO
   ENDDO

   ALLOCATE(ORBITAL_ROT(NW,NW,NSYM),RSHIFT(3,NW,NW,NSYM),D_WANN(NW,NW))
   ! RSHIFT depends on the operation even when the k-point permutation is identical.
   IO_NODE=IO%IU6>=0
   WP=W
   WP%WDES=>WDES_FULL_PEAD
   CALL wave_rotator_init(ROTATOR,WP%WDES%GRID,LATT_CUR,SYMM,.FALSE.,MLWF%ISPIN)
   ROTATOR%MAGROT=1
   DO IG=1,NSYM
      IS=ISYMOP(IG)
      CALL MLWF_DMN_ORBITALS(MLWF,LATT_CUR,SYMM%ISYMOP(:,:,IS),TRANS(:,IG), &
                            ROTATOR%RSSYMOP(:,:,IS),ORBITAL_ROT(:,:,IG),IS,RSHIFT(:,:,:,IG))
   ENDDO

   BANDS=PACK([(I,I=1,MLWF%NB_TOT)],.NOT.MLWF%EXCLUDE_BAND)
   ALLOCATE(S(MLWF%NB_TOT,MLWF%NB_TOT),D_BAND(NB,NB))
   DO ISP=1,MLWF%ISPIN
      IF (IO_NODE) THEN
         FILENAME=TRIM(MLWF%SEED_NAME)
         IF (MLWF%ISPIN==2) FILENAME=FILENAME//'.'//str(ISP)
         WRITE(IO%IU6,'(A,I0,A,I0,A,I0)') ' DMN spin ',ISP,': unitary operations = ',NSYM,', irreducible k points = ',NIR
         DO IG=1,NSYM
            WRITE(IO%IU6,'(A,I5,9I4,3F16.10)') ' DMN operation',IG,SYMM%ISYMOP(:,:,ISYMOP(IG)),TRANS(:,IG)
         ENDDO
         OPEN(NEWUNIT=IU,FILE=FILENAME//'.dmn',STATUS='REPLACE',ACTION='WRITE',IOSTAT=IOS)
         IF (IOS/=0) CALL vtutor%error('LWRITE_DMN cannot open '//FILENAME//'.dmn')
         WRITE(IU,'(A)') 'File generated by VASP: '//MLWF%SZNAM1
         WRITE(IU,'(4I12)') NB,NSYM,NIR,NK
         WRITE(IU,'(/,10I12)') K_TO_IR
         WRITE(IU,'(/,10I12)') IR_TO_K(:NIR)
         DO IR=1,NIR
            WRITE(IU,'(/,10I12)') K_IMAGES(:,IR)
         ENDDO

         DO IR=1,NIR
            IK=IR_TO_K(IR)
            DO IG=1,NSYM
               IS=ISYMOP(IG)
               JK=K_IMAGES(IG,IR)
               GK=MLWF%kpt_latt(:,JK)
               DO NP=1,NW
                  DO NI=1,NW
                     ANG=-TPI*DOT_PRODUCT(GK,RSHIFT(:,NI,NP,IG))
                     D_WANN(NI,NP)=ORBITAL_ROT(NI,NP,IG)*CMPLX(COS(ANG),SIN(ANG),q)
                  ENDDO
               ENDDO
               WRITE(IU,*)
               WRITE(IU,'(" (",ES23.15E3,",",ES23.15E3,")")') D_WANN
            ENDDO
         ENDDO
      ENDIF

      MAX_LEAKAGE=0
      DO IR=1,NIR
         IK=IR_TO_K(IR)
         CALL SETWDES(WP%WDES,WDES_SRC,IK)
         KET=ELEMENTS(WP,WDES_SRC,ISP)
         IF (IK/=MLWF%NEQUIV(IK)) THEN
            CALL NEWWAVA(KET,WDES_SRC,WDES_SRC%NBANDS)
            CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,KET)
         ENDIF
         DO IG=1,NSYM
            IS=ISYMOP(IG)
            IT=ITRANS(IG)
            JK=K_IMAGES(IG,IR)
            CALL SETWDES(WP%WDES,WDES_DST,JK)
            BRA=ELEMENTS(WP,WDES_DST,ISP)
            IF (JK/=MLWF%NEQUIV(JK)) THEN
               CALL NEWWAVA(BRA,WDES_DST,WDES_DST%NBANDS)
               CALL PEAD_WA_ROTATE(WP,P,LATT_CUR,ISP,BRA)
            ENDIF
            CALL NEWWAVA(ROTATED,WDES_DST,WDES_DST%NBANDS)
            ROTATOR%GTRANS(:,IS)=TRANS(:,IG)
            ROTATOR%ROTMAP(:,IS)=SYMM%ROTMAP(:,IS,IT)
            CALL wave_rotator_rotate_pw(ROTATOR,KET,ROTATED,IS,LATT_CUR,FULL_BLOCH=.TRUE.)
            IF (WDES_DST%LOVERL) THEN
               CALL wave_rotator_rotate_proj(ROTATOR,KET,KET%CPROJ,ROTATED,ROTATED%CPROJ,IS,LATT_CUR,P)
               GK=rotate_vkpt_irot(MLWF%kpt_latt(:,IK),SYMM%IGRPOP,IS)
               DK=REAL(NINT(GK-MLWF%kpt_latt(:,JK)),q)
               ! PAW projector coefficients carry exp(i G.r_atom), whereas the
               ! rotated plane waves carry exp(-i (Rk+RG).t). Fold both to k'.
               DO NI=1,WDES_DST%NIONS
                  NP=NI_GLOBAL(NI,WDES_DST%COMM_INB)
                  ANG=TPI*(DOT_PRODUCT(DK,T_INFO%POSION(:,NP))-DOT_PRODUCT(GK,TRANS(:,IG)))
                  PHASE=CMPLX(COS(ANG),SIN(ANG),q)
                  NPRO=WDES_DST%LMMAX(WDES_DST%ITYP(NI))
                  DO SPINOR=0,MLWF%NRSPINORS-1
                     I=WDES_DST%LMBASE(NI)+SPINOR*(WDES_DST%NPRO/MLWF%NRSPINORS)
                     ROTATED%CPROJ(I+1:I+NPRO,:)=PHASE*ROTATED%CPROJ(I+1:I+NPRO,:)
                  ENDDO
               ENDDO
            ENDIF
            CALL CALC_PAW_OVERLAP(BRA,ROTATED,CQIJ,S)
            D_BAND=S(BANDS,BANDS)
            IF (IO_NODE) THEN
               WRITE(IU,*)
               WRITE(IU,'(" (",ES23.15E3,",",ES23.15E3,")")') D_BAND
               MAX_LEAKAGE=MAX(MAX_LEAKAGE,MAXVAL(ABS(1._q-SUM(ABS(D_BAND)**2,DIM=1))))
            ENDIF
            CALL DELWAVA(ROTATED)
            IF (JK/=MLWF%NEQUIV(JK)) CALL DELWAVA(BRA)
         ENDDO
         IF (IK/=MLWF%NEQUIV(IK)) CALL DELWAVA(KET)
      ENDDO
      IF (IO_NODE) THEN
         CLOSE(IU)
         WRITE(IO%IU6,'(A,ES12.4)') ' DMN maximum band norm defect (including band truncation): ',MAX_LEAKAGE
      ENDIF
   ENDDO
   CALL wave_rotator_free(ROTATOR)
#else
   CALL vtutor%error('LWRITE_DMN requires the complex plane-wave build (vasp_std or vasp_ncl).')
#endif
 END SUBROUTINE MLWF_WANNIER90_DMN

!> Representation of g={R|t} in the angular and spin basis of the AMN projections.
 SUBROUTINE MLWF_DMN_ORBITALS(MLWF,LATT_CUR,ROT,TRANS,SPIN_ROT,D,IS,RSHIFT)
   USE asa, ONLY : SETUP_SYM_LL
   USE lattice, ONLY : latt
   USE locproj, ONLY : WANNIER90_ORBITAL_DEFINITIONS, SETROTYLM
   USE tutor, ONLY : vtutor
   USE string, ONLY : str

   TYPE(MLWF_type),INTENT(IN) :: MLWF
   TYPE(latt),INTENT(IN) :: LATT_CUR
   INTEGER,INTENT(IN) :: ROT(3,3),IS
   REAL(q),INTENT(IN) :: TRANS(3)
   COMPLEX(q),INTENT(IN) :: SPIN_ROT(2,2)
   COMPLEX(q),INTENT(OUT) :: D(:,:)
   REAL(q),INTENT(OUT) :: RSHIFT(:,:,:)

   INTEGER,PARAMETER :: LMAX=3,LM=(LMAX+1)**2
   REAL(q) :: SL(2*LMAX+1,2*LMAX+1,0:LMAX),C(LM,MLWF%NUM_WANN),RC(LM)
   GDEF :: AXIS_ROT(LM,LM)
   COMPLEX(q) :: SPIN(2,MLWF%NUM_WANN),U(2,2),CHECK(SIZE(D,1),SIZE(D,2))
   REAL(q) :: R(3,MLWF%NUM_WANN),RROT(3),DR(3),QAXIS(3),ALPHA,BETA
   INTEGER :: I,J,L,M1,M2,NS
   EXTERNAL EULER

   NS=MLWF%NRSPINORS
   DO I=1,MLWF%NUM_WANN
      ASSOCIATE(PROJ=>MLWF%LPRJ_functions(I))
         R(:,I)=PROJ%R
         CALL WANNIER90_ORBITAL_DEFINITIONS(PROJ%L,PROJ%M,C(:,I))
         CALL SETROTYLM(PROJ%PROJ_X,PROJ%PROJ_Z,LMAX,AXIS_ROT)
         C(:,I)=MATMUL(TRANSPOSE(REAL(AXIS_ROT,q)),C(:,I))
         SPIN(:,I)=[CMPLX(1._q,0._q,q),CMPLX(0._q,0._q,q)]
         IF (NS==2) THEN
            QAXIS=PROJ%SPIN_QAXIS
            CALL EULER(QAXIS,ALPHA,BETA)
            U(:,1)=[COS(BETA/2)*EXP(CMPLX(0._q,-ALPHA/2,q)),SIN(BETA/2)*EXP(CMPLX(0._q,ALPHA/2,q))]
            U(:,2)=[-SIN(BETA/2)*EXP(CMPLX(0._q,-ALPHA/2,q)),COS(BETA/2)*EXP(CMPLX(0._q,ALPHA/2,q))]
            SPIN(:,I)=U(:,PROJ%SPINOR)
         ENDIF
      END ASSOCIATE
   ENDDO

   CALL SETUP_SYM_LL(2*LMAX+1,LMAX,ROT,SL,LATT_CUR%A,LATT_CUR%B)
   D=0
   RSHIFT=0
   DO I=1,MLWF%NUM_WANN
      RROT=MATMUL(R(:,I),ROT)+TRANS
      DO L=0,LMAX
         M1=L**2+1
         M2=(L+1)**2
         RC(M1:M2)=MATMUL(SL(1:2*L+1,1:2*L+1,L),C(M1:M2,I))
      ENDDO
      IF (NS==2) U=TRANSPOSE(CONJG(SPIN_ROT))
      DO J=1,MLWF%NUM_WANN
         DR=RROT-R(:,J)
         IF (MAXVAL(ABS(DR-ANINT(DR)))>5E-5_q) CYCLE
         IF (MLWF%LPRJ_functions(I)%N/=MLWF%LPRJ_functions(J)%N) CYCLE
         IF (ABS(MLWF%LPRJ_functions(I)%ZA-MLWF%LPRJ_functions(J)%ZA)>1E-8_q) CYCLE
         D(J,I)=DOT_PRODUCT(C(:,J),RC)
         IF (NS==2) D(J,I)=D(J,I)*DOT_PRODUCT(SPIN(:,J),MATMUL(U,SPIN(:,I)))
         RSHIFT(:,J,I)=ANINT(DR)
      ENDDO
   ENDDO
   CHECK=MATMUL(TRANSPOSE(CONJG(D)),D)
   DO I=1,SIZE(D,1)
      CHECK(I,I)=CHECK(I,I)-1
   ENDDO
   IF (MAXVAL(ABS(CHECK))>1E-5_q) &
      CALL vtutor%error('LWRITE_DMN projections are not an orthonormal closed set under VASP operation '//str(IS))
 END SUBROUTINE MLWF_DMN_ORBITALS
