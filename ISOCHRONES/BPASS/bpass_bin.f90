PROGRAM BPASS_BIN

  ! routine to convert BPASS ascii spectral files to binary
  ! also smooths and down-samples the wavelength array and
  ! renormalizes the SSPs, whose native units are Lsun/A per 1E6 Msun stars
  ! must be compiled with BPASS compiler flag

  !To Do: 1) smooth the spectra to some common resolution
  !       2) extrapolate beyond 10um
  
  USE sps_vars; USE sps_utils
  IMPLICIT NONE
  INTEGER  :: z,dumi1,i,j,stat
  INTEGER, PARAMETER :: nspeci=100000
  REAL(SP) :: dumr1,d2,d3
  CHARACTER(3), DIMENSION(nz) :: zstype
  REAL(SP), DIMENSION(nspeci,nt,nz) :: ssp_in
  REAL(SP), DIMENSION(nt,nz)  :: mass_in
  REAL(SP), DIMENSION(nt)     :: agei
  REAL(SP), DIMENSION(nspeci) :: lami
  REAL(SP), DIMENSION(nspec)  :: lamo
  
  !----------------------------------------------------------------!

  !CALL GETENV('SPS_HOME',SPS_HOME)
  SPS_HOME='/Users/cconroy/sps/'
  
  zstype = (/'em4','001','002','003','004','006','008','010',&
       '014','020','030','040'/)
  
  DO z=1,nz
     
     OPEN(92,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/v2.2/'&
          //'spectra-bin-imf135all_100.z'//zstype(z)//'.dat',&
          FORM='FORMATTED',STATUS='OLD',ACTION='READ')
     DO i=1,nspeci
        READ(92,*) lami(i), ssp_in(i,:,z)
     ENDDO
     CLOSE(92)

     OPEN(92,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/v2.2/'&
          //'starmass-bin-imf135all_100.z'//zstype(z)//'.dat',&
          FORM='FORMATTED',STATUS='OLD',ACTION='READ')
      DO i=1,nt
         READ(92,*) agei(i), mass_in(i,z),d3
     ENDDO
     CLOSE(92)
     
  ENDDO

  mass_in = mass_in/1E6

  !new wavelength grid
  !2x undersampled at <1um; 10x undersampling to 10um then log sampling to 1E8)
  j=1
  DO i=1,10000,2
     lamo(j) = lami(i)
     j=j+1
  ENDDO
  DO i=10001,nspeci,10
     lamo(j) = lami(i)
     j=j+1
  ENDDO
  DO i=j,nspec
     lamo(i) = 10**((i-j+1.)/1000.*(8-LOG10(lamo(j-1))) + LOG10(lamo(j-1)))
  ENDDO


  bpass_spec_ssp=0.
  DO z=1,nz
     DO i=1,nt
        !convert to FSPS units (Lsun/Hz per 1 Msun of mass formed)
        ssp_in(:,i,z) = ssp_in(:,i,z)/1E6*lami**2/clight
        !interpolate to new wavelength grid
        bpass_spec_ssp(1:j-1,i,z) =  MAX(linterparr(lami,ssp_in(:,i,z),&
             lamo(1:j-1)),tiny_number)
        !BB extrapolation
        bpass_spec_ssp(j:,i,z) = bpass_spec_ssp(j-1,i,z)*&
             (lamo(j:)/lamo(j-1))**(-2)
     ENDDO
  ENDDO
  
  
  !write wavelength array
  OPEN(91,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass.lambda',&
       STATUS='REPLACE',iostat=stat,ACTION='WRITE')
  DO i=1,nspec
     WRITE(91,*) lamo(i)
  ENDDO
  CLOSE(91)

  !write mass array
  OPEN(91,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass.mass',&
       STATUS='REPLACE',iostat=stat,ACTION='WRITE')
  DO i=1,nt
     WRITE(91,'(20F9.3)') agei(i), mass_in(i,:)
  ENDDO
  CLOSE(91)

  !write BPASS SSP binary file
  OPEN(93,FILE=TRIM(SPS_HOME)//'/ISOCHRONES/BPASS/bpass_v2.2_salpeter100.ssp.bin',&
       FORM='UNFORMATTED',STATUS='REPLACE',access='direct',&
       recl=nspec*nt*nz*8)
  WRITE(93,rec=1) bpass_spec_ssp
  CLOSE(93)
     
  
END PROGRAM BPASS_BIN
