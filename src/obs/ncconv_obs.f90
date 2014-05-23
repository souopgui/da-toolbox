!> \file ncconv_obs.f90
!> Netcdf interface for conventional observations

!> Module defining the data structures and subroutines to sample and save observations to netcdf file
MODULE ncconv_obs
  USE conv_obs
  USE nctools
  USE general_tools
  USE debug_tools
  USE randmod
IMPLICIT NONE


  !> \brief User defined type for attributes of observations
  !<
  TYPE obs_structure_components
    CHARACTER(LEN=ip_cnl) ::&
      obs_fName  = "obs_fName",&
      ogap_fName = "ogap_fName",&
      aa_title   = "title",&
      aa_history = "history",&
      i_max_nobs = "i_max_nobs",&
      i_ndim = "i_ndim",&
      i_nobs   = "i_nobs",&
      r_dt     = "delta_t" ,&
      l_icoord_allocated = "l_icoord_allocated",&
      l_rcoord_allocated = "l_rcoord_allocated"

    CHARACTER(LEN=ip_cnl) :: ra_data = "data"
    CHARACTER(LEN=ip_cnl) ::&
      r_date   = "t",&
      i_ts     = "timeStep",&
      ra_obs    = "ra_obs",&
      ra_sigma   = "ra_sigma",&
      ra_Rm1    = "ra_Rm1",&
      ra_obsgap   = "ra_obsgap"
    CHARACTER(LEN=ip_cnl) ::&
      r_dateu  = "s",&
      i_tsu    = "N/A",&
      ra_obsu   = "N/A",&
      ra_sigmau  = "N/A",&
      ra_Rm1u   = "N/A",&
      ra_obsgapu  = "N/A"
    CHARACTER(LEN=ip_lnl) ::&
      r_dateln = "time",&
      i_tsln   = "time step ordinal number",&
      ra_obsln  = "observation data",&
      ra_sigmaln = "standard deviation of observation error",&
      ra_Rm1ln  = "inverse covariance (diag matrix)",&
      ra_obsgapln = "Diference between obs and model output"

    !wavelet related variables
    CHARACTER(LEN=ip_cnl) :: i_obs_level = "Obs_level"
    CHARACTER(LEN=ip_cnl) :: ia_M_vector = "M_vector", ia_M_vectoru = "N/A"
    CHARACTER(LEN=ip_lnl) :: ia_M_vectorln = "Size of the wavelet grid at the coarsest scale"
    !discretization variables
    CHARACTER(LEN=ip_cnl) :: ia_nobs = "nObs"
    CHARACTER(LEN=ip_cnl) :: ra_dx = "obs_step"
    CHARACTER(LEN=ip_cnl) :: ra_dxu = "m"
    CHARACTER(LEN=ip_lnl) :: ra_dxln = "step between successive observation points"

    CHARACTER(LEN=ip_cnl) ::&
      ia_icoord = "icoord",&
      ra_rcoord   = "rcoord",&
      l_rcoord   = "l_rcoord",&
      l_icoord   = "l_icoord"
    CHARACTER(LEN=ip_cnl) ::&
      ia_icoordu = "N/A",&
      ra_rcoordu  = "N/A",&
      l_rcoordu  = "N/A",&
      l_icoordu  = "l_icoord"
    CHARACTER(LEN=ip_lnl) ::&
      ia_icoordln = "Integer coordinates, indices in the discretization grid",&
      ra_rcoordln = "real coordinates in the computation domain",&
      l_rcoordln = "Status of the real coordinates",&
      l_icoordln = "Status of the integer coordinates"

  END TYPE obs_structure_components

  TYPE obsOutput
    CHARACTER(LEN=ip_snl) :: fileName
    CHARACTER(LEN=ip_snl) :: aa_title
    CHARACTER(LEN=ip_fnl) :: aa_history
    !wavelet related variables
    INTEGER :: i_obs_level = 0!< observation level, this variable changes the interpretation of integer coordinates
    INTEGER, DIMENSION(3) :: ia_M_vector = (/0, 0, 0/)!<wavelet grid at the coarsest scale (level =1)
    INTEGER, DIMENSION(3) :: ia_nobs = (/0, 0, 0/)!<Number of observations in each dimension
    REAL(dp), DIMENSION(3) :: ra_dx = (/-999, 999, 999/)

    INTEGER :: ncid      = -1
    INTEGER :: i_nobs    = -1
    INTEGER :: i_ndim    = -1
    REAL(dp):: r_dt      = -1.0d0
    INTEGER :: nb_record = -1 !<number of record in the file

    INTEGER :: ra_dataid   = -1
    INTEGER :: ra_sigmaid  = -1
    INTEGER :: ia_icoordid = -1
    INTEGER :: ra_rcoordid = -1
    INTEGER :: l_rcoordid  = -1
    INTEGER :: l_icoordid  = -1
    INTEGER :: i_tsid      = -1
    INTEGER :: r_dateid    = -1

    LOGICAL :: isOpened    = .FALSE.
    INTEGER :: i_nextRec   = -1
  END TYPE obsOutput

  INTERFACE initObsOutput
    MODULE PROCEDURE initObsOutput_with_scalars
    MODULE PROCEDURE initObsOutput_with_cop
  END INTERFACE

  !>
  INTERFACE ncObsRead
    MODULE PROCEDURE ncObsRead_from_DS
    MODULE PROCEDURE ncObsRead_snapshot
  END INTERFACE

  !>interface to read all the times steps associated with observation trajectory
  interface readObs_AllTimeStep
    module procedure readObs_AllTimeStep_DS
    module procedure readObs_AllTimeStep_fName
  end interface

  INTERFACE obs_sample_and_save
    MODULE PROCEDURE obs_sample_and_save_1D
    MODULE PROCEDURE obs_sample_and_save_2D
  END INTERFACE

  TYPE(obs_structure_components), PRIVATE, SAVE :: tm_oAtt
CONTAINS

  !> \brief initializes data structure that is used to saved conventional observations
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[in] ada_fileName file name for observation
  !! \param[in] nobs (optional) number of observation point in the physical space
  !! \param[in] ndim (optional) number of dimension of the physical space
  !! \param[in] dt (optional) time step
  !! \param[in] title (optional) title information for the file
  !! \param[in] hist (optional) history string, used to build hystory of the file
  !<
  SUBROUTINE initObsOutput_with_scalars(td_ncfile, ada_fileName, nobs, ndim, dt, title, hist)
    type(obsOutput), INTENT(INOUT) :: td_ncfile
    CHARACTER(LEN=*), INTENT(IN) :: ada_fileName
    INTEGER, INTENT(IN), OPTIONAL :: nobs, ndim
    REAL(dp), INTENT(IN), OPTIONAL:: dt
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: title, hist
    CHARACTER(LEN=ip_snl) :: ala_hist

    td_ncfile%fileName = ada_fileName
    IF( PRESENT(nobs) )THEN
      td_ncfile%i_nobs  = nobs
    ELSE
      td_ncfile%i_nobs  = -1
    END IF
    IF( PRESENT(ndim) )THEN
      td_ncfile%i_ndim  = ndim
    ELSE
      td_ncfile%i_ndim  = -1
    END IF
    IF( PRESENT(dt) )THEN
      td_ncfile%r_dt  = dt
    ELSE
      td_ncfile%r_dt  = -1.0d0
    END IF
    IF( PRESENT(title) )THEN
      td_ncfile%aa_title = title
    ELSE
      td_ncfile%aa_title = "Conventional observation for twin experiments"
    END IF
    IF( PRESENT(hist) )THEN
      ala_hist = hist
    ELSE
      ala_hist = "unnamed model for twin experiments"
    END IF
    td_ncfile%aa_history = make_ncHistory(ala_hist)
  END SUBROUTINE initObsOutput_with_scalars

  !> \brief initializes data structure that is used to saved conventional observations
  !! \param [in,out] td_ncfile data structure
  !! \param [in] td_cop data structure for conventional obs param
  !! \param [in] status (optional) status of the file (INPUT_FILE or OUTPUT_FILE), default is OUTPUT_FILE
  !<
  SUBROUTINE initObsOutput_with_cop(td_ncfile, td_cop, status)
    type(obsOutput), INTENT(INOUT)  :: td_ncfile
    TYPE(convObs_param), INTENT(IN) :: td_cop
    INTEGER, OPTIONAL, INTENT(IN)   :: status
    !local variables
    CHARACTER(len=ip_fnl) :: ala_fName
    INTEGER :: il_nobs, il_status

    IF(PRESENT(status))THEN
      il_status = status
    ELSE
      il_status = OUTPUT_FILE
    END IF

    ala_fName = make_fileName(OBS_DATA, il_status, prefix=td_cop%aa_prefix)
    il_nobs = regular_coord_count(ida_max=td_cop%ia_max, shifts=td_cop%ia_shifts, steps=td_cop%ia_steps)
    CALL initObsOutput_with_scalars(&
      td_ncfile, ala_fName,&
      nobs = il_nobs,&
      ndim = td_cop%i_ndim,&
      dt   = td_cop%r_dt,&
      title= td_cop%aa_title,&
      hist = td_cop%aa_hist &
    )
  END SUBROUTINE initObsOutput_with_cop

  SUBROUTINE printObsOutput(td_ncfile)
    type(obsOutput), INTENT(IN) :: td_ncfile
    PRINT*, 'fileName  = ', TRIM(td_ncfile%fileName)
    PRINT*, 'i_nobs    = ', td_ncfile%i_nobs
    PRINT*, 'nb_record = ', td_ncfile%nb_record
    PRINT*, 'r_dt      = ', td_ncfile%r_dt
  END SUBROUTINE printObsOutput

  !> Closes the observation trajectory file that is connected to \a td_ncfile
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !!
  !<
  SUBROUTINE ncObsClose(td_ncfile)
    type(obsOutput), INTENT(INOUT) :: td_ncfile

    CALL chkerr( nf90_close( td_ncfile%ncid), fname=td_ncfile%filename )
    td_ncfile%ncid = -1
    td_ncfile%i_nextRec = -1
    td_ncfile%isOpened = .FALSE.
  END SUBROUTINE ncObsClose

  !> Gets the number of records of observation trajectory stored in a netcdf file
  !! \param [in] fName name of the observation trajectory file
  !<
  function ncObs_nrec(fName) result(nrec)
    character(len=*), intent(in) :: fName
    !local variables
    type(ObsOutput) :: tl_oo
    integer :: nrec

    tl_oo%fileName = fName
    call ncObsOpen(tl_oo)
    nrec = tl_oo%nb_record
    call ncObsClose(tl_oo)
  end function ncObs_nrec

  !> Opens the netcdf file that stores observation trajectory
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !!
  !<
  SUBROUTINE ncObsOpen(td_ncfile)
    TYPE(obsOutput), INTENT(INOUT) :: td_ncfile
    INTEGER :: il_nDimensions, il_nVariables, il_nAttributes, il_unlimitedDimId, il_formatNum,&
      il_nobsid, il_ndimid
    character(len = nf90_max_name) :: lc_name

    CALL debug( TRIM(td_ncfile%fileName), 'In ncObsOpen, opening -- ', tag=dNETCDF )

    CALL chkerr(nf90_open(TRIM(td_ncfile%filename),NF90_NOCLOBBER,td_ncfile%ncid), fname=td_ncfile%filename)

    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%ra_data  , td_ncfile%ra_dataid   ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%ra_sigma , td_ncfile%ra_sigmaid  ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%ia_icoord, td_ncfile%ia_icoordid ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%ra_rcoord, td_ncfile%ra_rcoordid ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%l_icoord , td_ncfile%l_icoordid  ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%l_rcoord , td_ncfile%l_rcoordid  ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%r_date   , td_ncfile%r_dateid    ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_oAtt%i_ts     , td_ncfile%i_tsid      ) )

    CALL chkerr(nf90_inquire(td_ncfile%ncid, il_nDimensions, il_nVariables, il_nAttributes, &
        il_unlimitedDimId, il_formatNum))

    CALL chkerr(nf90_inquire_dimension(td_ncfile%ncid, il_unlimitedDimId, name = lc_name, len = td_ncfile%nb_record))
    CALL debug( td_ncfile%nb_record, 'In ncObsOpen, td_ncfile%nb_record = ', tag=dNETCDF )
    CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_oAtt%i_nobs, il_nobsid ) )
    CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_oAtt%i_ndim, il_ndimid ) )
    CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_nobsid, lc_name, td_ncfile%i_nobs ) )
    CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_ndimid, lc_name, td_ncfile%i_ndim ) )
    !PRINT*, "getting attributres"
    CALL readatt(td_ncfile)
    td_ncfile%isOpened = .TRUE.
    td_ncfile%i_nextRec = 1
  CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncObsOpen

  !> Reads attributes (parameters) of the observation trajectory stored in a netcdf file
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !!
  !<
  SUBROUTINE readatt( td_ncfile )
    type(ObsOutput), INTENT(INOUT) :: td_ncfile
    CALL debug("In readatt, reading attributes", tag=dNETCDF )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%r_dt       , td_ncfile%r_dt       ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%aa_title   , td_ncfile%aa_title   ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%aa_history , td_ncfile%aa_history ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%i_obs_level, td_ncfile%i_obs_level) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ia_M_vector, td_ncfile%ia_M_vector) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ia_nobs    , td_ncfile%ia_nobs    ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ra_dx      , td_ncfile%ra_dx      ) )
    CALL debug("End of attributes reading", tag=dNETCDF )
    !End of saving
  END SUBROUTINE readatt

  !> Creates a netcdf file for the observation trajectory described by \a td_ncfile
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !<
  SUBROUTINE ncObsCreate(td_ncfile)
    type(ObsOutput), INTENT(INOUT) :: td_ncfile
    INTEGER, DIMENSION(3) :: ila_dims3D
    INTEGER, DIMENSION(2) :: ila_dims2D
    INTEGER :: il_nobs, il_ndim, il_date

    CALL debug( TRIM(td_ncfile%fileName), 'In ncObsCreate, creating -- ', tag=dNETCDF )
    CALL chkerr(nf90_create(TRIM(td_ncfile%fileName),NF90_CLOBBER , td_ncfile%ncid), fname=td_ncfile%filename)
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_oAtt%i_nobs, td_ncfile%i_nobs, il_nobs))
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_oAtt%i_ndim, td_ncfile%i_ndim, il_ndim))
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_oAtt%r_date, NF90_UNLIMITED  , il_date))

    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%r_date  , NF90_DOUBLE, il_date, td_ncfile%r_dateid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%i_ts    , NF90_INT   , il_date, td_ncfile%i_tsid  ))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%l_icoord, NF90_INT   , il_date, td_ncfile%l_icoordid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%l_rcoord, NF90_INT   , il_date, td_ncfile%l_rcoordid))
    ila_dims2D(1)=il_nobs
    ila_dims2D(2)=il_date
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%ra_data , NF90_DOUBLE, ila_dims2D, td_ncfile%ra_dataid ))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%ra_sigma, NF90_DOUBLE, ila_dims2D, td_ncfile%ra_sigmaid))
    ila_dims3D(1)=il_ndim
    ila_dims3D(2)=il_nobs
    ila_dims3D(3)=il_date
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%ia_icoord, NF90_INT   , ila_dims3D, td_ncfile%ia_icoordid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_oAtt%ra_rcoord, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_rcoordid))
    CALL saveatt(td_ncfile)
    CALL chkerr(nf90_enddef(td_ncfile%ncid))
    td_ncfile%isOpened = .TRUE.
    td_ncfile%i_nextRec = 1
    CALL debug( '... ncObsCreate -> done', tag=dNETCDF )
  END SUBROUTINE ncObsCreate

  !> Saves the parameters of the observation trajectory
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !<
  SUBROUTINE saveatt(td_ncfile)
    type(ObsOutput), INTENT(INOUT) :: td_ncfile
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%r_dt       , td_ncfile%r_dt       ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%aa_title   , td_ncfile%aa_title   ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%aa_history , td_ncfile%aa_history ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%i_obs_level, td_ncfile%i_obs_level) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ia_M_vector, td_ncfile%ia_M_vector) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ia_nobs    , td_ncfile%ia_nobs    ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_oAtt%ra_dx      , td_ncfile%ra_dx      ) )
    CALL debug( "End of attributes recording", tag=dNETCDF )
    !End of saving
  END SUBROUTINE saveatt

  !> \brief Writes observation data to netcdf file
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[in,out] td_obs obs structure that contains the data to write
  !! \param[in] rec (optional) gives the ordinal number of the record to write
  !! \param[in] ogap (optional) this argument changes completely the behaviour of the subroutine. It value is not important. When it is absent, the file is supposed to contain observation data, otherwise, it contains obsgap data.
  !! When it is present, the file is supposed to contain the obsgap data, difference between observation and the state. When reading obsgap, only the obsgap field and the RM1 field are read
  !! \note the fields obs abd obsgap share the same place and name in the netcdf file. Similarly, the field Rm1 and sigma share the same place and name in the netcdf file.
  !<
  SUBROUTINE ncObsWrite(td_ncfile, td_obs, rec, ogap)
    type(ObsOutput), INTENT(INOUT) :: td_ncfile
    type(obs_structure), INTENT(IN) :: td_obs
    INTEGER, INTENT(IN), OPTIONAL  :: rec
    INTEGER, INTENT(IN), OPTIONAL  :: ogap

    INTEGER, DIMENSION(3) :: ila_count3D, ila_start3D
    INTEGER, DIMENSION(2) :: ila_count2D, ila_start2D
    INTEGER, DIMENSION(1) :: ila_start1D
    INTEGER :: il_rec, il_icoord, il_rcoord

    IF( PRESENT(rec) )THEN
      il_rec = rec
    ELSE
      il_rec = td_ncfile%i_nextRec
    END IF

    CALL debug(il_rec, 'In ncObsWrite, writing '//TRIM(td_ncfile%filename)//' - record', tag=dNETCDF )
    ila_start1D = (/il_rec/)
    ila_count2D=(/td_ncfile%i_nobs, 1 /)
    ila_start2D=(/ 1, il_rec /)
    ila_count3D= (/td_ncfile%i_ndim, td_ncfile%i_nobs, 1/)
    ila_start3D=(/ 1, 1, il_rec /)
    il_icoord = INT(td_obs%l_icoord)
    il_rcoord = INT(td_obs%l_rcoord)
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%l_icoordid, il_icoord    , start = ila_start1D))
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%l_rcoordid, il_rcoord    , start = ila_start1D))
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%r_dateid  , td_obs%r_date, start = ila_start1D))
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%i_tsid    , td_obs%i_ts  , start = ila_start1D))
    IF(td_obs%l_icoord)&
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ia_icoordid, td_obs%ia_icoord, start = ila_start3D, COUNT = ila_count3D))
    IF(td_obs%l_rcoord)&
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_rcoordid, td_obs%ra_rcoord, start = ila_start3D, COUNT = ila_count3D))

    IF( PRESENT(ogap) ) THEN
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_dataid , td_obs%ra_obsgap  , start = ila_start2D, COUNT = ila_count2D))
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_sigmaid, td_obs%ra_Rm1, start = ila_start2D, COUNT = ila_count2D))
    ELSE
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_dataid , td_obs%ra_obs  , start = ila_start2D, COUNT = ila_count2D))
      CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_sigmaid, td_obs%ra_sigma, start = ila_start2D, COUNT = ila_count2D))
    END IF
    td_ncfile%i_nextRec = il_rec + 1
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncObsWrite

  !> \brief Reads obs from netcdf file using obs output data structure
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[in,out] td_obs obs structure to contain the read data
  !! \param[in] rec (optional) gives the ordinal number of the record to read
  !! \param[in] ogap (optional) this argument changes completely the behaviour of the subroutine. It value is not important. When it is absent, the file is supposed to contain observation data and the subroutine allocates space if needed.
  !! When it is present, the file is suppose to contain the obsgap data, difference between observation and the state. When reading obsgap, only the obsgap field and the RM1 field are read. The space is suppose to be allocated.
  !!
  !! \note the fields obs and obsgap share the same place and name in the netcdf file. Similarly, the field Rm1 and sigma share the same place and name in the netcdf file.
  !<
  SUBROUTINE ncObsRead_from_DS(td_ncfile, td_obs, rec, ogap)
    type(ObsOutput), INTENT(INOUT) :: td_ncfile
    type(obs_structure), INTENT(INOUT) :: td_obs
    INTEGER, INTENT(IN), OPTIONAL  :: rec
    INTEGER, INTENT(IN), OPTIONAL  :: ogap

    INTEGER, DIMENSION(3) :: ila_count3D, ila_start3D
    INTEGER, DIMENSION(2) :: ila_count2D, ila_start2D
    INTEGER, DIMENSION(1) :: ila_start1D, ila_count1D
    !REAL(KIND=sp), DIMENSION(1)  :: rl_date
    INTEGER :: il_rec, il_icoord, il_rcoord
    LOGICAL :: ll_icoord, ll_rcoord

    IF( PRESENT(rec) )THEN
      il_rec = rec
    ELSE
      il_rec = td_ncfile%i_nextRec
    END IF
    CALL debug(il_rec, 'In ncObsRead, Reading '//TRIM(td_ncfile%filename)//' - record', tag=dNETCDF )
    ila_start1D = (/il_rec/)
    ila_count1D = (/1/)
    ila_count2D = (/td_ncfile%i_nobs, 1/)
    ila_start2D = (/ 1, il_rec /)
    ila_count3D = (/td_ncfile%i_ndim, td_ncfile%i_nobs, 1/)
    ila_start3D = (/ 1, 1, il_rec /)

    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%l_icoordid, il_icoord, start = ila_start1D ) )
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%l_rcoordid, il_rcoord, start = ila_start1D ) )
    ll_icoord = int2l(il_icoord)
    ll_rcoord = int2l(il_rcoord)
    CALL set_obsSize( td_obs, td_ncfile%i_nobs, id_coord_ndim=td_ncfile%i_ndim, ld_icoord=ll_icoord, ld_rcoord=ll_rcoord )
    !integer coordinates
    IF(td_obs%l_icoord)THEN
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ia_icoordid, td_obs%ia_icoord, start = ila_start3D, COUNT = ila_count3D))
    END IF
    !real coordinates
    IF(td_obs%l_rcoord)THEN
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_rcoordid, td_obs%ra_rcoord, start = ila_start3D, COUNT = ila_count3D))
    END IF
    CALL chkerr( nf90_get_var( td_ncfile%ncid, td_ncfile%r_dateid  , td_obs%r_date, start = ila_start1D ) )
    CALL chkerr (nf90_get_var( td_ncfile%ncid, td_ncfile%i_tsid    , td_obs%i_ts  , start = ila_start1D ) )

    IF( PRESENT(ogap) ) THEN
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_dataid , td_obs%ra_obsgap  , start = ila_start2D, COUNT = ila_count2D))
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_sigmaid, td_obs%ra_Rm1, start = ila_start2D, COUNT = ila_count2D))
    ELSE
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_dataid , td_obs%ra_obs  , start = ila_start2D, COUNT = ila_count2D))
      CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_sigmaid, td_obs%ra_sigma, start = ila_start2D, COUNT = ila_count2D))
      td_obs%ra_Rm1 = 1.0_cp/td_obs%ra_sigma**2
    END IF
    td_ncfile%i_nextRec = il_rec + 1
    td_obs%ra_dx = td_ncfile%ra_dx
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncObsRead_from_DS

  !> Read an obs snapshot from non opened netcdf file
  !! \param [in] fName name of the file
  !! \param [in,out] td_obs obs structure to contain the read data
  !! \param [in] rec gives the ordinal number of the record to read
  !! \param [in] ogap (optional) this argument changes completely the behaviour of the subroutine. \see ncObsRead
  !<
  subroutine ncObsRead_snapshot(fName, td_obs, rec, ogap)
    character(len=*), intent(in) :: fName
    type(obs_structure), INTENT(INOUT) :: td_obs
    integer, intent(in) :: rec
    INTEGER, INTENT(IN), OPTIONAL  :: ogap
    !local variables
    type(ObsOutput) :: tl_oo

    tl_oo%fileName = fName
    call ncObsOpen(tl_oo)
    call ncObsRead(tl_oo, td_obs, rec, ogap)
    call ncObsClose(tl_oo)
  end subroutine ncObsRead_snapshot


  !> Reads the dates associated with all the records of the observation trajectory
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[out] rda_date vector of date
  !!
  !<
  SUBROUTINE readObs_date(td_ncfile, rda_date)
    type(ObsOutput),INTENT(IN)      :: td_ncfile
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT)    :: rda_date

    CALL debug(TRIM(td_ncfile%filename), 'In readObs_date : reading date from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%r_dateid, rda_date))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readObs_date

  !> Reads a time step associated with a record of the observation trajectory
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[in] id_recNum Records of interest
  !! \param[out] id_ts time stept
  !!
  !<
  SUBROUTINE readObs_timeStep(td_ncfile, id_recNum, id_ts)
    type(ObsOutput),INTENT(IN):: td_ncfile
    INTEGER, INTENT(IN)         :: id_recNum
    INTEGER, INTENT(OUT)        :: id_ts

    CALL debug(TRIM(td_ncfile%filename), 'In readObs_timeStep : reading timeStep from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%i_tsid, id_ts, start=(/id_recNum/) ))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readObs_timeStep

  !> Reads the time steps associated with all the records of the observation trajectory, described by \a td_ncfile
  !! \param[in,out] td_ncfile data structure describing the observation trajectory
  !! \param[out] ida_ts vector of time steps
  !!
  !<
  SUBROUTINE readObs_AllTimeStep_DS(td_ncfile, ida_ts)
    type(ObsOutput),INTENT(IN)    :: td_ncfile
    INTEGER, DIMENSION(:), INTENT(OUT):: ida_ts

    CALL debug(TRIM(td_ncfile%filename), 'In readObs_AllTimeStep : reading timeStep from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%i_tsid, ida_ts))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readObs_AllTimeStep_DS


  !> Reads the time steps associated with all the records of the observation trajectory store in the file \a fName
  !! \param [in] fName name of the observation trajectory file
  !! \param [out] ida_ts vector of time steps
  !<
  subroutine readObs_AllTimeStep_fName(fName, ida_ts)
    character(len=*), intent(in) :: fName
    INTEGER, DIMENSION(:), INTENT(OUT):: ida_ts
    !local variables
    type(ObsOutput) :: tl_oo

    CALL debug(TRIM(fName), 'In readObs_AllTimeStep : reading timeStep from ', tag=dNETCDF)
    tl_oo%fileName = fName
    call ncObsOpen(tl_oo)
    call chkerr(nf90_get_var(tl_oo%ncid, tl_oo%i_tsid, ida_ts))
    call ncObsClose(tl_oo)
    CALL debug('... done', tag=dNETCDF)
  end subroutine readObs_AllTimeStep_fName

  !> \brief sample and save conventional observation
  !! \param[in] rda_state state of the system
  !! \param[in, out] td_obs data structure of the observation
  !! \param[in, out] td_obsOut data struction of the observation output file
  !! \param[in] td_cop data structure of the observation parameters
  !! \param[in] id_ts time step associated with observation to be saved
  !! \details rda_state can be only one of multiple variables describing the system state, that is the case that is automatically accounted for. If the user does not want to separates variables, he may have to do some work for the general case. Observation are sample from rda_state using coordinates in td_obs.
  !! This routine supposes that any dynamic array is allocated and have the right size. It is also assumed that the random number generator is initialized accordingly.
  !<
  SUBROUTINE obs_sample_and_save_1D( rda_state, td_obs, td_obsOut, td_cop, id_ts )
    TYPE(obs_structure), INTENT(IN OUT):: td_obs
    TYPE(convObs_param), INTENT(IN)    :: td_cop
    REAL(dp), DIMENSION(:), INTENT(IN) :: rda_state
    type(ObsOutput),INTENT(IN OUT)     :: td_obsOut
    INTEGER, INTENT(IN) :: id_ts
    !local variables
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rda_error
    INTEGER :: il_idx

    IF ( find( td_cop%ia_samplingTS, id_ts, il_idx ) )THEN
      IF(td_cop%l_sampleObs)THEN
        !Since there can be many observations with different error standard deviation, the randmod is initialize at every call
        ALLOCATE( rda_error( SIZE(td_obs%ra_obs) ) )
        IF(td_cop%r_sigmaR .GT. epsilon(0.0_cp))THEN
          CALL init_normal_rand(mu=0.0_cp, sigma=td_cop%r_sigmaR)
          CALL normal_random(rda_error)
          CALL debug(rda_error, "in obs_sample_and_save_1D, rda_error = ", tag=dTRACE)
          td_obs%ra_sigma = td_cop%r_sigmaR
        ELSE
          rda_error = 0.0_cp
          td_obs%ra_sigma = 1.0_cp
        END IF
        CALL subsample( rda_state, td_obs%ia_icoord(1,:), td_obs%ra_obs, error=rda_error )
        td_obs%l_icoord = .TRUE.
        td_obs%l_rcoord = .FALSE.
        td_obs%i_ts     = id_ts
        td_obs%r_date   = id_ts*td_cop%r_dt
        CALL debug( id_ts, "saving obs for time step :", tag=dTRACE )
        CALL ncObsWrite( td_obsOut, td_obs, rec=il_idx )
      END IF
    END IF
  END SUBROUTINE obs_sample_and_save_1D

  !> \brief Samples and save conventional observation fo two physical dimensions
  !! \param[in] rda_state state of the system
  !! \param[in, out] td_obs data structure of the observation
  !! \param[in, out] td_obsOut data struction of the observation output file
  !! \param[in] td_cop data structure of the observation parameters
  !! \param[in] id_ts time step associated with observation to be saved
  !! \details \a rda_state can be only one of multiple variables describing the system state, that is the case that is automatically accounted for. If the user does not want to separates variables, he has to write a wrapping subroutine. Observation are sample from \a rda_state using coordinates in td_obs.
  !!This routine supposes that any dynamic array is allocated and have the right size. It is also assumed that the random number generator is initialized accordingly.
  !<
  SUBROUTINE obs_sample_and_save_2D(rda_state, td_obs, td_obsOut, td_cop, id_ts)
    TYPE(obs_structure), INTENT(IN OUT)  :: td_obs
    TYPE(convObs_param), INTENT(IN)      :: td_cop
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_state
    type(ObsOutput),INTENT(IN OUT)       :: td_obsOut
    INTEGER, INTENT(IN) :: id_ts
    !local variables
    REAL(dp), DIMENSION(:), ALLOCATABLE :: rda_error
    INTEGER :: il_idx

    IF ( find( td_cop%ia_samplingTS, id_ts, il_idx ) )THEN
      CALL debug(SHAPE(rda_state), "In obs_sample_and_save_2D, SHAPE(rda_state) = ", tag=dTRACE)
      IF(td_cop%l_sampleObs)THEN
        !Since there can be many observations with different error standard deviation, the randmod is initialize at every call
        ALLOCATE( rda_error( SIZE(td_obs%ra_obs) ) )
        IF(td_cop%r_sigmaR .GT. epsilon(0.0_cp))THEN
          CALL init_normal_rand(mu=0.0_cp, sigma=td_cop%r_sigmaR)
          CALL normal_random(rda_error)
          td_obs%ra_sigma = td_cop%r_sigmaR
        ELSE
          rda_error = 0.0_cp
          td_obs%ra_sigma = 1.0_cp
        END IF
        !CALL debug(td_obs%ia_icoord, 'td_obs%ia_icoord = ', tag=dALLWAYS)
        !CALL dpause('In obs_sample_and_save_2D')
        CALL subsample( rda_state, td_obs%ia_icoord, td_obs%ra_obs, error=rda_error )
        td_obs%l_icoord = .TRUE.
        td_obs%l_rcoord = .FALSE.
        td_obs%i_ts     = id_ts
        td_obs%r_date   = id_ts*td_cop%r_dt
        CALL debug(id_ts, "In obs_sample_and_save_2D, saving obs for time step :", tag=dTRACE)
        CALL debug(td_obsOut%filename, "  td_obsOut%filename = ", tag=dTRACE)
        CALL ncObsWrite(td_obsOut, td_obs, rec=il_idx)
      END IF
    END IF
  END SUBROUTINE obs_sample_and_save_2D

END MODULE ncconv_obs