!!Special module for boundary conditions (netcdf part)
!!@author Innocent Souopgui
MODULE ncbc_tools
  USE general_constant
  USE general_tools
  USE bc_tools
  USE nctools
  USE debug_tools
IMPLICIT NONE

  !> \brief User defined type for attributes of boundary condition
  !<
  TYPE bc_structure_components
    CHARACTER(LEN=ip_cnl) ::&
      aa_title   = "title",&
      aa_history = "history",&
      l_by_axis  = "group_by_axis",&
      i_ndim = "ndim",&!nx     = "nx",&
      nxm2   = "nxm2",&
      ny     = "ny",&
      nym2   = "nym2",&
      nz     = "nz",&
      r_dt   = "delta_t"

    CHARACTER(LEN=ip_cnl) ::&
      r_date   = "t",&
      i_ts   = "timeStep"
    CHARACTER(LEN=ip_cnl) ::&
      r_dateu  = "s",&
      i_tsu  = "N/A"
    CHARACTER(LEN=ip_lnl) ::&
      r_dateln = "time",&
      i_tsln = "time step ordinal number"

    CHARACTER(LEN=ip_cnl) ::&
      ra_data   = "data",&
      ra_xmin_bc   = "xmin_data",&
      ra_xmax_bc   = "xmax_data",&
      ra_ymin_bc   = "ymin_data",&
      ra_ymax_bc   = "ymax_data",&
      ra_zmin_bc   = "zmin_data",&
      ra_zmax_bc   = "zmax_data"
    CHARACTER(LEN=ip_cnl) ::&
      ra_datau  = "N/A",&
      ra_xmin_bcu  = "N/A",&
      ra_xmax_bcu  = "N/A",&
      ra_ymin_bcu  = "N/A",&
      ra_ymax_bcu  = "N/A",&
      ra_zmin_bcu  = "N/A",&
      ra_zmax_bcu  = "N/A"
    CHARACTER(LEN=ip_lnl) ::&
      ra_dataln = "Boundary condition data, when not grouped by axis",&
      ra_xmin_bcln = "Boundary condition data for xmin face",&
      ra_xmax_bcln = "Boundary condition data for xmax face",&
      ra_ymin_bcln = "Boundary condition data for ymin face",&
      ra_ymax_bcln = "Boundary condition data for ymax face",&
      ra_zmin_bcln = "Boundary condition data for zmin face",&
      ra_zmax_bcln = "Boundary condition data for zmax face"
  END TYPE bc_structure_components

  !>
  !! Do not the fields nx, ny and nz for other purposes as their value may not be accurate after reading from file,
  !! particulary if the dimensionality of the problem is less than 3
  !<
  TYPE bcOutput
    CHARACTER(LEN=ip_snl) :: fileName
    CHARACTER(LEN=ip_snl) :: aa_title
    CHARACTER(LEN=ip_fnl) :: aa_history

    INTEGER :: ncid    = -1
    !basically, each of the following number is equal to the number of point in the corresponding axis of the computational domain plus 2(corners). Remember that the boundary data are set in the ghost points.
    INTEGER, PRIVATE :: nx = -1!number of data on the x boundary (include the corners)
    INTEGER, PRIVATE :: ny = -1!number of data on the y boundary (include the corners)
    INTEGER, PRIVATE :: nz = -1!number of data on the z boundary (include the corners)
    INTEGER :: i_ndim  = -1
    REAL(dp):: r_dt    = -1.0d0
    LOGICAL :: l_by_axis
    INTEGER :: nb_record = -1 !number of record in the file

    INTEGER :: ra_dataid    = -1
    INTEGER :: ra_xmin_bcid = -1
    INTEGER :: ra_xmax_bcid = -1
    INTEGER :: ra_ymin_bcid = -1
    INTEGER :: ra_ymax_bcid = -1
    INTEGER :: ra_zmin_bcid = -1
    INTEGER :: ra_zmax_bcid = -1
    INTEGER :: i_tsid       = -1
    INTEGER :: r_dateid     = -1

    LOGICAL :: isOpened  = .FALSE.
    INTEGER :: i_nextRec = -1
  END TYPE bcOutput

  TYPE(bc_structure_components), PRIVATE, SAVE :: tm_bcAtt
  PRIVATE check_conformance
CONTAINS

  !> \brief print data structure for boundary conditions ntcdf file
  !! @param [in] td_ncfile boundary condition netcdf file data structure
  !<
  SUBROUTINE print_bcOutput(td_ncfile)
    type(bcOutput), INTENT(IN) :: td_ncfile

    CALL debug('', 'printing bcOutput---------------------------------', tag=dALLWAYS)
    CALL debug(td_ncfile%isOpened  , '  isOpened   = ', tag=dALLWAYS)
    CALL debug(td_ncfile%fileName  , '  fileName   = ', tag=dALLWAYS)
    CALL debug(td_ncfile%aa_title  , '  aa_title   = ', tag=dALLWAYS)
    CALL debug(td_ncfile%aa_history, '  aa_history = ', tag=dALLWAYS)
    CALL debug(td_ncfile%nx        , '  nx         = ', tag=dALLWAYS)
    CALL debug(td_ncfile%ny        , '  ny         = ', tag=dALLWAYS)
    CALL debug(td_ncfile%nz        , '  nz         = ', tag=dALLWAYS)
    CALL debug(td_ncfile%i_ndim    , '  i_ndim     = ', tag=dALLWAYS)
    CALL debug(td_ncfile%r_dt      , '  r_dt       = ', tag=dALLWAYS)
    CALL debug(td_ncfile%l_by_axis , '  l_by_axis  = ', tag=dALLWAYS)
    CALL debug(td_ncfile%nb_record , '  nb_record  = ', tag=dALLWAYS)
    CALL debug(td_ncfile%i_nextRec , '  i_nextRec  = ', tag=dALLWAYS)
    CALL debug('', '.......................................................', tag=dALLWAYS)
  END SUBROUTINE print_bcOutput

  !> \brief initializes data structure that is used to saved boundary conditions
  !! @param [in,out] td_ncfile boundary condition netcdf file data structure
  !! @param [in] td_bcp parameters to use initialized
  !! @param [in] status (optional) status of the file (INPUT_FILE or OUTPUT_FILE), default is OUTPUT_FILE
  !<
  SUBROUTINE initBCOutput(td_ncfile, td_bcp, status)
    type(bcOutput), INTENT(INOUT) :: td_ncfile
    TYPE(bc_param), INTENT(IN)    :: td_bcp
    INTEGER, OPTIONAL, INTENT(IN) :: status
    !local variables
    !CHARACTER(len=ip_fnl) :: ala_fName
    INTEGER :: il_status
    IF( .NOT.td_bcp%l_by_axis ) CALL stop_program('In init_bcs, only grouping by axis is supported for now')
    IF(PRESENT(status))THEN
      il_status = status
    ELSE
      il_status = OUTPUT_FILE
    END IF

    td_ncfile%fileName = make_fileName(BC_DATA, il_status, prefix=td_bcp%aa_prefix)
    !td_ncfile% = td_bcp%
    td_ncfile%aa_history = td_bcp%aa_hist
    td_ncfile%aa_title   = td_bcp%aa_title
    td_ncfile%l_by_axis  = td_bcp%l_by_axis
    td_ncfile%i_ndim     = td_bcp%i_ndim
    td_ncfile%r_dt       = td_bcp%r_dt
    IF( td_bcp%i_ndim>=1 ) THEN
      td_ncfile%nx = td_bcp%ia_nxyz(1)
    ELSE
      td_ncfile%nx = 1
    END IF
    IF( td_bcp%i_ndim>=2 ) THEN
      td_ncfile%ny = td_bcp%ia_nxyz(2)
    ELSE
      td_ncfile%ny = 1
    END IF
    IF( td_bcp%i_ndim>=3 ) THEN
      td_ncfile%nz = td_bcp%ia_nxyz(3)
    ELSE
      td_ncfile%nz = 1
    END IF
    !
    td_ncfile%ncid = -1
    td_ncfile%i_nextRec = -1
    td_ncfile%isOpened = .FALSE.
  END SUBROUTINE initBCOutput

  SUBROUTINE ncBCClose(td_ncfile)
    type(bcOutput), INTENT(INOUT) :: td_ncfile

    CALL chkerr( nf90_close( td_ncfile%ncid), fname=td_ncfile%filename )
    td_ncfile%ncid = -1
    td_ncfile%i_nextRec = -1
    td_ncfile%isOpened = .FALSE.
  END SUBROUTINE ncBCClose

  SUBROUTINE ncBCOpen(td_ncfile)
    TYPE(bcOutput), INTENT(INOUT) :: td_ncfile
    !local variables
    INTEGER :: il_nDimensions, il_nVariables, il_nAttributes, il_unlimitedDimId, il_formatNum,&
      il_nxm2id, il_nyid, il_nzid!, il_ndimid
    INTEGER :: il_nxm2
    character(len = nf90_max_name) :: lc_name

    CALL debug( TRIM(td_ncfile%fileName), 'In ncBCOpen, opening -- ', tag=dNETCDF )

    !CALL debug(100, 'In ncBCOpen * ', tag=dALLWAYS)
    CALL chkerr(nf90_open(TRIM(td_ncfile%filename),NF90_NOCLOBBER,td_ncfile%ncid), fname=td_ncfile%filename)
    !CALL debug(200, 'In ncBCOpen * ', tag=dALLWAYS)
    CALL read_bcAtt(td_ncfile)

    !CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_data   , td_ncfile%ra_dataid    ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_xmin_bc, td_ncfile%ra_xmin_bcid ) )
    !CALL debug(300, 'In ncBCOpen * ', tag=dALLWAYS)
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_xmax_bc, td_ncfile%ra_xmax_bcid ) )
    !CALL debug(400, 'In ncBCOpen * ', tag=dALLWAYS)
    IF(td_ncfile%i_ndim >= 2)THEN
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_ymin_bc, td_ncfile%ra_ymin_bcid ) )
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_ymax_bc, td_ncfile%ra_ymax_bcid ) )
    ELSE
      td_ncfile%ra_ymin_bcid = -1
      td_ncfile%ra_ymax_bcid = -1
    END IF
    !CALL debug(500, 'In ncBCOpen * ', tag=dALLWAYS)
    IF(td_ncfile%i_ndim == 3)THEN
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_zmin_bc, td_ncfile%ra_zmin_bcid ) )
      CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%ra_zmax_bc, td_ncfile%ra_zmax_bcid ) )
    ELSE
      td_ncfile%ra_zmin_bcid = -1
      td_ncfile%ra_zmax_bcid = -1
    END IF
    !CALL debug(600, 'In ncBCOpen * ', tag=dALLWAYS)
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%r_date, td_ncfile%r_dateid ) )
    CALL chkerr( nf90_inq_varid( td_ncfile%ncid, tm_bcAtt%i_ts  , td_ncfile%i_tsid   ) )

    CALL chkerr(nf90_inquire(td_ncfile%ncid, il_nDimensions, il_nVariables, il_nAttributes, &
        il_unlimitedDimId, il_formatNum))

    CALL chkerr(nf90_inquire_dimension(td_ncfile%ncid, il_unlimitedDimId, name = lc_name, len = td_ncfile%nb_record))
    !CALL debug( td_ncfile%nb_record, 'In ncBCOpen, td_ncfile%nb_record = ', tag=dNETCDF )
    CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_bcAtt%ny  , il_nyid   ) )
    CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_nyid  , lc_name, td_ncfile%ny ) )
    CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_bcAtt%nz  , il_nzid   ) )
    CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_nzid  , lc_name, td_ncfile%nz ) )
    !CALL debug(700, 'In ncBCOpen * ', tag=dALLWAYS)
    IF(td_ncfile%i_ndim >= 2)THEN
      CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_bcAtt%nxm2, il_nxm2id ) )
      CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_nxm2id, lc_name, il_nxm2      ) )
      td_ncfile%nx = il_nxm2 +2
    ELSE
      td_ncfile%nx = 1
    END IF
    !CALL debug(800, 'In ncBCOpen * ', tag=dALLWAYS)
    !CALL chkerr( nf90_inq_dimid( td_ncfile%ncid, tm_bcAtt%i_ndim, il_ndimid ) )
    !CALL chkerr( nf90_inquire_dimension( td_ncfile%ncid, il_ndimid, lc_name, td_ncfile%i_ndim ) )

    td_ncfile%isOpened  = .TRUE.
    td_ncfile%i_nextRec = 1
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncBCOpen

  SUBROUTINE read_bcAtt( td_ncfile )
    !read des attributs (param�tres) de la trajectoire
    type(bcOutput), INTENT(INOUT) :: td_ncfile
    !local variables
    INTEGER :: il_tmp
    CALL debug("In read_bcAtt, reading attributes", tag=dNETCDF )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%r_dt      , td_ncfile%r_dt       ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%aa_title  , td_ncfile%aa_title   ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%aa_history, td_ncfile%aa_history ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%i_ndim    , td_ncfile%i_ndim     ) )
    CALL chkerr( nf90_get_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%l_by_axis , il_tmp               ) )
    td_ncfile%l_by_axis = int2l(il_tmp)
    CALL debug("End of attributes reading", tag=dNETCDF )
    !End of saving
  END SUBROUTINE read_bcAtt

  SUBROUTINE ncBCCreate(td_ncfile)
    type(bcOutput), INTENT(INOUT) :: td_ncfile
    INTEGER, DIMENSION(3) :: ila_dims3D
    !INTEGER, DIMENSION(2) :: ila_dims2D
    INTEGER :: il_nxm2, il_ny, il_nym2, il_nz, il_date!, il_ndim

    CALL debug( TRIM(td_ncfile%fileName), 'In ncBCCreate, creating -- ', tag=dNETCDF )
    CALL chkerr(nf90_create (TRIM(td_ncfile%fileName), NF90_CLOBBER, td_ncfile%ncid   ), fname=td_ncfile%filename)
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%ny    , td_ncfile%ny  , il_ny    ))
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%nz    , td_ncfile%nz  , il_nz    ))
    IF(td_ncfile%i_ndim >= 2)&
      CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%nxm2  , td_ncfile%nx-2, il_nxm2  ))
    IF(td_ncfile%i_ndim == 3)&
      CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%nym2  , td_ncfile%ny-2, il_nym2  ))
    !CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%i_ndim, td_ncfile%i_ndim, il_ndim))
    CALL chkerr(nf90_def_dim(td_ncfile%ncid, tm_bcAtt%r_date, NF90_UNLIMITED  , il_date))

    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%r_date, NF90_DOUBLE, il_date, td_ncfile%r_dateid))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%i_ts  , NF90_INT   , il_date, td_ncfile%i_tsid  ))
    ila_dims3D = (/il_ny, il_nz, il_date/)
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_xmin_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_xmin_bcid ))
    CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_xmax_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_xmax_bcid ))
    IF(td_ncfile%i_ndim >= 2)THEN
      ila_dims3D = (/il_nxm2, il_nz, il_date/)
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_ymin_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_ymin_bcid ))
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_ymax_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_ymax_bcid ))
    END IF
    IF(td_ncfile%i_ndim == 3)THEN
      ila_dims3D = (/il_nxm2, il_nym2, il_date/)
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_zmin_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_zmin_bcid ))
      CALL chkerr(nf90_def_var(td_ncfile%ncid, tm_bcAtt%ra_zmax_bc, NF90_DOUBLE, ila_dims3D, td_ncfile%ra_zmax_bcid ))
    END IF
    CALL save_bcAtt(td_ncfile)
    CALL chkerr(nf90_enddef(td_ncfile%ncid))
    td_ncfile%isOpened  = .TRUE.
    td_ncfile%i_nextRec = 1
    CALL debug( '... ncBCCreate -> done', tag=dNETCDF )
  END SUBROUTINE ncBCCreate

  SUBROUTINE save_bcAtt(td_ncfile)
    !sauvegarde des attributs (parametres) de la trajectoire
    type(bcOutput), INTENT(INOUT) :: td_ncfile
    !local variables
    INTEGER :: il_tmp
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%r_dt      , td_ncfile%r_dt      ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%aa_title  , td_ncfile%aa_title  ) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%aa_history, td_ncfile%aa_history) )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%i_ndim    , td_ncfile%i_ndim    ) )
    il_tmp = logical2int ( td_ncfile%l_by_axis )
    CALL chkerr( nf90_put_att( td_ncfile%ncid, NF90_GLOBAL, tm_bcAtt%l_by_axis, il_tmp ) )
    CALL debug( "End of attributes recording", tag=dNETCDF )
    !End of saving
  END SUBROUTINE save_bcAtt

  !> \brief write boundary conditions data to netcdf file
  !! @param [in,out] td_ncfile file to be written
  !! @param [in,out] td_bc data structure that contains the data to write
  !! @param [in] rec (optional) gives the ordinal number of the record to write
  !<
  SUBROUTINE ncBCWrite(td_ncfile, td_bc, rec)
    TYPE(bcOutput), INTENT(INOUT)  :: td_ncfile
    TYPE(bc_structure), INTENT(IN) :: td_bc
    INTEGER, INTENT(IN), OPTIONAL  :: rec
    !local variables
    INTEGER, DIMENSION(3) :: ila_count3D, ila_start3D
    INTEGER, DIMENSION(2) :: ila_start2D!, ila_count2D
    INTEGER, DIMENSION(1) :: ila_start1D
    INTEGER               :: il_rec!, il_icoord, il_rcoord
    REAL(cp) :: rl_date

    CALL check_conformance(td_bc, td_ncfile)
    IF( PRESENT(rec) )THEN
      il_rec = rec
    ELSE
      il_rec = td_ncfile%i_nextRec
    END IF
    rl_date = td_bc%i_ts*td_ncfile%r_dt
    CALL debug(il_rec, 'In nBCWrite, writing '//TRIM(td_ncfile%filename)//' - record', tag=dNETCDF )
    ila_start1D = (/il_rec/)
    ila_start2D = (/ 1, il_rec /)
    ila_start3D = (/ 1, 1, il_rec /)
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%r_dateid, rl_date   , start = ila_start1D))
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%i_tsid  , td_bc%i_ts, start = ila_start1D))
    ila_count3D = (/td_ncfile%ny, td_ncfile%nz, 1/)
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_xmin_bcid, td_bc%ra_xmin_bc, start = ila_start3D, COUNT = ila_count3D))
    CALL chkerr(nf90_put_var(td_ncfile%ncid, td_ncfile%ra_xmax_bcid, td_bc%ra_xmax_bc, start = ila_start3D, COUNT = ila_count3D))
    IF(td_ncfile%i_ndim >= 2)THEN
      ila_count3D = (/td_ncfile%nx-2, td_ncfile%nz, 1/)
      CALL chkerr(&
        nf90_put_var(td_ncfile%ncid, td_ncfile%ra_ymin_bcid, td_bc%ra_ymin_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
      CALL chkerr(&
        nf90_put_var(td_ncfile%ncid, td_ncfile%ra_ymax_bcid, td_bc%ra_ymax_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
    END IF
    IF(td_ncfile%i_ndim == 3)THEN
      ila_count3D = (/td_ncfile%nx-2, td_ncfile%ny-2, 1/)
      CALL chkerr(&
        nf90_put_var(td_ncfile%ncid, td_ncfile%ra_zmin_bcid, td_bc%ra_zmin_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
      CALL chkerr(&
        nf90_put_var(td_ncfile%ncid, td_ncfile%ra_zmax_bcid, td_bc%ra_zmax_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
    END IF

    td_ncfile%i_nextRec = il_rec + 1
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncBCWrite

  !> \brief read boundary conditions from netcdf file
  !! @param [in,out] td_ncfile file to be read from
  !! @param [in,out] td_bc data structure to contain the read data
  !! @param [in] rec (optional) gives the ordinal number of the record to read
  !<
  SUBROUTINE ncBCRead(td_ncfile, td_bc, rec)
    type(bcOutput), INTENT(INOUT) :: td_ncfile
    type(bc_structure), INTENT(INOUT) :: td_bc
    INTEGER, INTENT(IN), OPTIONAL  :: rec

    INTEGER, DIMENSION(3) :: ila_count3D, ila_start3D
    INTEGER, DIMENSION(2) :: ila_start2D!, ila_count2D
    INTEGER, DIMENSION(1) :: ila_start1D, ila_count1D
    !REAL(KIND=sp), DIMENSION(1)  :: rl_date
    INTEGER :: il_rec!, il_rcoord, il_icoord
    !LOGICAL :: ll_icoord, ll_rcoord

    CALL check_conformance(td_bc, td_ncfile)
    !CALL print_bcs(td_bc)
    IF( PRESENT(rec) )THEN
      il_rec = rec
    ELSE
      il_rec = td_ncfile%i_nextRec
    END IF
    CALL debug(il_rec, 'In ncBCRead, Reading '//TRIM(td_ncfile%filename)//' - record', tag=dNETCDF )
    ila_start1D = (/il_rec/)
    ila_count1D = (/1/)
    ila_start2D = (/ 1, il_rec /)
    ila_start3D = (/ 1, 1, il_rec /)

    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%r_dateid, td_bc%r_date, start = ila_start1D))
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%i_tsid  , td_bc%i_ts  , start = ila_start1D))
    ila_count3D = (/td_ncfile%ny, td_ncfile%nz, 1/)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_xmin_bcid, td_bc%ra_xmin_bc, start = ila_start3D, COUNT = ila_count3D))
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%ra_xmax_bcid, td_bc%ra_xmax_bc, start = ila_start3D, COUNT = ila_count3D))
    IF(td_ncfile%i_ndim >= 2)THEN
      ila_count3D = (/td_ncfile%nx-2, td_ncfile%nz, 1/)
      CALL chkerr(&
        nf90_get_var(td_ncfile%ncid, td_ncfile%ra_ymin_bcid, td_bc%ra_ymin_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
      CALL chkerr(&
        nf90_get_var(td_ncfile%ncid, td_ncfile%ra_ymax_bcid, td_bc%ra_ymax_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
    END IF
    IF(td_ncfile%i_ndim == 3)THEN
      ila_count3D = (/td_ncfile%nx-2, td_ncfile%ny-2, 1/)
      CALL chkerr(&
        nf90_get_var(td_ncfile%ncid, td_ncfile%ra_zmin_bcid, td_bc%ra_zmin_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
      CALL chkerr(&
        nf90_get_var(td_ncfile%ncid, td_ncfile%ra_zmax_bcid, td_bc%ra_zmax_bc, start = ila_start3D, COUNT = ila_count3D)&
      )
    END IF

    td_ncfile%i_nextRec = il_rec + 1
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE ncBCRead

  SUBROUTINE readBC_date(td_ncfile, rda_date)
    type(bcOutput),INTENT(IN)              :: td_ncfile
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT):: rda_date

    CALL debug(TRIM(td_ncfile%filename), 'In readBC_date : reading date from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%r_dateid, rda_date))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readBC_date

  SUBROUTINE readBC_timeStep(td_ncfile, id_recNum, id_ts)
    type(bcOutput),INTENT(IN):: td_ncfile
    INTEGER, INTENT(IN)       :: id_recNum
    INTEGER, INTENT(OUT)      :: id_ts

    CALL debug(TRIM(td_ncfile%filename), 'In readBC_timeStep : reading timeStep from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%i_tsid, id_ts, start=(/id_recNum/) ))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readBC_timeStep

  SUBROUTINE readBC_AllTimeStep(td_ncfile, ida_ts)
    type(bcOutput),INTENT(IN)        :: td_ncfile
    INTEGER, DIMENSION(:), INTENT(OUT):: ida_ts

    CALL debug(TRIM(td_ncfile%filename), 'In readBC_AllTimeStep : reading timeStep from ', tag=dNETCDF)
    CALL chkerr(nf90_get_var(td_ncfile%ncid, td_ncfile%i_tsid, ida_ts))
    CALL debug('... done', tag=dNETCDF)
  END SUBROUTINE readBC_AllTimeStep

  !> \brief this routine checks the conformance of the bd data structure and the netcdf data structure
  !! @param [in] td_bc BC data structure
  !! @param [in] td_ncfile netcdf data structure for BC
  !! The program is stopped if there is conformance
  !<
  SUBROUTINE check_conformance(td_bc, td_ncfile)
    TYPE(bc_structure), INTENT(IN) :: td_bc
    type(bcOutput), INTENT(IN)     :: td_ncfile
    !local variables
    INTEGER, DIMENSION(2) :: ila_shape

!       ila_shape(1) = SIZE(td_ncfile%ra_ymin_bcid,1)
!       ila_shape(2) = SIZE(td_ncfile%ra_ymin_bcid,2)
    CALL debug('Entering check_conformance', tag=dNETCDF)
    !check only the xmin since xmax is allocate at the same time
    ila_shape = SHAPE(td_bc%ra_xmin_bc)
!     CALL debug ( ila_shape, ' ila_shape = ', tag=dALLWAYS )
!     CALL debug ( (/td_ncfile%ny, td_ncfile%nz/), ' (/td_ncfile%ny, td_ncfile%nz/) = ', tag=dALLWAYS )
    !CALL dpause( 'Pause in check_conformance' )
    IF( ANY( ila_shape/=(/td_ncfile%ny, td_ncfile%nz/) ) )THEN
      CALL debug(ila_shape, 'ila_shape = ', tag=dNETCDF)
      CALL stop_program("In check_conformance: non conformance of the xmin_bc with the netcdf file")
    END IF
    IF(td_ncfile%i_ndim >= 2)THEN
      ila_shape = SHAPE(td_bc%ra_ymin_bc)
      IF( ANY( ila_shape/=(/td_ncfile%nx-2, td_ncfile%nz/) ) )THEN
        CALL debug(ila_shape, 'ila_shape = ', tag=dNETCDF)
        CALL stop_program("In check_conformance: non conformance of the ymin_bc with the netcdf file")
      END IF
    END IF
    IF(td_ncfile%i_ndim == 3)THEN
      ila_shape = SHAPE(td_bc%ra_zmin_bc)
      IF( ANY( ila_shape/=(/td_ncfile%nx-2, td_ncfile%ny-2/) ) )THEN
        CALL debug(ila_shape, 'ila_shape = ', tag=dNETCDF)
        CALL stop_program("In check_conformance: non conformance of the zmin_bc with the netcdf file")
      END IF
    END IF
    CALL debug('Exiting check_conformance', tag=dNETCDF)
  END SUBROUTINE check_conformance

END MODULE ncbc_tools