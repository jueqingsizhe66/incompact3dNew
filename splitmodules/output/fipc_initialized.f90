  Subroutine FIPC_initialized( flag, error )

    ! Test whether FreeIPC has been initialized.
    ! On success return .True. in FLAG if is has, otherwise .FALSE., and
    ! set ERROR to FIPC_SUCCESS
    ! On error FLAG is undefined and ERROR is any value but FIPC_SUCCESS

    Logical, Intent( Out ) :: flag
    Integer, Intent( Out ) :: error

    flag = Associated( FIPC_ctxt_world%initialized )

    error = FIPC_success

  End Subroutine FIPC_initialized
