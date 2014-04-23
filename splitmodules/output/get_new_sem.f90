  Subroutine get_new_sem( semid )

    ! Get a semaphore in exclusive mode
    !
    ! On success the semaphore id is returned in SEMID
    ! On error SEMID is negative. Possible values are
    ! the negative of the values that errno can be set to
    ! by shmget ( see
    !
    ! http://www.opengroup.org/onlinepubs/009695399/functions/shmget.html
    !
    ! ) if the shmget failed, or - Huge( 1 ) if despite trying
    ! multiple times we failed to create a seg

    Integer( c_int ), Intent(   Out ) :: semid

    Integer, Parameter :: max_tries = 200

    Integer            :: tries

    tries = 1

    Do
       ! Try to create a sem in exlusive mode. Note that fipc_seg_create
       ! tries to use a different value of semid each time it is called.
       semid = fipc_get_sem( SEG_CREATE, SEG_EXCLUDE, SEG_UREAD + SEG_UWRITE )
       If( semid /= EEXIST ) Then
          ! Either fipc_get_sem was succesful, or it failed for a reason other
          ! than the sem already exisiting
          Exit
       End If
       ! The sem we tried to create already exists. Try again ...
       tries = tries + 1
       If( tries > max_tries ) Then
          ! ... unless we've got bored. Make sure we can't go on for ever;
          ! As in the old geek joke we're sort of putting a known elephant in Cairo.
          ! http://en.wikipedia.org/wiki/Elephant_in_Cairo
          semid = - Huge( semid )
          Exit
       End If
    End Do

  End Subroutine get_new_sem
