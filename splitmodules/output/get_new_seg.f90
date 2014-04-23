  Subroutine get_new_seg( n, shmid )

    ! Get a segment of size N bytes in exclusive mode
    !
    ! On success the segment id is returned in SHMID
    ! On error SHMID is negative. Possible values are
    ! the negative of the values that errno can be set to
    ! by shmget ( see
    !
    ! http://www.opengroup.org/onlinepubs/009695399/functions/shmget.html
    !
    ! ) if the shmget failed, or - Huge( 1 ) if despite trying
    ! multiple times we failed to create a seg

    Integer( c_long ), Intent( In    ) :: n
    Integer( c_int  ), Intent(   Out ) :: shmid

    Integer, Parameter :: max_tries = 200

    Integer            :: tries

    tries = 1

    Do
       ! Try to create a seg in exlusive mode. Note that fipc_seg_create
       ! tries to use a different value of shmid each time it is called.
       shmid = fipc_get_seg( n, SEG_CREATE, SEG_EXCLUDE, SEG_UREAD + SEG_UWRITE )
       If( shmid /= EEXIST ) Then
          ! Either fipc_get_seg was succesful, or it failed for a reason other
          ! than the seg already exisiting
          Exit
       End If
       ! The seg we tried to create already exists. Try again ...
       tries = tries + 1
       If( tries > max_tries ) Then
          ! ... unless we've got bored. Make sure we can't go on for ever;
          ! As in the old geek joke we're sort of putting a known elephant in Cairo.
          ! http://en.wikipedia.org/wiki/Elephant_in_Cairo
          shmid = - Huge( tries )
          Exit
       End If
    End Do

  End Subroutine get_new_seg
