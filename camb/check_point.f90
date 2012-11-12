module check_point
   use precision
implicit none
private

     !! This is the number of times the code will check point for each cl 
     !! computation (so 2x this total check points if you are doing the 
     !! scalar and tensor cls).
   integer :: NUM_CHECK_POINTS = 10

     !! Use these to determine the number of significant figures to keep in 
     !! checkpoint files.  If there are excessive unexplained validator failures
     !! you might try the 26.16 format.  This will increase the size of the 
     !! checkpoint file.
   !character(len=50) :: realfmt = 'E15.5'
   !character(len=50) :: realfmt = 'E16.7'
   character(len=50) :: realfmt = 'E20.10'
   !character(len=50) :: realfmt = 'E26.16'

     !! Turn on a few comments.
   logical :: Verbose = .false.
   
     !! Whether to write checkpoint files in Ascii or Binary.
   logical, parameter :: WriteAscii = .false.

   character(len=18) :: chkpnt_file
   logical :: resume_camb

   public :: NUM_CHECK_POINTS
   public :: check_point_init, check_point_resume, check_point_record, &
             check_point_clean

contains

   !!===================================================================================!!

   subroutine check_point_init(cl_type)

         character(len=*), intent(in) :: cl_type

         select case(trim(cl_type))
            case('scalar')
               chkpnt_file = 'camb_scalarcls.chk'
            case('tensor')
               chkpnt_file = 'camb_tensorcls.chk'
            case('vector')
               chkpnt_file = 'camb_vectorcls.chk'
         end select

         inquire(file=chkpnt_file,exist=resume_camb)

         if (Verbose) write(*,'(2A)') ' >> Checkpoint File = ', trim(chkpnt_file)
         if (Verbose) write(*,'(A,L)') ' >> Resuming = ', resume_camb

   end subroutine check_point_init

   !!===================================================================================!!

   subroutine check_point_resume(num_ks_resume, perc_resume, delta_p_l_k)

         integer, intent(out) :: num_ks_resume, perc_resume
         real(dl), dimension(:,:,:), intent(out) :: delta_p_l_k
  
         integer :: nsrcs, nls, nks, ios, i, j, ltmp, bnds(1:2)
         character(len=1024) :: str, frmt
  

         num_ks_resume = 0
         perc_resume = 0
         delta_p_l_k = 0.0_dl
         if (.not. resume_camb) return

         nsrcs = size(delta_p_l_k,1)
         nls = size(delta_p_l_k,2)
         nks = size(delta_p_l_k,3)
         
         if (WriteAscii) then
            open(unit=1,file=chkpnt_file,form='formatted')
         else
            open(unit=1,file=chkpnt_file,form='unformatted')
         end if
         do 
            if (WriteAscii) then
               read(1,'(A)',iostat=ios) str
            else
               read(1,iostat=ios) ltmp
            end if

            if (ios /= 0) then
               num_ks_resume = 0
               perc_resume = 0
               if (Verbose) write(*,'(A)') ' >> Unable to recover checkpoint file.'
               return
            end if

            if (WriteAscii) then
               if (trim(adjustl(str)) == 'finished') exit
            else
               if (ltmp == -999) exit
            end if
            num_ks_resume = num_ks_resume + 1
         end do
         rewind(1)

         num_ks_resume = num_ks_resume / (2*nsrcs)
         write(*,'(A,I5,A,I5,A)') ' >> Recovering values from ', num_ks_resume, &
                                  ' of ', nks, ' integrations.'
         
         do i = 1, num_ks_resume
            do j = 1, nsrcs
               if (WriteAscii) then
                  read(1,'(2I4)') bnds(1:2)
                  write(str,*) bnds(2) - bnds(1) + 1
                  frmt = '('//trim(adjustl(str))//trim(realfmt)//')'
                  read(1,frmt) delta_p_l_k(j,bnds(1):bnds(2),i)
               else
                  read(1) bnds(1:2)
                  read(1) delta_p_l_k(j,bnds(1):bnds(2),i)
               end if
            end do
         end do
         close(1)

         if (Verbose) write(*,'(A,I7)') ' >> Num_ks_Resume = ', num_ks_resume
         perc_resume = get_resume_percent(num_ks_resume,nks)

   end subroutine check_point_resume
   
   !!===================================================================================!!
   
   subroutine check_point_record(num_ks_record,delta_p_l_k)

         integer, intent(in) :: num_ks_record
         real(dl), dimension(:,:,:), intent(in) :: delta_p_l_k

         integer :: nsrcs, nls, i, j, bnds(1:2)
         character(len=1024) :: str, frmt


         if (Verbose) write(*,'(A)') ' >> Checkpointing'

         nsrcs = size(delta_p_l_k,1)
         nls = size(delta_p_l_k,2)
         
         if (WriteAscii) then
            open(unit=1,file=chkpnt_file,form='formatted')
         else
            open(unit=1,file=chkpnt_file,form='unformatted')
         end if
         do i = 1, num_ks_record
            do j = 1, nsrcs
               bnds = get_num_entries(delta_p_l_k(j,:,i))
               if (WriteAscii) then
                  write(1,'(2I4)') bnds(:)
                  write(str,*) bnds(2) - bnds(1) + 1
                  frmt = '('//trim(adjustl(str))//trim(realfmt)//')'
                  write(1,frmt) delta_p_l_k(j,bnds(1):bnds(2),i)
               else
                  write(1) bnds(:)
                  write(1) delta_p_l_k(j,bnds(1):bnds(2),i)
               end if
            end do
         end do
         if (WriteAscii) then
            write(1,'(A)') 'finished'
         else
            write(1) -999
         end if
         close(1)

         if (Verbose) write(*,'(A)') ' >> Finished Checkpointing'

   end subroutine check_point_record

   !!===================================================================================!!

   subroutine check_point_clean

         character(len=18) :: fname
         logical :: lex

         fname = 'camb_scalarcls.chk'
         inquire(file=fname,exist=lex)
         if (lex) then
            open(unit=1,file=fname)
            close(1,status='delete')
         end if
         
         fname = 'camb_tensorcls.chk'
         inquire(file=fname,exist=lex)
         if (lex) then
            open(unit=1,file=fname)
            close(1,status='delete')
         end if
         
         fname = 'camb_vectorcls.chk'
         inquire(file=fname,exist=lex)
         if (lex) then
            open(unit=1,file=fname)
            close(1,status='delete')
         end if

   end subroutine check_point_clean

   !!===================================================================================!!
 
   function get_num_entries(vec) result(bnds)
     !! Returns the number n of such that if k>n vec(k) = 0
         real(dl), dimension(:), intent(in) :: vec
         integer :: bnds(1:2)

         integer :: i, dimen, n
 
         dimen = size(vec)

         n = 1
         do i = 1, dimen
            if (vec(i) /= 0.0_dl) exit
            n = n + 1
         end do
         bnds(1) = n

         n = dimen
         do i = dimen, 1, -1
            if (vec(i) /= 0.0_dl) exit
            n = n - 1
         end do
         bnds(2) = n

   end function get_num_entries
   
   !!===================================================================================!!
 
   function get_resume_percent(num_ks_resume,num_ks_total) result(iperc)

         integer, intent(in) :: num_ks_resume, num_ks_total
         integer :: iperc

         integer :: itmp


         itmp = (num_ks_resume*100)/num_ks_total
         select case (itmp)
            case(:5)
               iperc = 0
            case(6:15)
               iperc = 10
            case(16:25)
               iperc = 20
            case(26:35)
               iperc = 30
            case(36:45)
               iperc = 40
            case(46:55)
               iperc = 50
            case(56:65)
               iperc = 60
            case(66:75)
               iperc = 70
            case(76:85)
               iperc = 80
            case(86:95)
               iperc = 90
            case(96:)
               iperc = 100
         end select 

   end function get_resume_percent

   !!===================================================================================!!

end module check_point
