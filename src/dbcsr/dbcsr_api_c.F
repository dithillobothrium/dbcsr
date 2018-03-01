module dbcsr_api_c
   use, intrinsic :: ISO_C_BINDING, only : c_loc, c_ptr, c_double, c_f_pointer, c_int, c_char, c_null_ptr
   USE dbcsr_api
   USE machine,                         ONLY: default_output_unit
   
   implicit none
   private
   
   contains

      function c_f_string(c_str) result(f_str)
        use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char
        type(c_ptr), intent(in) :: c_str
        character(:,kind=c_char), pointer :: f_str
        character(kind=c_char), pointer :: arr(:)
        interface
          ! steal std c library function rather than writing our own.
          function strlen(s) bind(c, name='strlen')
            use, intrinsic :: iso_c_binding, only: c_ptr, c_size_t
            implicit none
            !----
            type(c_ptr), intent(in), value :: s
            integer(c_size_t) :: strlen
          end function strlen
        end interface
        !****
        call c_f_pointer(c_str, arr, [strlen(c_str)])
        call get_scalar_pointer(size(arr), arr, f_str)
      end function c_f_string
      subroutine get_scalar_pointer(scalar_len, scalar, ptr)
        use, intrinsic :: iso_c_binding, only: c_char
        integer, intent(in) :: scalar_len
        character(kind=c_char,len=scalar_len), intent(in), target :: scalar(1)
        character(:,kind=c_char), intent(out), pointer :: ptr
        !***
        ptr => scalar(1)
      end subroutine get_scalar_pointer

      subroutine c_dbcsr_init_lib() bind(C, name="c_dbcsr_init_lib")
        print *,"init"
        call dbcsr_init_lib()
      end subroutine

      
      subroutine c_dbcsr_finalise_lib(fcomm) bind(C, name="c_dbcsr_finalize_lib_aux")
        integer, intent(in), value :: fcomm
        integer :: rank
        integer:: io_stream
        io_stream = 0
        print *,"finalize"
        call mpi_comm_rank(fcomm,rank)
        print *,"rank ", rank
        if (rank .eq. 0) then
            io_stream = default_output_unit
        endif
        call dbcsr_finalize_lib(fcomm, io_stream)
      end subroutine

      
      subroutine c_dbcsr_distribution_new(c_dist, group, c_row_dist, row_dist_size, &
                                    c_col_dist, col_dist_size) &
                                    bind(C, name="c_dbcsr_distribution_new")
        type(c_ptr), intent(out)                :: c_dist
        integer, intent(in), value              :: group
        integer, intent(in), value              :: row_dist_size, col_dist_size
        integer, intent(in), target             :: c_row_dist(row_dist_size), c_col_dist(col_dist_size)
        type(dbcsr_distribution_type), pointer  :: dist
        integer, pointer                        :: row_dist(:), col_dist(:)
        
        row_dist => c_row_dist
        col_dist => c_col_dist
        allocate(dist)
        call dbcsr_distribution_new(dist, group=group, row_dist=row_dist, &
                                    col_dist=col_dist, reuse_arrays=.false.)
        c_dist = c_loc(dist)
      end subroutine


      subroutine c_dbcsr_distribution_release(c_dist) bind(C, name="c_dbcsr_distribution_release")
        type(c_ptr), intent(inout)      :: c_dist
        type(dbcsr_distribution_type),  pointer          :: dist

        call c_f_pointer(c_dist, dist) 

        call dbcsr_distribution_release(dist)

        c_dist = c_null_ptr
      end subroutine


      subroutine test(c_name) bind(C, name="test")
        use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer, c_char
        type(c_ptr), value, intent(in) :: c_name
        print *, c_f_string(c_name)                           
      end subroutine
      
      
      subroutine c_dbcsr_create_new_d(c_matrix, c_name, c_dist, c_matrix_type, &
                                      c_row_blk_sizes, c_row_blk_sizes_length, &
                                      c_col_blk_sizes, c_col_blk_sizes_length) &
                                      bind(C, name="c_dbcsr_create_new_d")
        type(c_ptr), intent(out)            :: c_matrix
        type(c_ptr), value, intent(in)      :: c_name
        type(c_ptr), value, intent(in)      :: c_dist  
        character, value, intent(in)        :: c_matrix_type
        integer, intent(in), value          :: c_row_blk_sizes_length, c_col_blk_sizes_length
        integer, intent(in), target         :: c_row_blk_sizes(c_row_blk_sizes_length)
        integer, intent(in), target         :: c_col_blk_sizes(c_col_blk_sizes_length) 
        
        type(dbcsr_type),  pointer              :: matrix
        type(dbcsr_distribution_type), pointer  :: dist
        integer, pointer                        :: row_blk_sizes(:), col_blk_sizes(:) 
        !character(len=*),pointer                :: fname

        character(:,kind=c_char), pointer :: fname
        
        ! associate pointers
        row_blk_sizes => c_row_blk_sizes
        col_blk_sizes => c_col_blk_sizes

        ! convert normal distribution pointer to fortran pointer
        call c_f_pointer(c_dist, dist)

        ! allocate heap memory for matrix structure
        allocate(matrix)

        ! convert normal name to fortran name
        fname = c_f_string(c_name)

        print *,fname
        call dbcsr_create(matrix = matrix, name = fname, dist = dist, &
                          matrix_type = c_matrix_type, row_blk_size = row_blk_sizes, &
                          col_blk_size = col_blk_sizes,  data_type=dbcsr_type_real_8)

        c_matrix = c_loc(matrix)
      end subroutine


      subroutine c_dbcsr_release(c_matrix) bind(C, name="c_dbcsr_release")
        type(c_ptr), intent(inout)      :: c_matrix
        type(dbcsr_type),  pointer          :: matrix

        call c_f_pointer(c_matrix, matrix) 

        call dbcsr_release(matrix)

        c_matrix = c_null_ptr
      end subroutine


end module
