#pragma once
#include "voxlap5.h"

extern long loadsxl (const char *, char **, char **, char **);
extern void savesxl (char *, vx5sprite[] , long[]);

extern char *parspr (vx5sprite *, char **);
extern void loadnul (dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);
extern long loaddta (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);
extern long loadpng (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);
extern void loadbsp (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);

extern void pngoutputpixel (long);
extern void pngoutopenfile (const char *, long, long);

extern long loadsky (const char *);
extern long loadvxl (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);
extern long savevxl (const char *, dpoint3d *, dpoint3d *, dpoint3d *, dpoint3d *);

extern void relpathinit (char *);
extern char *stripdir (char *filnam);

extern void getspr (vx5sprite *, const char *);
extern kfatype *getkfa (const char *);

#if defined(_WIN32) && defined(USE_MEMMAPPED_IO)
	enum e_open_mode {
	    if_exists_fail_if_not_exists_create,
	    if_exists_keep_if_dont_exists_fail,
	    if_exists_keep_if_dont_exists_create,
	    if_exists_truncate_if_not_exists_fail,
	    if_exists_truncate_if_not_exists_create,
	};
	typedef struct{
	    char* data_;
	    size_t size_;
	    size_t capacity_;
	    size_t vbit_size_;
	    size_t vbit_capacity_;

	#if defined(__unix__)
		int file_handle_;
	} mem_mapped_file_info;// Note closure here 
	#elif defined(_WIN32)
	    typedef void * HANDLE;
	    HANDLE file_handle_;
	    HANDLE file_mapping_handle_;
	} mem_mapped_file_info;// Note closure here 
	#else
	#error Only Posix or Windows systems can use memory-mapped files.
	#endif
	extern void * allocate_mem_mapped_file( mem_mapped_file_info, char *, const size_t);
#endif

