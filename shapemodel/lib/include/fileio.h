#ifndef	_fileio_
# define	_fileio_
//
//	fileio.h	- Input/Output from/to files
//
//	Oliver Grau
//

#include	<stdio.h>

/*! 
  \class meta_file fileio.h
  \brief base class for fileio
  */

class	meta_file {
	public:
		int  error() const;
		void set_errno(int e);
		void clear_errno();
	protected:
		meta_file();
		int	m_errno;
};

/*!
  \class meta_read_file fileio.h
  */

class	meta_read_file: public meta_file	{
	public:
		virtual	unsigned	read( char &)=0;
		virtual	unsigned	read( int &)=0;
		virtual	unsigned	read( double &)=0;
		virtual	unsigned	read( char *, int max )=0;
		virtual	unsigned	read( void *, unsigned, unsigned )=0;
		virtual	long	fpos()=0;
		virtual	void	seek(long)=0;
	protected:
		meta_read_file () {};
		virtual ~meta_read_file() {};
};

/*!
  \class meta_write_file fileio.h
  */

class	meta_write_file: public meta_file	{
	public:
		virtual void	write( char )=0;
		virtual void	write( int )=0;
		virtual void	write( double )=0;
		virtual void	write( const char * )=0;
		virtual void	write( const void *, unsigned, unsigned )=0;
		virtual long	fpos()=0;
		virtual void	seek(long)=0;
		virtual void	flush()=0;
	protected:
		meta_write_file () {};
		virtual ~meta_write_file() {};
};


/*!
  \class read_file fileio.h
  \brief Class for file read access
  */

class	read_file: public meta_read_file	{
         public:
                //! create file from existing FILE
		read_file (FILE *fp);

		//! open file with the name fn for reading
		read_file (const char *fn);

		virtual ~read_file ();

		/*! read a single character. return 1 on success
		  and 0 otherwise */
		virtual	unsigned	read( char &);

		/*! read an integer and leading skipps whitespaces.
		return 1 on success and 0 otherwise */
		virtual	unsigned	read( int &);

		/*! read an double and leading skipps whitespaces.
		return 1 on success and 0 otherwise */
		virtual	unsigned	read( double &);

		/*! read a character string. The leading '\n' will be
		  stripped. returns the string length */
		virtual	unsigned	read( char *, int max );

		/*! read binary data to buffer b, of size siz*n.
		  returns n on success */
		virtual	unsigned	read( void *, unsigned, unsigned );

		//! returns actual position of the filepointer
		virtual	long	fpos();

		//! sets filepointer to pos
		virtual	void	seek(long);
	protected:
		FILE	*fp;
		int	ext;
};



/*!
  \class write_file fileio.h
  \brief Class for file write access
  */

class	write_file: public meta_write_file	{
	public:
                //! create file from existing FILE
		write_file (FILE *fp);

		//! open file with the name fn for writing
		write_file (const char *fn, int append=0 );

		virtual ~write_file ();

		//! write a single character.
		virtual	void	write( char );

		//! write an integer
		virtual	void	write( int );

		//! 	write a double
		virtual	void	write( double );
		
		//! write a character string
		virtual	void	write( const char * );

		//! write binary data to buffer b, of size siz*n
		virtual	void	write( const void *, unsigned, unsigned );

		//! returns actual position of the filepointer
		virtual	long	fpos();

		//! sets filepointer to pos
		virtual	void	seek(long);

		//! flush file buffer
		virtual	void	flush();
	protected:
		FILE	*fp;
		int	ext;
};


#endif
