//
//	o_meta.h	- Header for oGeM
//
//	Oliver Grau, Jan.1993
//	Oliver Grau, Dec.1995
//	Oliver Grau, May-1996
//

#ifndef	_o_meta_incl
#  define	_o_meta_incl

typedef	unsigned long	UniqueIdentifier;
typedef	double		o_Time;

/*!
  \class o_MetaClass o_MetaClass.h
  \brief Abstract class. Keeps objectname, class name and handles <A HREF="http://www.tnt.uni-hannover.de/project/eu/panorama/wpg1/data-structures/memory-management.html">memory-management</A>
  
  Super class of the most classes in the library. Contains some <A HREF="http://www.tnt.uni-hannover.de/project/eu/panorama/wpg1/data-structures/common.html">common</A> features and mechanisms for <A HREF="http://www.tnt.uni-hannover.de/project/eu/panorama/wpg1/data-structures/memory-management.html">memory-management</A>
  */

class o_MetaClass {
private:
  //!test
	int     refcnt;
	char	*name;
	UniqueIdentifier	id;
	o_Time	timestamp;
protected:
	o_MetaClass (char *nam);
 	virtual	~o_MetaClass ();
	int ReqDelete( o_MetaClass * optr );
	void SetIdentifier(UniqueIdentifier idin);
public:
        /*!    set object time stamp to 'intim'. Usually the time stamp
	       is set by stream i/o classes (\sa o_StreamIn ) */
	void	SetTimeStamp(o_Time intim);

	/*! 	returns time stamp. The time stamp is set during restoration
		with o_StreamIn or by operators dealing with temporal
		changes. The SetTimeStamp() method is not public */
	o_Time	GetTimeStamp();

        /*!	copy method. Each child class should implement an
		own method if nessessary. */
        void	Copy( const o_MetaClass &org );

	/*!  returns identification key of the object p */
	o_key	ident(void *);	

        /*!  return identifier number. The identifier is a unique handle is usually automatically set. */
  	UniqueIdentifier	GetIdentifier() const;

	/*!  returns name of class */
	virtual	const	char	*GetClassName() const;

	/*!  returns name of object */
        const	char	*GetName() const;

        /*!  sets optional name of object */
	void	SetName(  const char * nam);

        /*!  creates a link from a oGeM-Object to the point */
	void Connect( o_MetaClass * optr );

        /*!  deletes a link from a oGeM-Object to the point */
	void DisConnect( o_MetaClass * optr );

        /*!  returns number of links from other objects */
	int	RefCount() const { return refcnt; };

};

#endif
