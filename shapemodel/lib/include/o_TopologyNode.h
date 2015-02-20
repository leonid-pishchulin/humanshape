#ifndef	_o_TopologyNode
#  define	_o_TopologyNode
//
//	o_TopologyNode.h	- Header for oGeM TopologyNode class
//
//	Oliver Grau, May 1994
//

class   o_MetaClass;

/*!
  \class o_TopologyNode o_TopologyNode.h
  \brief abstract class

  o_TopologyNode is an abstract class - do not instantiate
*/

class	o_TopologyNode {
	public:
		virtual	o_MetaClass	&Get_ParentNode();
		virtual	o_MetaClass	&Get_TopNode();

		/*! return state of the topology update flag. The flag is used
		for keeping objects or objects hierarchies up to date.
		\sa o_StreamOut . */
		int GetTopologyUpdateFlag() const ;

		/*! set topology update flag. The topology update flag must be
		set if any topological changes are made. It is typically 
		set automaticly in list access functions (like o_World::AddBody() ) */
		void SetTopologyUpdateFlag();

		void	Set_ParentNode( o_MetaClass &parentnode );
		void	Set_TopNode( o_MetaClass &topnode );
	protected:

		o_MetaClass	*parent;
		o_MetaClass	*top;

		// methods
		o_TopologyNode( );
		virtual ~o_TopologyNode( );

		void ResetTopologyUpdateFlag();
	private:
		int topology_update_flag;
};


#endif
