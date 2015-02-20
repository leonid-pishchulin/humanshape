#ifndef _o_SaveCOLLADA_h
#define _o_SaveCOLLADA_h

#include "o_World.h"
#include "o_MetaCamera.h"
#include "o_TextureMap.h"
#include "o_SaveJPG.h"
#include <vector>
#include <string>
#include <map>

typedef std::map<o_Surface*,std::string> surfMap;
typedef std::map<o_Body*, surfMap> bodyMap; 

class o_SaveCOLLADA{
public:
	o_SaveCOLLADA(std::string name):texmapno(0){filename = name.substr(0,name.rfind('.'));}
	virtual ~o_SaveCOLLADA(){}

	int Save(std::vector<o_MetaCamera* >& list, float framerate, o_World& world);
	
protected:
	void SaveTextures(o_World& world, ofstream& out);
	void SaveTextures(o_Surface& s, ofstream& out, map<o_Surface*,string>& m);
	string SaveTextures(o_LocalCoordinateSystem& s, ofstream& out);

	std::string SaveTextureMap( o_TextureMap *map );
	std::string filename;
	std::map<o_Body*,std::map<o_Surface*,std::string> >surftexture;
	//std::map<o_LocalCoordinateSystem*,std::string> surfmaterial;

	int texmapno;
};

#endif
