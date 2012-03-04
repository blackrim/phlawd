#ifndef _GENEDB_H
#define _GENEDB_H

#include <string>

using namespace std;

class GeneDB{
private:
    string name;//filename
    
public:
    GeneDB();
    GeneDB(string name);
    
    void initialize(bool overwrite);
};

#endif
