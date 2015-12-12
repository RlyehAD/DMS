#include "dms_headers.h"

bool readAngles(std::vector< std::vector<PetscInt> >& indices, std::fstream& fp) {

        PetscFunctionBegin;
        std::string line;

        bool fAngles = false;
        bool success = false;

	PetscInt index[2];

        while( std::getline(fp, line) ) {

                std::size_t found = line.find("angles");

                if (found != std::string::npos) {
                        fAngles = true;
                        success	= true;
			continue;
                }

		if(fAngles)
                        if(line.empty())
                                fAngles = false;

                if(fAngles)
                        if( line.find(";") == std::string::npos ) {

                                std::istringstream iss(line);

                                std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                        std::istream_iterator<std::string>{}};


				if(tokens.size() > 0) {
					std::istringstream ( tokens[0] ) >> index[0];
					std::istringstream ( tokens[2] ) >> index[1];

					index[0] -= 1;
                                        index[1] -= 1;

					std::vector<PetscInt> tmp(index, index + 2);

                                	indices.push_back(tmp);
				}
                        }                          	
        }

        PetscFunctionReturn(success);
}

bool readBonds(std::vector< std::vector<PetscInt> >& indices, std::fstream& fp) {

	PetscFunctionBegin;
	std::string line;

	bool fBonds = false;
	bool success = false;

	PetscInt index[2];

	while( std::getline(fp, line) ) {
	
		std::size_t found = line.find("bonds");

		if (found != std::string::npos) {
			fBonds = true;
			success = true;
			continue;
		}

		if(fBonds)
                        if(line.empty())
                                fBonds = false;

		if(fBonds)
			if( line.find(";") == std::string::npos ) {

				std::istringstream iss(line);

				std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                      			std::istream_iterator<std::string>{}};

				if(tokens.size() > 0) {
                                        std::istringstream ( tokens[0] ) >> index[0];
                                        std::istringstream ( tokens[1] ) >> index[1];

					index[0] -= 1;
					index[1] -= 1;

                                        std::vector<PetscInt> tmp(index, index + 2);

                                        indices.push_back(tmp);
                                }

			}				
	}

	PetscFunctionReturn(success);
}

bool vecSort (std::vector<PetscInt> i, std::vector<PetscInt> j) { 
	return (i[0] < j[0]); 
}

std::vector< std::vector<PetscInt> > dmsReadTop(char* fname2) {

        std::fstream fp;
	std::vector<char*> fname;
	fname.push_back("TOP/topol_Protein_chain_A.itp");
	fname.push_back("TOP/topol_Protein_chain_B.itp");
	fname.push_back("TOP/topol_Protein_chain_C.itp");
	fname.push_back("TOP/topol_Protein_chain_D.itp");
	fname.push_back("TOP/topol_Protein_chain_E.itp");
	fname.push_back("TOP/topol_Protein_chain_F.itp");
	fname.push_back("TOP/topol_Protein_chain_G.itp");
	fname.push_back("TOP/topol_Protein_chain_H.itp");
	fname.push_back("TOP/topol_Protein_chain_I.itp");
	fname.push_back("TOP/topol_Protein_chain_J.itp");
	fname.push_back("TOP/topol_Protein_chain_K.itp");
	fname.push_back("TOP/topol_Protein_chain_L.itp");

        std::vector< std::vector<PetscInt> > indices;

	for(auto i = 0; i < 10; i++) {

	        fp.open(fname[i], std::ios::in);

        	readBonds(indices, fp);

		fp.clear();
		fp.seekg(0, std::ios::beg);
        	readAngles(indices, fp);

		fp.close();
	}

	std::sort (indices.begin(), indices.end(), vecSort); 
	// sorting makes memory misses fewer when assembling the matrices

	auto nCons = indices.size();
	auto nAtoms = 28745;
	indices.resize(nCons);

	for(auto i = nCons; i < 1 * nCons; i++) {
		indices[i].push_back(indices[i - nCons][0] + nAtoms); 
		indices[i].push_back(indices[i - nCons][1] + nAtoms);
	}

        PetscFunctionReturn(indices);
}

