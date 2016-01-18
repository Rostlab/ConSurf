// $Id: recognizeFormat.cpp 391 2005-06-12 12:44:57Z ninio $

#include "recognizeFormat.h"
#include "maseFormat.h"
#include "sequenceContainer.h"
#include "molphyFormat.h"
#include "phylipFormat.h"
#include "nexusFormat.h"
#include "fastaFormat.h"
#include "clustalFormat.h"
#include "nexusFormat.h"


sequenceContainer recognizeFormat::read(istream &infile, const alphabet* alph) {
	sequenceContainer mySeqData = readUnAligned(infile, alph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}

sequenceContainer recognizeFormat::readUnAligned(istream &infile, const alphabet* alph) {
	// recognize a format and returns the sequence container of it.
	sequenceContainer sc;
	if (!infile){
		string tmp = "error unable to open sequence input file ";
		errorMsg::reportError(tmp);
	}

	// this part eats spaces, tabs and such.
	char check = infile.peek();
	while ((check==' ') || (check == '\n') || (check == '\t')) {
		infile.get();
		check = infile.peek();
	}

	switch (check){
	case '#':
		sc=nexusFormat::readUnAligned(infile,alph);
		break;
	case '>':
		sc=fastaFormat::readUnAligned(infile,alph);
		break;
	case 'C':
		sc=clustalFormat::readUnAligned(infile,alph);
		break;
	case ';':
		sc=maseFormat::readUnAligned(infile,alph);
		break;	

	default:
		if (isdigit(check)){ 
			// here it can be either MOLPHY format of PHYLIP format
			// in PHYLIP format there are lines that are not empty, but the first 10 characters
			// are space.
			string s;
			getline(infile,s, '\n' ); // read the first line which are numbers in both formats
			getline(infile,s, '\n' ); // read the second line
			bool phylipFormat = false;
			int r = s.find_first_of(' ');
			if ((r==(s.size()-1)) || (r==-1)) phylipFormat = false;
			else phylipFormat = true;
			
			infile.seekg(0, ios::beg);
			if (phylipFormat == false) {
				sc=molphyFormat::readUnAligned(infile,alph);
			} else {
				sc=phylipFormat::readUnAligned(infile,alph);
			}
		}
		else{
			string line;
			getline(infile, line, '\n');
			string tmp2 = "The program can't recognise your format!";
			tmp2+="\nThis is the first line in your format:\n";
			tmp2+=line;
			errorMsg::reportError(tmp2);
		}
		break;
	}
	return sc;
}
