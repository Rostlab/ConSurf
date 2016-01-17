// $Id: nexusFormat.cpp 391 2005-06-12 12:44:57Z ninio $

#include "nexusFormat.h"
#include "someUtil.h"
#include "errorMsg.h"

sequenceContainer nexusFormat::read(istream &infile, const alphabet* pAlph) {
	sequenceContainer mySeqData = readUnAligned(infile, pAlph);
	mySeqData.makeSureAllSeqAreSameLengthAndGetLen();
	return mySeqData;
}

sequenceContainer nexusFormat::readUnAligned(istream &infile, const alphabet* pAlph) {
	if (!infile) {
		errorMsg::reportError("unable to read mase format, could not open file");
	}
	sequenceContainer mySeqData;;

	vector<string> seqFileData;
	putFileIntoVectorStringArray(infile,seqFileData);

	vector<string>::const_iterator it1 = seqFileData.begin();
	// make sure that the first 6 chars in the first line is #NEXUS
	if (it1->size()<6) errorMsg::reportError("first word in a nexus sequence file format must be #NEXUS",1);
	if (  ((*it1)[0] != '#') 
		|| (((*it1)[1] != 'N') && ((*it1)[1] != 'n'))
		|| (((*it1)[2] != 'E') && ((*it1)[2] != 'e'))
		|| (((*it1)[3] != 'X') && ((*it1)[3] != 'x'))
		|| (((*it1)[4] != 'U') && ((*it1)[4] != 'u'))
		|| (((*it1)[5] != 'S') && ((*it1)[5] != 's')) ) {
		errorMsg::reportError("first word in a nexus sequence file format must be #NEXUS",1);
	}
	it1++;

	while ( ( (*it1).find("matrix")  == -1) && (it1!= seqFileData.end()))
	{ //check for the word matrix
		++it1;
	}
	
	int localid=0;
	if ((*it1).find("matrix") != -1)
	{

		for (++it1; it1 != seqFileData.end() ; ++it1)
		{
			if ((*it1).find("end;") != -1)
				break;
			if (it1->empty() || ((*it1).find(';') != -1)) 
			{ // empty line constinue
				continue;
			}
			sequence seq(pAlph);
			
			string taxonName;
			string remark;
			string stringSeq;
			bool beforeName = true;
			string::const_iterator stringIt = (it1)->begin();
			for (; stringIt != (it1)->end(); ++stringIt)
			{ //first loop finds the taxon name
				if ( ((*stringIt) == ' ') || ((*stringIt) == '\t'))
					if (beforeName == true)
						continue; //spaces before taxon name are legal
					else
						break; //A space marks the end of the taxon name 
				else 
				{
					taxonName += (*stringIt);
					beforeName = false;
				}
			}

			for (; stringIt != (it1)->end(); ++stringIt) 
			{//second loop finds the sequecne
				if ( ((*stringIt)==' ') || 	((*stringIt) == '\t'))
					continue;
				else stringSeq += (*stringIt);
			}

			mySeqData.add(sequence(stringSeq, taxonName, remark, localid, pAlph));
			localid++;
		}
	}
	else
	{
		errorMsg::reportError("no sequence data in nexus file - no matrix keyword found");
	}
	
	return mySeqData;
}

void nexusFormat::write(ostream &out, const sequenceContainer& sc) {
	//vector<string> gfr = sd.getGeneralFileRemarks();
	//if (gfr.empty()) out<<";;\n;;\n";
	//for (vector<string>::const_iterator k=gfr.begin() ; k !=  gfr.end() ; ++k )
	//	out<<(*k)<<endl;
	out<<"#NEXUS"<<endl;
	out<<"begin data;"<<endl;
	out<<"dimensions ntax="<<sc.numberOfSeqs()<<" nchar="<<sc.seqLen() <<";"<<endl;
	if (sc.alphabetSize() == 4) 
		out<<"format datatype=dna gap=-;"<<endl;
	else
		out<<"format datatype=protein gap=-;"<<endl;
	out<<"matrix"<<endl;

	for (sequenceContainer::constTaxaIterator itSeq=sc.constTaxaBegin();itSeq!=sc.constTaxaEnd();++itSeq) {
		out<<"\t"<<itSeq->name()<<"\t"<<itSeq->toString()<<endl;
	}
	out<<";"<<endl;
	out<<"end;"<<endl;
}

