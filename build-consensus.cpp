#include <iterator>
#include <unistd.h>
#include "family.h"
#include "utils.h"
#include "readUtils.h"
#include "faidx.h"
#include "defines.h"

using namespace std;
using namespace BamTools;

set<Site> dbSNPSet;
set<Site> germlineSet;
vector<Interval> coverageIntervals;
typedef map<string, class ReadFamily> ReadMap;
ReadMap readMap;
map<string, long int> chrMap;

bool germlineFlag = false;

void fillMutHeaders(std::string inputBamPrefix);
void call_mutations(std::string intRef, int leftPosition, int rightPosition, const char *refchars);
void closeMutFiles();

void debug_pileups()
{
	copy(pileupSS[0].begin(), pileupSS[0].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSS[1].begin(), pileupSS[1].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSS[2].begin(), pileupSS[2].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSS[3].begin(), pileupSS[3].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSingleton[0].begin(), pileupSingleton[0].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSingleton[1].begin(), pileupSingleton[1].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSingleton[2].begin(), pileupSingleton[2].end(), ostream_iterator<int>(cout,",")); cout << endl;
	copy(pileupSingleton[3].begin(), pileupSingleton[3].end(), ostream_iterator<int>(cout,",")); cout << endl;
}

void printCoverage(ofstream &coverageFile, string intervalRef, int intRef, int leftPosition, int rightPosition)
{
	static int intervalIterator;
	static int backupSinglePileup, backupDuplex, backupSingleton;
	//cout << "Enter coverage with interval " << intervalRef << ":" << leftPosition << "-" << rightPosition << endl;
	while(intervalIterator<coverageIntervals.size()) {
		Interval interval = coverageIntervals[intervalIterator];
		Interval window(chrRank(intervalRef), leftPosition, rightPosition);
		//cout << "Enter while loop with iterator value: " << intervalIterator << " and interval " << interval.chr << ":" << interval.start << "-" << interval.end << endl;
		if(interval.overlapsWith(window)) {
			//cout << "Enter overlaps\n";
			string target = intervalRef+":"+to_string((long long int)(interval.start+1))+"-"+to_string((long long int)(interval.end+1));
			int regionLength = interval.end - interval.start + 1;
			int totalDuplex=backupDuplex, totalSinglePileup=backupSinglePileup, totalSingleton=backupSingleton;
			//debug_pileups();
			for(int i=1; i<=regionLength; i++) {
				int sit = interval.start+1 -leftPosition + i;
				if(sit<0)
					continue;
				else if(sit>=10000)
					break;
				//cout << sit << " ";
				totalDuplex += (pileup[0][sit]+pileup[1][sit]+pileup[2][sit]+pileup[3][sit]+pileup[4][sit]);
				totalSinglePileup += (pileupSS[0][sit]+pileupSS[1][sit]+pileupSS[2][sit]+pileupSS[3][sit]+pileupSS[4][sit]);
				totalSingleton += (pileupSingleton[0][sit]+pileupSingleton[1][sit]+pileupSingleton[2][sit]+pileupSingleton[3][sit]+pileupSingleton[4][sit]);
			}
			cout << endl;
			if(interval.endsInside(window)) {
				coverageFile << target << "\t" << totalSinglePileup << "\t" << totalSinglePileup/regionLength << "\t" << totalDuplex << "\t" << totalDuplex/regionLength << "\t" << totalSingleton << "\t" << totalSingleton/regionLength << "\t" << (totalSinglePileup+totalSingleton)/regionLength << endl;
				intervalIterator++;
				backupSinglePileup = 0;
				backupDuplex = 0;
				backupSingleton = 0;
			} else {
				backupSinglePileup = totalSinglePileup;
				backupDuplex = totalDuplex;
				backupSingleton = totalSingleton;
				return;
			}
		} else {
			//cout << "Failed if condition\n";
			break;
		}
	}
	return;
}

bool get_alt_family(ReadFamily &family, ReadFamily &family_alt)
{
	string barcode_alt, barcode, key_alt, key;
	key = family.getKey();
	vector<string> keyVec = split(key,'-');

	if(keyVec.size()<4 || keyVec[0].empty() || keyVec[1].empty() || keyVec[2].empty() || keyVec[3].empty()) {
		cerr << "ERROR: Family with key:" << key << " seems invalid" << endl;
		exit(1);
	}

	barcode = keyVec[2]+"-"+keyVec[3];
	if(!generate_alt_barcode(barcode, barcode_alt, false)) {
#ifdef DEBUG
		cerr << "WARNING: Skipping Family: " << key << endl;
#endif
		return false;
	}
	key_alt = keyVec[1]+"-"+keyVec[0]+"-"+barcode_alt;

	if(readMap.find(key_alt)==readMap.end())
		return false;
	else
		family_alt = readMap[key_alt];
	return true;
}

bool getNextInterval(ifstream &intervalFile, int &intChr, int &intStart, int &intEnd)
{
	char delim = '\t';
	string intervalLine;

	while(getline(intervalFile, intervalLine)) {
		if(intervalLine[0] == '@')
			continue;
		else
			break;
	}

	vector<string> intervalRegion = split(intervalLine, delim);
	if(intervalRegion.size()<3 || intervalRegion[0].empty() || intervalRegion[1].empty() || intervalRegion[2].empty()) {
		cerr << "ERROR: Skipping invalid interval line " << intervalLine << endl;
		return false;
	}
	intChr = chrRank(intervalRegion[0]);
	intStart = stoi(intervalRegion[1]);
	intEnd = stoi(intervalRegion[2]);

	return true;
}

bool getNextdbSNPSite(ifstream &dbSNPFile, int &dbSNPChr, int &dbSNPPos)
{
	char delim = '\t';
	string dbSNPLine;

	while(getline(dbSNPFile,dbSNPLine)) {
		if(dbSNPLine[0] == '#')
			continue;
		else
			break;
	}

	vector<string> dbSNPRegion = split(dbSNPLine, delim);
	if(dbSNPRegion.size()<2 || dbSNPRegion[0].empty() || dbSNPRegion[1].empty()) {
		cerr << "ERROR: Skipping invalid dbSNP site " << dbSNPLine << endl;
		return false;
	}
	dbSNPChr = chrRank(dbSNPRegion[0]);
	dbSNPPos = stoi(dbSNPRegion[1]);

	return true;
}

bool readdbSNPFile(string dbSNPList, string intervalList)
{
	int intChr, intStart, intEnd, dbSNPChr, dbSNPPos;
	ifstream dbSNPFile(dbSNPList);
	ifstream intervalFile(intervalList);
	bool flag = true;

	if(!intervalFile.is_open()) {
		cerr << "ERROR: Could not open intervals list: " << intervalList << endl;
		return false;
	}
	if(!dbSNPFile.is_open()) {
		cerr << "ERROR: Could not open dbSNP file: " << dbSNPList << endl;
		return false;
	}

	if(getNextInterval(intervalFile, intChr, intStart, intEnd)==false) {
		cerr << "ERROR: Could not get next interval" << endl;
		flag = false;
	}
	if(getNextdbSNPSite(dbSNPFile, dbSNPChr, dbSNPPos)==false) {
		cerr << "ERROR: Could not get next dbSNP site" << endl;
		flag = false;
	}

	bool retStat;
	while(flag==true) {
		retStat = true;
		if( (intChr==dbSNPChr && intStart>dbSNPPos) || (intChr>dbSNPChr) ) {
			retStat = getNextdbSNPSite(dbSNPFile, dbSNPChr, dbSNPPos);
		} else if(intChr==dbSNPChr && intStart<=dbSNPPos && intEnd>=dbSNPPos) {
			Site site(dbSNPChr-1,dbSNPPos);
			dbSNPSet.insert(site);
			retStat = getNextdbSNPSite(dbSNPFile, dbSNPChr, dbSNPPos);
		} else if( (intChr==dbSNPChr&&intEnd<dbSNPPos) || dbSNPChr>intChr) {
			retStat = getNextInterval(intervalFile, intChr, intStart, intEnd);
		}
		if(retStat==false) {
			break;
		}
	}
	intervalFile.close();
	dbSNPFile.close();
	cout << "DBSNP set size = " << dbSNPSet.size() << endl;
	return true;
}

bool indexBamFile(string outputPath)
{
	BamReader reader;
	if (!reader.Open(outputPath))
		return false;
	reader.CreateIndex(BamIndex::STANDARD);
	reader.Close();
	return true;
}

// REVISIT: Output number of reads read and thrown out
bool getNextRead(BamReader &reader, BamAlignment &al)
{
	bool cigarFlag;
	while(true) {
		if(reader.GetNextAlignment(al)==false) {
			return false;
		}
		if(!al.IsMapped() && al.RefID==-1) {
			cout << "Read: " << al.Name << " onwards skipping" << endl;
			return false;
		}
#ifdef DEBUG
		char charName[100];
		strcpy(charName,al.Name.c_str());
		cout << "Debugging:" << charName << endl;
#endif
		if(!al.IsMapped() || al.MapQuality<MinMapQual || abs(al.InsertSize)>=MAX_INSERT_RANGE || al.RefID!=al.MateRefID) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " fails one of MQ/Insize checks. Skipping..\n";
#endif
			continue;
		}
		cigarFlag = true;
		for(vector<CigarOp>::iterator it=al.CigarData.begin(); it!=al.CigarData.end(); it++) {
			if((*it).Type=='S' || (*it).Type=='H') {
				cigarFlag = false;
				break;
			}
		}

		if(cigarFlag==false) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " contains clipped region. Skipping..\n";
#endif
			continue;
		}

		if(!al.IsPrimaryAlignment()) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " not primary alignment. Skipping..\n";
#endif
			continue;
		}
		if(al.IsReverseStrand()&&al.IsMateReverseStrand()) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " is reverse Strand pair. Skipping..\n";
#endif
			continue;
		}
		if(!al.IsReverseStrand()&&!al.IsMateReverseStrand()) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " is forward Strand pair. Skipping..\n";
#endif
			continue;
		}

		int8_t editdist = (int8_t)0;
		al.GetTag("NM",editdist);
		if(editdist-indelLen(al)>MAX_NM) {
#ifdef DEBUG
			cerr << "WARNING: Skipping read: " << al.Name << " with edit distance of " << +editdist << endl;
#endif
			continue;
		}
		if(al.BuildCharData()==false) {
			return false;
		}
// REVISIT: Soft-clipped and N bases next
		if(al.QueryBases.find('N')!=string::npos) {
#ifdef DEBUG
			cerr << "WARNING: Read: " << al.Name << " contains N chars... " << al.QueryBases << endl;
#endif
			continue;
		}
		return true;
	}
}

int maxRefCoord(string chr)
{
	if(chrMap.find(chr)==chrMap.end()) {
		cerr << "Chromosome Map does not contain chr: " << chr << ". Exiting..." << endl;
		exit(1);
	} else
		return chrMap[chr];
}

int prepAlignmentIntervals(RefVector& referenceSequences, vector<string> & bamChr, int *bamStart, int *bamEnd)
{
	int i = 0;
	for(vector<RefData>::iterator it=referenceSequences.begin(); it!=referenceSequences.end(); it++) {
		chrMap.insert(std::pair<string, long int> ((*it).RefName, (long int)(*it).RefLength));
		for(int pos=0; pos<(*it).RefLength; ) {
			int start = pos;
			int end = pos + 10000;
			if(end > (*it).RefLength)
				end = (*it).RefLength;
			bamChr[i] = (*it).RefName;
			bamStart[i] = start;
			bamEnd[i++] = end - 1;
			pos = end;
		}
	}
	return i;
}

void populateCoverageIntervals(ifstream &intervalFile)
{
	char delim = '\t';
	string intervalLine;

	while(getline(intervalFile, intervalLine)) {
		if(intervalLine[0] == '@')
			continue;
		vector<string> intervalRegion = split(intervalLine, delim);
		if(intervalRegion.size()<3 || intervalRegion[0].empty() || intervalRegion[1].empty() || intervalRegion[2].empty()) {
			cerr << "ERROR: Skipping invalid interval line " << intervalLine << endl;
			continue;
		}
		int intChr = chrRank(intervalRegion[0]);
		int intStart = stoi(intervalRegion[1])-1;
		int intEnd = stoi(intervalRegion[2])-1;
		Interval interval(intChr, intStart, intEnd);
		coverageIntervals.insert(coverageIntervals.end(), interval);
	}
	return;
}

bool buildFamilies(string inputBamFile, string outputPrefix, string intervalList, string refFasta)
{
	BamReader reader;
	faidx_t *fai;
	string intervalLine;
	ifstream intervalFile(intervalList);
	int maxIntervals = 500000;
	vector<string> bamChr(maxIntervals);
	int bamStart[maxIntervals], bamEnd[maxIntervals];

	if(!reader.Open(inputBamFile)) {
		cerr << "ERROR: Could not open BAM file:" << inputBamFile << endl;
		return false;
	}
	if(!reader.OpenIndex(inputBamFile+".bai")) {
		cerr << "ERROR: Could not open BAM index file:" << inputBamFile+".bai" << endl;
		return false;
	}
	if(!intervalFile.is_open()) {
		cerr << "ERROR: Could not open intervals list: " << intervalList << endl;
		return false;
	}
	populateCoverageIntervals(intervalFile);
	intervalFile.close();

	string consensusSortedBam(outputPrefix+".consensus.bam");
	//string nonConsensusSortedBam(inputBamFile+".nonConsensus.cleaned.bam");
	SamHeader header = reader.GetHeader();
	RefVector referenceSequences = reader.GetReferenceData();
	BamWriter ssWriterSorted; //, bcWriterSorted;
	if(!ssWriterSorted.Open(consensusSortedBam, header, referenceSequences)) {
		cerr << "ERROR: Could not open output bam file " << consensusSortedBam << endl;
		return false;
	}
/*
	if(!bcWriterSorted.Open(nonConsensusSortedBam, header, referenceSequences)) {
		cerr << "Could not open output bam file " << nonConsensusSortedBam << endl;
		return false;
	}
*/
	ofstream coverageFile(outputPrefix+".coverage.tsv", ofstream::out);
	coverageFile << "Target\tTotalSingleStrandCoverage\tAverageSingleStrandCoverage\tTotalDuplexCoverage\tAverageDuplexCoverage\tTotalSingletonCoverage\tAverageSingletonCoverage\tAverageCoverage\n";

	if(CallMutations==true) fillMutHeaders(outputPrefix);

	int intIdx = -0, bkIntRef = -1;
	int bamIntervalLength = prepAlignmentIntervals(referenceSequences, bamChr, bamStart, bamEnd);
	int intRef, leftPosition, rightPosition;
	do {
		char delim = '\t';
		intRef = reader.GetReferenceID(bamChr[intIdx]);
		leftPosition = bamStart[intIdx];
		rightPosition = bamEnd[intIdx];
		cout << "Processing interval: " << intRef << ":" << leftPosition << "-" << rightPosition << " now...\n";

		if(!reader.SetRegion(intRef, leftPosition, intRef, rightPosition)) {
#ifdef DEBUG
			cerr << "WARNING: Skipping interval: Could not Set region " << intRef << ":" << leftPosition << "-" << rightPosition << endl;
#endif
			continue;
		}

		BamAlignment al;
		int count = 0, len = 0;
		int leftInsert, rightInsert;
		leftInsert = leftPosition-MAX_INSERT_RANGE < 1 ? 1 : leftPosition-MAX_INSERT_RANGE;
		rightInsert = rightPosition+MAX_INSERT_RANGE > maxRefCoord(bamChr[intIdx])-1 ? maxRefCoord(bamChr[intIdx])-1 : rightPosition+MAX_INSERT_RANGE;
		string strreg = bamChr[intIdx] + ":" + to_string((long long int)leftInsert) + "-" + to_string((long long int)rightInsert);
		string overlapReg = to_string((long long int)intRef) + ":" + to_string((long long int)leftPosition) + "-" + to_string((long long int)rightPosition);
		len = (rightInsert) - (leftInsert) + 1;
		fai  = fai_load((refFasta.c_str()));
		if(fai==NULL) {
			cerr << "ERROR: Failed to load fasta.fai file: " << refFasta << endl;
			return false;
		}
		char *refchars = fai_fetch(fai, strreg.c_str(), &len);
		refStr = refchars;

		if(bkIntRef!=intRef) {
			intervalVector.clear();
			bkIntRef = intRef;
		}

		while (true) {
			BamAlignment bkal;
			bkal = al;

			string barcode, barcode_alt;
			string leftID, leftPos, rightID, rightPos;
			string leftRefSeq, rightRefSeq;

			if(!getNextRead(reader, al)) {
				cout << "No more reads in interval: " << intRef << ":" << leftPosition << "-" << rightPosition << endl;
				cout << "Last read read: " << bkal.Name << " " << bkal.Position << endl;
				break;
			}
			bool pair = al.IsPaired();
			int order = al.IsFirstMate() ? 1 : 2;
			int orientation = al.IsReverseStrand() ? 2 : 1;
			int readLen = al.Length;

			al.GetTag(UMIT,barcode);
#ifdef DEBUG
			cout << "BARCODE: " << barcode << endl;
#endif
			if(!generate_alt_barcode(barcode, barcode_alt, true)) {
#ifdef DEBUG
				cerr << "WARNING: Duplex Pair Barcode not generated. Skipping Read: " << al.Name << endl;
#endif
				continue;
			}
			int startPos = al.Position;
			int insSize = al.InsertSize;
			if(insSize>0)
				insSize += (inslen(al));
			else if(insSize<0)
				insSize -= (inslen(al));
			leftID = to_string((long long int)al.RefID);
			rightID = to_string((long long int)al.MateRefID);
			if(order==1) {
				if(orientation==1) {
					leftPos = to_string((long long int)startPos);
					rightPos = to_string((long long int)(startPos+insSize-1));
					if((startPos+insSize-al.Length)-(leftInsert)<0 || (startPos+insSize-al.Length)-(leftInsert)>len)
						continue;
					rightRefSeq = refStr.substr((startPos+insSize-al.Length)-(leftInsert)+1, readLen*2);
				} else if(orientation==2) {
					leftPos = to_string((long long int)(startPos+al.Length-1));
					rightPos = to_string((long long int)startPos+insSize+al.Length);
					if((startPos+insSize+al.Length)-(leftInsert)<0 || (startPos+insSize+al.Length)-(leftInsert)>len)
						continue;
					if(startPos<=al.MatePosition)
						continue;
					rightRefSeq = refStr.substr((startPos+insSize+al.Length)-(leftInsert)+1, readLen*2);
				} else {
					cerr << "ERROR: In order 1, unexpected orientation value: " << orientation << endl;
					return false;
				}
				leftRefSeq = refStr.substr(startPos-(leftInsert)+1, readLen*2);
			} else if(order==2){
				if(orientation==1) {
					leftPos = to_string((long long int)startPos+insSize-1);
					rightPos = to_string((long long int)startPos);
					if((startPos+insSize-al.Length)-(leftInsert)<0 || (startPos+insSize-al.Length)-(leftInsert)>len)
						continue;
					if(al.MatePosition<=startPos)
						continue;
					leftRefSeq = refStr.substr((startPos+insSize-al.Length)-(leftInsert)+1, readLen*2);
				} else if(orientation==2) {
					leftPos = to_string((long long int)startPos+insSize+al.Length);
					rightPos = to_string((long long int)(startPos+al.Length-1));
					if((startPos+insSize+al.Length)-(leftInsert)<0 || (startPos+insSize+al.Length)-(leftInsert)>len)
						continue;
					leftRefSeq = refStr.substr((startPos+insSize+al.Length)-(leftInsert)+1, readLen*2);
				} else {
					cerr << "ERROR: In order 2, unexpected orientation value: " << orientation << endl;
					return false;
				}
				rightRefSeq = refStr.substr(startPos-(leftInsert)+1, readLen*2);
			} else {
				cerr << "ERROR: Unexpected order value: " << order << endl;
				return false;
			}
			string key = leftID+":"+leftPos+"-"+rightID+":"+rightPos+"-"+barcode;
			string key_alt = rightID+":"+rightPos+"-"+leftID+":"+leftPos+"-"+barcode_alt;
			PairedRead rp(pair, orientation, order, al.Name, al.QueryBases, al.CigarData, al.Qualities, al.MapQuality, leftRefSeq, rightRefSeq, al);
#ifdef DEBUG
			char keyStr[100], keyAltStr[100];
			strcpy(keyStr,key.c_str());
			strcpy(keyAltStr,key_alt.c_str());
#endif

			// Build consensus for each pair of duplex independently
			// Traverse readMap and match duplex families
			if(readMap.find(key)==readMap.end()) {
				ReadFamily family(key, rp, stoi(leftID), stoi(leftPos), stoi(rightID), stoi(rightPos), order, abs(al.InsertSize), leftRefSeq, rightRefSeq, readLen);
				readMap.insert(std::pair<string,ReadFamily> (key,family));
			} else {
				if(readMap[key].checkReadPair(al.Name))
					readMap[key].addReadToPair(al.Name, al.QueryBases, orientation, order, al.CigarData, al.Qualities, al.MapQuality, al);
				else
					readMap[key].addNewReadPair(rp, order);
			}
		}

		// tumor_lod.fill(0.0); errPile.fill(0.0); refPile.fill(0); altPile.fill(0);
		for(int ait=0; ait<MAX_INTERVAL_LEN; ait++) {
			pileupCluster[ait].clear();
			pileupSSCluster[ait].clear();
			pileupSingletonCluster[ait].clear();
		}
		for(int ait=0; ait<5; ait++) {
			pileup[ait].fill(0);
			pileupSS[ait].fill(0);
			pileupSingleton[ait].fill(0);
		}
		//for(int ait=0; ait<9; ait++) { strandArray[ait].fill(0); }

		map<string,ReadFamily>::iterator bkit;
		for (map<string,ReadFamily>::iterator it=readMap.begin(); it!=readMap.end(); ) {
			if(germlineFlag==false)
				it->second.callConsensus(leftPosition, rightPosition);
			else if(it->second.checkBarcodeCollision()==true)
				it->second.callConsensus(leftPosition, rightPosition);
			//it->second.printFamily();
			it++;
		}

		for (map<string,ReadFamily>::iterator it=readMap.begin(); it!=readMap.end(); ) {
			ReadFamily family_alt;
			if(get_alt_family(it->second, family_alt)==true) {
				it->second.matchConsensus(family_alt);
				it->second.addPileup(leftPosition, rightPosition);
			}
			//if(it->second.getFamilySize()>=2) {
			it->second.addSingleStrandedPileup(leftPosition, rightPosition);
			//} else {
			if(UseSingletonReads)
				it->second.addSingletonPileup(leftPosition, rightPosition);
			//}
			//it->second.dumpNonConsensusBAM();
			it->second.dumpConsensusBAM();
			it->second.clearFamily();
			bkit = it;
			it++;
			readMap.erase(bkit);
		}

		intervalVector.push_back(overlapReg);
		printCoverage(coverageFile, bamChr[intIdx], intRef, leftPosition, rightPosition);
		sort(consensusAlignmentVector.begin(), consensusAlignmentVector.end(), Algorithms::Sort::ByPosition());
		for(vector<BamAlignment>::iterator al=consensusAlignmentVector.begin(); al!=consensusAlignmentVector.end(); al++) {
			ssWriterSorted.SaveAlignment(*al);
		}
		consensusAlignmentVector.clear();
/*
		sort(nonConsensusAlignmentVector.begin(), nonConsensusAlignmentVector.end(), Algorithms::Sort::ByPosition());
		for(vector<BamAlignment>::iterator al=nonConsensusAlignmentVector.begin(); al!=nonConsensusAlignmentVector.end(); al++) {
			bcWriterSorted.SaveAlignment(*al);
		}
		nonConsensusAlignmentVector.clear();
*/
		if(CallMutations==true) call_mutations(bamChr[intIdx], leftPosition, rightPosition, refchars);

#ifdef DEBUG
		for(int refct=0; refct<strlen(refchars)&&refct<(rightPosition+1-leftPosition); refct++)
			cout << refchars[MAX_INSERT_RANGE+refct] << ",";
		cout << endl;
		for(int refct=0; refct<strlen(refchars)&&refct<(rightPosition+1-leftPosition); refct++)
			cout << leftPosition+refct << ",";
		cout << endl;
		debug_pileups();
#endif
		free(refchars);
		if(fai!=NULL)
			fai_destroy(fai);
	} while(++intIdx < bamIntervalLength);

	cout << "No more reads found in bam file: " << inputBamFile  << endl;

	reader.Close();
	ssWriterSorted.Close();
	//bcWriterSorted.Close();
	coverageFile.close();

	if(CallMutations==true) closeMutFiles();
/*
	if(indexBamFile(nonConsensusSortedBam) == false) {
		cerr << "Bam file indexing failed for: " << nonConsensusSortedBam << endl;
		return 1;
	}
*/
	if(indexBamFile(consensusSortedBam) == false) {
		cerr << "ERROR: Bam file indexing failed for: " << consensusSortedBam << endl;
		return 1;
	}
	return true;
}

int main(int argc, char* argv[])
{
	string inputBamFile, outputPrefix, intervalList, refFasta, dbSNPFile, germlineFile, umiTag, useIntervalList;
	int min_bq, min_mq, max_umi, con_unique, min_con_size, padding_length = 0;
	extern char *optarg;
	extern int optind;
	int c, err = 0;
	bool pad_intervals = false;
	static char usage[] = " -i <input_bam> -o <output_prefix> -l <interval_list> -r <ref_fasta> -d <dbsnp_file>\n [-b min_base_qual] [-m min_map_qual] [-U UMI Tag] [-e max_umi_edits] [-n consensus_unique_count] [-c min_consensus_size] [-s use_singleton_reads] [-p padding_length] [-v]\n\nParameter Description:\n\n -i Path to input bam file\n -o Prefix to output files including directory path (Will be used for consensus bam and coverage file)\n -l Interval file in bed format\n -r Path to reference fasta file\n -d Path to dbSNP file\n -b Minimum base quality to consider for consensus building [Default 20]\n -m Minimum mapping quality to consider for consensus building [Default 30]\n -U UMI Tag name [Default BC]\n -e Max edits permitted in UMI for consensus building [Default 2]\n -n Number of consensus reads in a family to be marked as unique. Used to downweight singleton reads [Default 1]\n -c Minimum number of reads in a family to call consensus [Default 2]\n -s Use singleton reads towards coverage and variant calling [Default true]\n -v Call Variants from consensus building. Switched off for now [Default false]\n -p Pad intervals for mutation calling [Default 0]\n\n";

	while ((c = getopt(argc, argv, "i:o:l:r:d:b:m:u:e:n:c:s:vp:")) != -1) {
		switch(c) {
		case 'i':
			inputBamFile = optarg;
			break;
		case 'o':
			outputPrefix = optarg;
			break;
		case 'l':
			intervalList = optarg;
			break;
		case 'r':
			refFasta = optarg;
			break;
		case 'd':
			dbSNPFile = optarg;
			break;
/*
		case 'g':
			germlineFlag = true;
			germlineFile = optarg;
			break;
*/
		case 'b':
			min_bq = atoi(optarg);
			if(min_bq < 0 || min_bq > MAX_POSSIBLE_BASE_QUAL)
				cerr << "Invalid minimum base quality of " << min_bq << ". Ignoring..." << endl;
			else
				MinBaseQual = min_bq;
			break;
		case 'm':
			min_mq = atoi(optarg);
			if(min_mq < 0 || min_mq > MAX_POSSIBLE_MAP_QUAL)
				cerr << "Invalid minimum mapping quality of " << min_mq << ". Ignoring..." << endl;
			else
				MinMapQual = min_mq;
			break;
		case 'e':
			max_umi = atoi(optarg);
			if(max_umi < 0 || max_umi > MAX_FEASIBLE_UMI_EDIT_DISTANCE)
				cerr << "Invalid maximum umi edit distance of " << max_umi << ". Ignoring..." << endl;
			else
				MaxUMIEditDist = max_umi;
			break;
		case 'n':
			con_unique = atoi(optarg);
			if(con_unique < 1)
				cerr << "Invalid consensus family unique count of " << con_unique << ". Ignoring..." << endl;
			else
				ConsensusUniqueCount = con_unique;
			break;
		case 's':
			if(strcmp(optarg,"true")!=0 && strcmp(optarg,"false")!=0)
				cerr << "Incorrect assignment. Please specify either 'true' or 'false'. Ignoring..." << endl;
			else
				istringstream(optarg) >> std::boolalpha >> UseSingletonReads;
			break;
		case 'c':
			min_con_size = atoi(optarg);
			if(min_con_size < MIN_POSSIBLE_CON_SIZE)
				cerr << "Invalid minimum consensus family size of " << min_con_size << ". Ignoring..." << endl;
			else
				MinConSize = min_con_size;
			break;
		case 'U':
			umiTag = optarg;
			UMIT = umiTag;
			break;
		case 'v':
			cout << "Unoperational option: -v (call_mutations). Skipping for now..." << endl;
			CallMutations = false;
			break;
		case 'p':
			cout << "Getting padding length" << endl;
			padding_length = atoi(optarg);
			if(padding_length < 0 || padding_length > MAX_POSSIBLE_PAD_LENGTH) {
				cerr << "Invalid padding length of " << padding_length << ". Ignoring..." << endl;
				padding_length = 0;
			}
			break;
		case '?':
			err = true;
			break;
		}
	}
	if(inputBamFile.empty() || outputPrefix.empty() || intervalList.empty() || refFasta.empty() || dbSNPFile.empty()) {
		cerr << argv[0] << ": missing command-line arguments" << endl;
		cerr << "Usage: " << argv[0] <<  usage;
		exit(1);
	} else if(err==true || optind < argc) {
		cerr << argv[0] << ": Ignoring unrecognized command-line parameters\n";
		cerr << "Usage: " << argv[0] <<  usage;
	}

	time_t start_time;
	time(&start_time);
	cout << "Start time: " << start_time << endl;

	if(verifyCoordinates(refFasta, inputBamFile)==false) {
		cerr << "ERROR: Inconsistent coordinates between reference/bam/interval list files. Exiting..." << endl;
		return 1;
	}
	cout << "Coordinates verified...\n";

	if(padding_length > 0) {
		string paddedIntervalList = getPaddedIntervalList(intervalList, outputPrefix, padding_length);
		useIntervalList = paddedIntervalList;
	} else {
		useIntervalList = intervalList;
	}


	if(readdbSNPFile(dbSNPFile, useIntervalList)==false) {
		cerr << "ERROR: dbSNP read failed" << endl;
		return 1;
	}

	cout << "Finished reading dbSNP file...\n";
/*
	if(germlineFlag==true) {
		if(readGermlineFile(germlineFile)==false) {
			cerr << "WARNING: Failed to read germline file. Ignoring..." << endl;
		}
		cout << "Finished reading germline file...\n";
	}
*/
	if(buildFamilies(inputBamFile, outputPrefix, useIntervalList, refFasta)==false) {
		cerr << "ERROR: Failed to build consensus families for snv analysis" << endl;
		return 1;
	}
	cout << "Finished consensus building...\n";

	time_t end_time;
	time(&end_time);
	cout << "End time: " << end_time << endl;
	cout << "Total time (seconds): " << difftime(end_time, start_time) << endl;
	cout << "Consensus building completed successfully. Exiting now..." << endl;
	return 0;
}

