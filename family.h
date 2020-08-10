#ifndef FAMILY_h
#define FAMILY_h

#include <array>
#include <set>
#include <cmath>
#include "api/BamWriter.h"
#include "api/algorithms/Sort.h"
#include "ds.h"
#include "defines.h"

extern int MinBaseQual;
extern int MinMapQual;
extern int ConsensusUniqueCount;
extern int MinConSize;
extern bool CallMutations;
extern bool UseSingletonReads;
extern std::string UMIT;
extern std::string OUMI;

extern std::set<Site> germlineSet;
extern std::string refStr;
extern std::vector<std::string> intervalVector;
// extern std::vector<BamTools::BamAlignment> nonConsensusAlignmentVector;
extern std::vector<BamTools::BamAlignment> consensusAlignmentVector;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileup;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileupSS;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileupSingleton;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupCluster;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupSSCluster;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupSingletonCluster;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

class PairedRead {
private:
	void baseSeqToArray(std::string seq, int orient, int order, std::vector<BamTools::CigarOp> &cigarData, std::string qual, uint16_t mapq, std::string leftRefSeq, std::string rightRefSeq);
public:
	bool paired;
	bool marked;
	bool leftPair;
	bool rightPair;
	int firstOrient, secondOrient;
	uint16_t firstMapq, secondMapq;
	std::string readID;
	std::string firstBaseQual;
	std::string secondBaseQual;
	std::array<std::array<int,MAX_READ_LEN>,4> firstSeq;
	std::array<std::array<int,MAX_READ_LEN>,4> secondSeq;
	BamTools::BamAlignment firstAlignment;
	BamTools::BamAlignment secondAlignment;
	PairedRead *nextPair;

	PairedRead(bool pair, int orient, int order, std::string ID, std::string seq, std::vector<BamTools::CigarOp> &cigarData, std::string qual, uint16_t mapq, std::string leftRefSeq, std::string rightRefSeq, BamTools::BamAlignment al);

	PairedRead(const PairedRead &rp, int order);

	void addRead(int orient, int order, std::string seq, std::vector<BamTools::CigarOp> cigarData, std::string qual, uint16_t mapq, std::string leftRefSeq, std::string rightRefSeq, BamTools::BamAlignment al);
	void dumpTags(BamTools::BamAlignment &al, std::string barcode, int edit, int pair, std::string key, int ct);
};

class ReadFamily {
private:
	bool leftReads;
	bool rightReads;
	bool report;
	int count;
	int leftCount;
	int rightCount;
	int leftRef;
	int leftPos;
	int rightRef;
	int rightPos;
	int insertSize;
	int readLen;
	int leftEdit;
	int rightEdit;
	std::string key;
	std::string leftRefSeq;
	std::string rightRefSeq;
public:
	std::string leftSeq;
	std::string rightSeq;
	std::string leftMatch;
	std::string rightMatch;
	PairedRead *readPair, *lastPair;
	ReadFamily() {}
	ReadFamily(std::string key, PairedRead &rp, int lref, int lpos, int rref, int rpos, int order, int insertSize, std::string leftRefSeq, std::string rightRefSeq, int rl);
	std::string getKey();
	void getContig(int *ref, int *pos);
	bool checkReadPair(std::string ID);
	int missingReadPair();
	void addReadToPair(std::string ID, std::string seq, int orient, int order, std::vector<BamTools::CigarOp> cigarData, std::string qual, uint16_t mapq, BamTools::BamAlignment);
	void addNewReadPair(PairedRead &RP, int order);
	void callConsensus(int startPosition, int endPosition);
	bool checkBarcodeCollision();
	void matchConsensus(ReadFamily family);
	void addPileup(int, int);
	void addSingleStrandedPileup(int, int);
	void addSingletonPileup(int startPosition, int endPosition);
	void clearFamily();
	void printFamily();
	void getSingleStrandedCoverage(int leftPosition, int rightPosition);
	int getFamilySize() {return count;}
	int getFirstPairCount() {return leftCount;}
	int getSecondPairCount() {return rightCount;}
	//int dumpNonConsensusBAM();
	int dumpConsensusBAM();
};

void addPileupCounts(int pileupPos, char ls, char lrs, double E, int leftPos, int rightPos, int A, int C, int G, int T, int ind);
void addLOD(PairedRead *rpit, int pileupPos, char leftSeq, char leftRefSeq, int sit, int ind);

#endif

