#ifndef MUTATIONS_h
#define MUTATIONS_h

#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <set>
#include <cmath>
#include "ds.h"
#include "defines.h"
#include "readUtils.h"
#include "family.h"

extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileup;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileupSS;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,5> pileupSingleton;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupCluster;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupSSCluster;
extern std::array<std::vector<int>,MAX_INTERVAL_LEN> pileupSingletonCluster;
extern std::array<std::array<int,MAX_INTERVAL_LEN>,9> strandArray;
extern std::array<int,MAX_INTERVAL_LEN> refPile;
extern std::array<int,MAX_INTERVAL_LEN> altPile;
extern std::array<double,MAX_INTERVAL_LEN> errPile;
extern std::array<double,MAX_INTERVAL_LEN> tumor_lod;
extern std::set<Site> dbSNPSet;

enum MUTATION_STATUS {NOT_MUTATION=0, ADAPTER_ARTIFACT=1, ILLUMINA_EDGE_ARTIFACT=2, MUTATION=64};

void addPileupCounts(int pileupPos, char ls, char lrs, double E, int leftPos, int rightPos, int A, int C, int G, int T, int ind);
void addLOD(PairedRead *rpit, int pileupPos, char leftSeq, char leftRefSeq, int sit, int ind);
void fillMutHeaders(std::string inputBamPrefix);
void call_mutations(std::string intRef, int leftPosition, int rightPosition, const char *refchars);
void closeMutFiles();

class Mutation {
private:
	bool indbSNP;
	char ref, alt;
	int position;
	int A, C, G, T;
	int fwdRef, revRef, fwdAlt, revAlt;
	int refDup, refSS, refSG;
	int altDup, altSS, altSG;
	int artifactCount;
	int mutDepth;
	int flag;
	double AFDup, AFSS, AFSG;
	double lod, ep, ratio;
	std::string refID;
	void calculateSOR();
public:
	Mutation() {}
	Mutation (char ref, std::string refID, int position, int A, int ssA, int singA, int C, int ssC, int singC, int G, int ssG, int singG, int T, int ssT, int singT, int artifactCountSS, int mutSiteSS, int fwdRef, int revRef, int fwdAlt, int revAlt, double fstar_lod, double, bool);
	void getMutationStatus();
	void printSite(std::ofstream &mutFile);
	void printVcf(std::ofstream &mutFile);
	void printIfMutation(std::ofstream &mutFile);
};

#endif

