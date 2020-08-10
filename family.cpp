#include "family.h"

using namespace std;
using namespace BamTools;

int MinBaseQual = 20;
int MinMapQual = 30;
int ConsensusUniqueCount = 1;
int MinConSize = 2;
bool CallMutations = false;
bool UseSingletonReads = true;
string UMIT = "BC";
string OUMI = "OB";

string refStr;
vector<string> intervalVector;
vector<BamAlignment> consensusAlignmentVector;
//vector<BamAlignment> nonConsensusAlignmentVector;
// These arrays store unique coverage values
array<array<int,MAX_INTERVAL_LEN>,5> pileup;
array<array<int,MAX_INTERVAL_LEN>,5> pileupSS;
array<array<int,MAX_INTERVAL_LEN>,5> pileupSingleton;
// These arrays store position of mutation for clustering-based errors
array<vector<int>,MAX_INTERVAL_LEN> pileupCluster;
array<vector<int>,MAX_INTERVAL_LEN> pileupSSCluster;
array<vector<int>,MAX_INTERVAL_LEN> pileupSingletonCluster;

bool intervalOverlap(string intv, int refID, int start, int end)
{
	vector<string> intervalRegion = split(intv, ':');
	vector<string> coords = split(intervalRegion[1],'-');
	
	if(intervalRegion[0].empty() || coords[0].empty() || coords[1].empty()) {
		cerr << "ERROR: Invalid interval line provided for overlap check" << intv << endl;
		return false;
	}
	int intRef = stoi(intervalRegion[0]);
	int leftPosition = stoi(coords[0]);
	int rightPosition = stoi(coords[1]);

	if(intRef!=refID || (start<leftPosition&&end<leftPosition) || (start>rightPosition&&end>rightPosition))
		return false;
	else
		return true;
}

int distanceToInterval(string intv, int start)
{
	vector<string> intervalRegion = split(intv, ':');
	vector<string> coords = split(intervalRegion[1],'-');
	
	if(intervalRegion[0].empty() || coords[0].empty() || coords[1].empty()) {
		cerr << "ERROR: Invalid interval line provided for overlap check" << intv << endl;
		return false;
	}
	int leftPosition = stoi(coords[0]);
	int rightPosition = stoi(coords[1]);

	return start-rightPosition+1;
}

bool priorOverlap(int refID, int start, int end, int readLen)
{
	for(vector<string>::reverse_iterator intv=intervalVector.rbegin(); intv!=intervalVector.rend(); intv++) {
		if(distanceToInterval(*intv, start)>readLen)
			return false;
		if(intervalOverlap(*intv, refID, start, end)) {
			return true;
		}
	}
	return false;
}

PairedRead::PairedRead(bool pair, int orient, int order, string ID, string seq, vector<CigarOp> &cigarData, string qual, uint16_t mapq, string leftRefSeq, string rightRefSeq, BamAlignment al)
{
	paired = pair;
	marked = paired == true ? false : true;
	leftPair = false;
	rightPair = false;
	firstOrient = secondOrient = 0;
	firstMapq = secondMapq = 0;
	readID = ID;
	for(int ait=0; ait<4; ait++) {
		firstSeq[ait].fill(0);
		secondSeq[ait].fill(0);
	}
	baseSeqToArray(seq, orient, order, cigarData, qual, mapq, leftRefSeq, rightRefSeq);
	if(order==1)
		firstAlignment = al;
	else if(order==2)
		secondAlignment = al;
	else
		cerr << "ERROR: Incorrect segment number " << order << " for read " << al.Name << endl;
	nextPair = NULL;
}

PairedRead::PairedRead(const PairedRead &rp, int order)
{
	paired = rp.paired;
	marked = paired == true ? false : true;
	firstOrient = rp.firstOrient;
	secondOrient = rp.secondOrient;
	leftPair = rp.leftPair;
	rightPair = rp.rightPair;
	readID = rp.readID;
	firstBaseQual = rp.firstBaseQual;
	secondBaseQual = rp.secondBaseQual;
	firstMapq = rp.firstMapq;
	secondMapq = rp.secondMapq;
	firstSeq = rp.firstSeq;
	secondSeq = rp.secondSeq;
	firstAlignment = rp.firstAlignment;
	secondAlignment = rp.secondAlignment;
	nextPair = NULL;
}

void PairedRead::baseSeqToArray(string seq, int orient, int order, vector<CigarOp> &cigarData, string qual, uint16_t mapq, string leftRefSeq, string rightRefSeq)
{
	int itSeq = 0, i = 0, readSum = 0;
	array <int, MAX_READ_LEN> insString, delString;

	//cigarString.fill(0);
	delString.fill(0);
	insString.fill(0);

	for(vector<CigarOp>::iterator it=cigarData.begin(); it!=cigarData.end(); it++) {
		if((*it).Type == 'T') // REVISIT: Should probably be removing the readSum check
			readSum += (*it).Length;
		if((*it).Type == 'M') {
			for(int j=0; i<seq.length()-readSum && j<(*it).Length; i++,j++) {
			//	cigarString[i] = 1;
			}
		} else if((*it).Type == 'I') {
			for(int j=0; i<seq.length()-readSum && j<(*it).Length; i++,j++) {
			//	cigarString[i] = 1;
				insString[i] = 1;
			}
		} else if((*it).Type == 'D') {
			for(int j=0; i<seq.length() && j<(*it).Length; i++,j++) {
			//	cigarString[i] = 1;
				delString[i] = 1;
			}
		}
	}
	char base;
	int delSum, insSum;
	if(order==1) {
		firstOrient = orient;
		firstMapq = mapq;
		delSum = 0;
		insSum = 0;
		for(int sit=0; sit<seq.length(); sit++) {
			//REVISIT: Probably put back this qual check for QC
			if(insString[sit]==1) {insSum += 1; continue;}
			if(delString[sit]==1) {delSum += 1; firstBaseQual += '#'; continue;}
			firstBaseQual += qual[sit-delSum];
			base = seq[sit-delSum];
			//base = cigarString[sit]==1 ? seq[sit-delSum] : leftRefSeq[sit-insSum];
			switch(base) {
				case 'A': firstSeq[0][sit-insSum]++;break;
				case 'C': firstSeq[1][sit-insSum]++;break;
				case 'G': firstSeq[2][sit-insSum]++;break;
				case 'T': firstSeq[3][sit-insSum]++;break;
				default: break;
			}
		}
		leftPair = true;
	} else if(order==2) {
		secondOrient = orient;
		secondMapq = mapq;
		delSum = 0;
		insSum = 0;
		for(int sit = 0; sit<seq.length(); sit++) {
			if(insString[sit]==1) {insSum += 1; continue;}
			if(delString[sit]==1) {delSum += 1; secondBaseQual += '#'; continue;}
			secondBaseQual += qual[sit-delSum];
			base = seq[sit-delSum];
			//base = cigarString[sit]==1 ? seq[sit-delSum] : leftRefSeq[sit-insSum];
			switch(base) {
				case 'A': secondSeq[0][sit-insSum]++;break;
				case 'C': secondSeq[1][sit-insSum]++;break;
				case 'G': secondSeq[2][sit-insSum]++;break;
				case 'T': secondSeq[3][sit-insSum]++;break;
				default: break;
			}
		}
		rightPair = true;
	}
	return;
}

void PairedRead::addRead(int orient, int order, string seq, vector<CigarOp> cigarData, string qual, uint16_t mapq, string leftRefSeq, string rightRefSeq, BamAlignment al)
{
	baseSeqToArray(seq, orient, order, cigarData, qual, mapq, leftRefSeq, rightRefSeq);
	if(order==1)
		firstAlignment = al;
	else if(order==2)
		secondAlignment = al;
	else
		cerr << "ERROR: Incorrect segment number " << order << " for read " << al.Name << endl;
	marked = true;
	return;
}

void PairedRead::dumpTags(BamAlignment &al, string barcode, int edit, int pair, string key, int ct)
{
	string tagStr;
	int8_t tagVal = (int8_t)0;
	BamAlignment baseAl;
	vector<string> intTagList = {"OP"};
	vector<string> numTagList = {"NM", "MQ", "AS", "XS"};
	vector<string> stringTagList = {UMIT, "MD", "RG", "OQ", "OC", "XP"};

	if(pair==1)
		baseAl = firstAlignment;
	else if(pair==2)
		baseAl = secondAlignment;
	for(vector<string>::iterator tag=stringTagList.begin(); tag!=stringTagList.end(); tag++) {
		if(!baseAl.HasTag(*tag)) {
			if(*tag=="RG") {
				cerr << "ERROR: RG Tag missing from read: " << readID << ". Exiting..." << endl;
				exit(1);
			}
			if(*tag==UMIT||barcode.empty()) {
				cerr << "ERROR: UMI Tag missing from read: " << readID << ". Exiting..." << endl;
				exit(1);
			}
		} else {
			if(*tag==UMIT) {
				baseAl.GetTag(*tag,tagStr);
				al.AddTag(OUMI,string("Z"),tagStr);
				al.AddTag(UMIT,string("Z"),barcode);
				continue;
			}
			baseAl.GetTag(*tag,tagStr);
			al.AddTag(*tag,string("Z"),tagStr);
		}
	}
	al.AddTag("XK", string("Z"), key);
	for(vector<string>::iterator tag=numTagList.begin(); tag!=numTagList.end(); tag++) {
		if(baseAl.HasTag(*tag)) {
			if(*tag=="NM") {
				baseAl.GetTag(*tag,tagVal);
				al.AddTag("ON",string("i"),(int)tagVal);
				al.AddTag(*tag,string("i"),edit);
				continue;
			}
			baseAl.GetTag(*tag,tagVal);
			al.AddTag(*tag,string("i"),(int)tagVal);
		}
	}
	al.AddTag("XC", string("i"), (int)ct);
	for(vector<string>::iterator tag=intTagList.begin(); tag!=intTagList.end(); tag++) {
		if(baseAl.HasTag(*tag)) {
			//char type;
			//baseAl.GetTagType(*tag,type);
			//cout << "Tag: " << *tag << " has type: " << type << endl;
			int32_t tagInt = (int32_t)0;
			baseAl.GetTag(*tag,tagInt);
			al.AddTag(*tag,string("i"),(int)tagInt);
		}
	}
}

ReadFamily::ReadFamily(string key, PairedRead &rp, int lref, int lpos, int rref, int rpos, int order, int insSize, string leftSeq, string rightSeq, int rl)
{
	this->leftReads = false;
	this->rightReads = false;
	this->report = false;
	this->count = 1;
	if(order==1) {
		leftCount=1; rightCount=0;
	} else if(order==2) {
		leftCount=0; rightCount=1;
	}
	this->key = key;
	this->readPair = new PairedRead(rp, order);
	this->lastPair = readPair;
	this->leftRef = lref;
	this->leftPos = lpos;
	this->rightRef = rref;
	this->rightPos = rpos;
	this->insertSize = insSize;
	this->leftRefSeq = leftSeq;
	this->rightRefSeq = rightSeq;
	this->readLen = rl;
	this->leftEdit = 0;
	this->rightEdit = 0;
}

void ReadFamily::clearFamily()
{
	PairedRead *rpit = readPair;
	PairedRead *bkit = readPair;
	while(rpit!=NULL) {
		rpit = rpit->nextPair;
		delete bkit;
		bkit = rpit;
	}
}

void ReadFamily::printFamily()
{
	PairedRead *rpit = readPair;
	cout << "Family with key: " << key << " and contig: " <<  leftPos << ":" << rightPos << " and counts " << count << "=" << leftCount << "+" << rightCount << endl;
	cout << "leftRefSeq:" << endl << leftRefSeq << endl << "rightRefSeq:" << endl << rightRefSeq << endl;
	while(rpit!=NULL) {
		cout << rpit->readID << " " << readLen << endl << leftSeq << endl << rightSeq << endl;
		rpit = rpit->nextPair;
	}
}

string ReadFamily::getKey()
{
	return key;
}

void ReadFamily::getContig(int *ref, int *pos)
{
	*ref = rightRef;
	*pos = rightPos;
	return;
}

bool ReadFamily::checkReadPair(string ID)
{
	PairedRead *rpit = readPair;
	while(rpit!=NULL) {
		if(rpit->readID==ID)
			return true;
		rpit = rpit->nextPair;
	}
	return false;
}

int ReadFamily::missingReadPair()
{
	if(leftCount==0)
		return 1;
	else if(rightCount==0)
		return 2;
	else
		return 0;
}

void ReadFamily::addReadToPair(string ID, string seq, int orient, int order, vector<CigarOp> cigarData, string qual, uint16_t mapq, BamAlignment al)
{
	PairedRead *rpit = readPair;
	while(rpit!=NULL) {
		if(rpit->readID==ID) {
			rpit->addRead(orient, order, seq, cigarData, qual, mapq, leftRefSeq, rightRefSeq, al);
			if(order==1) {
				leftCount++;
			} else if(order==2) {
				rightCount++;
			}
			return;
		}
		rpit = rpit->nextPair;
	}
	cerr << "ERROR: Adding read " << ID << " to existing read pair: Unexpected end of function" << endl;
	exit(1);
}

void ReadFamily::addNewReadPair(PairedRead &RP, int order)
{
	PairedRead *rp = new PairedRead(RP, order);
	lastPair->nextPair = rp;
	lastPair = rp;
	count++;
	if(order==1) {
		leftCount++;
	} else if(order==2) {
		rightCount++;
	}
	return;
}

// REVISIT: If one read end got filtered out, remove from analysis here.
bool ReadFamily::checkBarcodeCollision()
{
	array <int,MAX_READ_LEN> A, C, G, T;
	PairedRead *rpit = readPair;

	A.fill(0);C.fill(0);G.fill(0);T.fill(0);
	while(rpit!=NULL) {
		if(rpit->leftPair) {
			leftReads = true;
			for(int sit=0; sit<readLen; sit++) {
				if(rpit->firstBaseQual[sit]-33 >= MinBaseQual) {
					A[sit] += rpit->firstSeq[0][sit];
					C[sit] += rpit->firstSeq[1][sit];
					G[sit] += rpit->firstSeq[2][sit];
					T[sit] += rpit->firstSeq[3][sit];
				}
			}
		}
		rpit = rpit->nextPair;
	}
	int thresh, ctSum = 0;
	if(leftReads) {
		for(int sit=0; sit<readLen; sit++) {
			int sitePos;

			if(leftPos<=rightPos) {
				sitePos = leftPos+sit;
			} else {
				sitePos = leftPos -(readLen-1)+sit;
			}

			ctSum = A[sit]+C[sit]+G[sit]+T[sit];
			thresh = ceil(0.9*ctSum);
			// If either of the counts are not = total count, check if present in dbSNP too.
			Site site(leftRef,sitePos);
			if(!(ctSum==A[sit]||ctSum==C[sit]||ctSum==G[sit]||ctSum==T[sit])) {
				if(germlineSet.find(site)!=germlineSet.end())
					return false;
			}
		}
	}

	A.fill(0);C.fill(0);G.fill(0);T.fill(0);
	rpit = readPair;
	while(rpit!=NULL) {
		if(rpit->rightPair) {
			rightReads = true;
			for(int sit=0; sit<readLen; sit++) {
				if(rpit->secondBaseQual[sit]-33 >= MinBaseQual) {
					A[sit] += rpit->secondSeq[0][sit];
					C[sit] += rpit->secondSeq[1][sit];
					G[sit] += rpit->secondSeq[2][sit];
					T[sit] += rpit->secondSeq[3][sit];
				}
			}
		}
		rpit = rpit->nextPair;
	}
	if(rightReads) {
		for(int sit=0; sit<readLen; sit++) {
			int sitePos;

			if(leftPos<=rightPos)
				sitePos = rightPos-(readLen-1)+sit;
			else
				sitePos = rightPos+sit;

			ctSum = A[sit]+C[sit]+G[sit]+T[sit];
			thresh = ceil(0.9*ctSum);

			Site site(rightRef,sitePos);
			if(!(ctSum==A[sit]||ctSum==C[sit]||ctSum==G[sit]||ctSum==T[sit])) {
				if(germlineSet.find(site)!=germlineSet.end())
					return false;
			}
		}
	}

	return true;
}

void ReadFamily::callConsensus(int startPosition, int endPosition)
{
	array <int,MAX_READ_LEN> A, C, G, T;
	array <int,MAX_READ_LEN> E, Es;
	PairedRead *rpit = readPair;

	A.fill(0);C.fill(0);G.fill(0);T.fill(0);E.fill(0);
	while(rpit!=NULL) {
		if(rpit->leftPair) {
			leftReads = true;
#ifdef DEBUG
			char charName[100];
			strcpy(charName,rpit->readID.c_str());
			cout << "Debugging:" << charName << endl;
#endif
			for(int sit=0; sit<readLen; sit++) {
				if(rpit->firstSeq[0][sit]+rpit->firstSeq[1][sit]+rpit->firstSeq[2][sit]+rpit->firstSeq[3][sit] == 0)
					continue;
				if(rpit->firstBaseQual[sit]-33 >= MinBaseQual) {
					A[sit] += rpit->firstSeq[0][sit];
					C[sit] += rpit->firstSeq[1][sit];
					G[sit] += rpit->firstSeq[2][sit];
					T[sit] += rpit->firstSeq[3][sit];
					E[sit] += (rpit->firstBaseQual[sit]-33);
				}
			}
		}
		rpit = rpit->nextPair;
	}
	int thresh, ctSum = 0;
	if(leftReads) {
		for(int sit=0; sit<readLen; sit++) {
			int pileupPos;

			if(leftPos<=rightPos) {
				pileupPos = leftPos+sit-startPosition+1;
			} else {
				pileupPos = leftPos-(readLen-1)+sit-startPosition+1;
			}

			ctSum = A[sit]+C[sit]+G[sit]+T[sit];
			thresh = ceil(0.9*ctSum);

			//if((ctSum<2 && leftCount>1) || ctSum==0) {
			if(ctSum==0) {
				leftSeq += 'N';
				//leftSeq += leftRefSeq[sit];
			} else if(C[sit]>=thresh) {
				leftSeq += 'C';
			} else if(A[sit]>=thresh) {
				leftSeq += 'A';
			} else if(G[sit]>=thresh) {
				leftSeq += 'G';
			} else if(T[sit]>=thresh) {
				leftSeq += 'T';
			} else {
				leftSeq += leftRefSeq[sit];
			}

			if(CallMutations==true) addPileupCounts(pileupPos, *leftSeq.rbegin(), leftRefSeq[sit], E[sit], leftPos, rightPos, A[sit], C[sit], G[sit], T[sit], (int)1);
		}
	}
	A.fill(0);C.fill(0);G.fill(0);T.fill(0), Es.fill(0.0);
	rpit = readPair;
	while(rpit!=NULL) {
		if(rpit->rightPair) {
			rightReads = true;
#ifdef DEBUG
		char charName[100];
		strcpy(charName,rpit->readID.c_str());
		cout << "Debugging:" << charName << endl;
#endif
			for(int sit=0; sit<readLen; sit++) {
				if(rpit->secondSeq[0][sit]+rpit->secondSeq[1][sit]+rpit->secondSeq[2][sit]+rpit->secondSeq[3][sit] == 0)
					continue;
				if(rpit->secondBaseQual[sit]-33 >= MinBaseQual) {
					A[sit] += rpit->secondSeq[0][sit];
					C[sit] += rpit->secondSeq[1][sit];
					G[sit] += rpit->secondSeq[2][sit];
					T[sit] += rpit->secondSeq[3][sit];
					Es[sit] += (rpit->secondBaseQual[sit]-33);
				}
			}
		}
		rpit = rpit->nextPair;
	}
	if(rightReads) {
		for(int sit=0; sit<readLen; sit++) {
			int pileupPos;

			if(leftPos<=rightPos)
				pileupPos = rightPos-(readLen-1)+sit-startPosition+1;
			else
				pileupPos = rightPos+sit-startPosition+1;

			int pileupMin, pileupMax;
			if(leftPos<=rightPos) {
				pileupMin = rightPos-(readLen-1)-startPosition;
				pileupMax = leftPos+(readLen-1)-startPosition;
			} else {
				pileupMin = leftPos-(readLen-1)-startPosition;
				pileupMax = rightPos+(readLen-1)-startPosition;
			}

			ctSum = A[sit]+C[sit]+G[sit]+T[sit];
			thresh = ceil(0.9*ctSum);
			//if((ctSum<2 && rightCount>1) || ctSum==0) {
			if(ctSum==0) {
				rightSeq += 'N';
				//rightSeq += rightRefSeq[sit];
			} else if(C[sit]>=thresh) {
				rightSeq += 'C';
			} else if(A[sit]>=thresh) {
				rightSeq += 'A';
			} else if(G[sit]>=thresh) {
				rightSeq += 'G';
			} else if(T[sit]>=thresh) {
				rightSeq += 'T';
			} else {
				rightSeq += rightRefSeq[sit];
			}

			if(leftReads==true&&leftCount>0 && pileupPos>=pileupMin&&pileupPos<=pileupMax)
				continue;
			if(CallMutations==true) addPileupCounts(pileupPos, *rightSeq.rbegin(), rightRefSeq[sit], Es[sit], leftPos, rightPos, A[sit], C[sit], G[sit], T[sit], (int)2);

		}
	}
	return;
}

void ReadFamily::matchConsensus(ReadFamily family_alt)
{
	ReadFamily &family = *this;

	if(family.leftReads) {
		if(!family_alt.rightReads) {
#ifdef DEBUG
			cout << "Family:" << family_alt.getKey() << " has no consensus right read for family:" << family.getKey() << endl;
#endif
		} else
		for(int i=0; i<family.leftSeq.length(); i++) {
			if(family.leftSeq[i]=='N' || family_alt.rightSeq[i]=='N')
				family.leftMatch += 'N';
			else if(family.leftSeq[i]==family_alt.rightSeq[i])
				family.leftMatch += family.leftSeq[i];
			else {
				family.leftMatch += family.leftRefSeq[i]; 
				if(family.leftSeq[i]!=family.leftRefSeq[i]) {
					char refb = family.leftRefSeq[i];
					char altb = family.leftSeq[i];
				}
			}
		}
	}
	if(family.rightReads) {
		if(!family_alt.leftReads) {
#ifdef DEBUG
			cout << "Family:" << family_alt.getKey() << " has no consensus left read for family:" << family.getKey() << endl;
#endif
		} else
		for(int i=0; i<family.rightSeq.length(); i++) {
			if(family.rightSeq[i]=='N' || family_alt.leftSeq[i]=='N')
				family.rightMatch += 'N';
			else if(family.rightSeq[i]==family_alt.leftSeq[i])
				family.rightMatch += family.rightSeq[i];
			else {
				family.rightMatch += family.rightRefSeq[i]; 
				if(family.rightSeq[i]!=family.rightRefSeq[i]) {
					char refb = family.rightRefSeq[i];
					char altb = family.rightSeq[i];
				}
			}
		}
	}
	return;
}

void ReadFamily::addPileup(int startPosition, int endPosition)
{
	if(leftReads==true&&leftCount>0) {
		for(int sit = 1; sit<leftMatch.length(); sit++) {
			int pileupPos;

			if(leftPos<=rightPos)
				pileupPos = leftPos+sit-startPosition;
			else
				pileupPos = leftPos-(readLen-1)+sit-startPosition;
			if(pileupPos<0)
				continue;

			bool artFlag = false;
			if(leftMatch[sit-1] != leftRefSeq[sit-1] && leftMatch[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>leftMatch.length()-depth) && insertSize<=readLen) {
					pileup[4][pileupPos]+=2;
					artFlag = true;
				} else {
					// NOTE: If considering taking median of position inside read for clustered read check, store position in pileup[4] instead of count
					int clusDepth = 0;
					if(sit<readLen/2)
						clusDepth = sit;
					else if(sit>=readLen/2)
						clusDepth = leftMatch.length()-sit;
					else 
						cerr << "ERROR: Not expecting sit = " << sit << endl;
					pileupCluster[pileupPos].push_back(clusDepth);
				}
			}

			PairedRead *rpit = readPair;

			if(CallMutations==true) addLOD(rpit, pileupPos, leftSeq[sit-1], leftRefSeq[sit-1], sit, (int)1);
			if(artFlag) {
				artFlag = false;
				continue;
			}
			switch(leftMatch[sit-1]) {
				case 'A': pileup[0][pileupPos]+=2;break;
				case 'C': pileup[1][pileupPos]+=2;break;
				case 'G': pileup[2][pileupPos]+=2;break;
				case 'T': pileup[3][pileupPos]+=2;break;
				default: break;
			}
		}
		}
		if(rightReads==true&&rightCount>0) {
			int pileupMin, pileupMax;
			if(leftPos<=rightPos) {
				pileupMin = rightPos-(readLen-1)-startPosition;
				pileupMax = leftPos+(readLen-1)-startPosition;
			} else {
				pileupMin = leftPos-(readLen-1)-startPosition;
				pileupMax = rightPos+(readLen-1)-startPosition;
			}
		for(int sit = 1; sit<rightMatch.length(); sit++) {
			int pileupPos;
			if(leftPos<=rightPos)
				pileupPos = rightPos-(readLen-1)+sit-startPosition;
			else
				pileupPos = rightPos+sit-startPosition;
			if(pileupPos<0)
				continue;

			bool artFlag = false;
			if(rightMatch[sit-1] != rightRefSeq[sit-1] && rightMatch[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>rightMatch.length()-depth) && insertSize<=readLen) {
					pileup[4][pileupPos]+=2;
					artFlag = true;
				} else {
					int clusDepth = 0;
					if(sit<readLen/2)
						clusDepth = sit;
					else if(sit>=readLen/2)
						clusDepth = rightMatch.length()-sit;
					else
						cerr << "ERROR: Not expecting sit = " << sit << endl;
					pileupCluster[pileupPos].push_back(clusDepth);
				}
			}

			if(leftReads==true&&leftCount>0 && pileupPos>=pileupMin&&pileupPos<=pileupMax)
				continue;

			PairedRead *rpit = readPair;
			if(CallMutations==true) addLOD(rpit, pileupPos, rightSeq[sit-1], rightRefSeq[sit-1], sit, (int)2);
			if(artFlag) {
				artFlag = false;
				continue;
			}
			switch(rightMatch[sit-1]) {
				case 'A': pileup[0][pileupPos]+=2;break;
				case 'C': pileup[1][pileupPos]+=2;break;
				case 'G': pileup[2][pileupPos]+=2;break;
				case 'T': pileup[3][pileupPos]+=2;break;
				default: break;
			}
		}
	}
	return;
}

void ReadFamily::addSingleStrandedPileup(int startPosition, int endPosition)
{
	if(count<0)
		return;
	if(leftReads==true&&leftCount>1) {
		for(int sit = 1; sit<=leftSeq.length(); sit++) {
			if(leftSeq[sit-1] != leftRefSeq[sit-1])
				leftEdit++;

			int pileupPos;
			if(leftPos<=rightPos)
				pileupPos = leftPos+sit-startPosition;
			else
				pileupPos = leftPos-(readLen-1)+sit-startPosition;
			if(pileupPos<0)
				continue;

			bool artFlag = false;
			if(leftSeq[sit-1] != leftRefSeq[sit-1] && leftSeq[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>leftSeq.length()-depth) && insertSize<=readLen) {
					pileupSS[4][pileupPos]++;
					artFlag = true;
				}
					int clusDepth = 0;
					if(sit<=readLen/2)
						clusDepth = sit;
					else if(sit>readLen/2)
						clusDepth = leftSeq.length()-sit+1;
					else 
						cerr << "ERROR: Not expecting sit = " << sit << endl;
					pileupSSCluster[pileupPos].push_back(clusDepth);
			}
			if(artFlag) {
				artFlag = false;
				pileupSS[4][pileupPos]++;
				continue;
			}
			switch(leftSeq[sit-1]) {
				case 'A': pileupSS[0][pileupPos]++;break;
				case 'C': pileupSS[1][pileupPos]++;break;
				case 'G': pileupSS[2][pileupPos]++;break;
				case 'T': pileupSS[3][pileupPos]++;break;
				default: break;
			}
			//readPileupArray[sit] += leftCount;
		}
	}
	if(rightReads==true&&rightCount>1) {
		int pileupMin, pileupMax;
		if(leftPos<=rightPos) {
			pileupMin = rightPos-(readLen-1)-startPosition;
			pileupMax = leftPos+(readLen-1)-startPosition;
		} else {
			pileupMin = leftPos-(readLen-1)-startPosition;
			pileupMax = rightPos+(readLen-1)-startPosition;
		}
		for(int sit = 1; sit<=rightSeq.length(); sit++) {
			if(rightSeq[sit-1] != rightRefSeq[sit-1])
				rightEdit++;

			int pileupPos;
			if(leftPos<=rightPos)
				pileupPos = rightPos-(readLen-1)+sit-startPosition;
			else
				pileupPos = rightPos+sit-startPosition;
			if(pileupPos<0)
				continue;
	
			bool artFlag = false;
			if(rightSeq[sit-1] != rightRefSeq[sit-1] && rightSeq[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>rightSeq.length()-depth) && insertSize<=readLen) {
					pileupSS[4][pileupPos]++;
					artFlag = true;
				}
				int clusDepth = 0;
				if(sit<=readLen/2)
					clusDepth = sit;
				else if(sit>readLen/2)
					clusDepth = rightSeq.length()-sit+1;
				else
					cerr << "ERROR: Not expecting sit = " << sit << endl;
				pileupSSCluster[pileupPos].push_back(clusDepth);
			}
			if(artFlag) {
				artFlag = false;
				pileupSS[4][pileupPos]++;
				continue;
			}
			if(leftReads==true&&leftCount>0 && pileupPos>=pileupMin&&pileupPos<=pileupMax)
				continue;
	
			switch(rightSeq[sit-1]) {
				case 'A': pileupSS[0][pileupPos]++;break;
				case 'C': pileupSS[1][pileupPos]++;break;
				case 'G': pileupSS[2][pileupPos]++;break;
				case 'T': pileupSS[3][pileupPos]++;break;
				default: break;
			}
			// readPileupArray[sit] += rightCount;
		}
	}
	return;
}

// Family size expected to be == 1 and 1 only
void ReadFamily::addSingletonPileup(int startPosition, int endPosition)
{
	if(count<0)
		return;
	if(leftReads==true&&leftCount==1) {
		for(int sit = 1; sit<=leftSeq.length(); sit++) {
			if(leftSeq[sit-1] != leftRefSeq[sit-1])
				leftEdit++;

			int pileupPos;
			if(leftPos<=rightPos)
				pileupPos = leftPos+sit-startPosition;
			else
				pileupPos = leftPos-(readLen-1)+sit-startPosition;
			if(pileupPos<0)
				continue;

			bool artFlag = false;
			if(leftSeq[sit-1] != leftRefSeq[sit-1] && leftSeq[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>leftSeq.length()-depth) && insertSize<=readLen) {
					pileupSingleton[4][pileupPos]++;
					artFlag = true;
				}
				int clusDepth = 0;
				if(sit<=readLen/2)
					clusDepth = sit;
				else if(sit>readLen/2)
					clusDepth = leftSeq.length()-sit+1;
				else 
					cerr << "ERROR: Not expecting sit = " << sit << endl;
				pileupSingletonCluster[pileupPos].push_back(clusDepth);
			}
			if(artFlag) {
				artFlag = false;
				pileupSingleton[4][pileupPos]++;
				continue;
			}
			switch(leftSeq[sit-1]) {
				case 'A': pileupSingleton[0][pileupPos]++;break;
				case 'C': pileupSingleton[1][pileupPos]++;break;
				case 'G': pileupSingleton[2][pileupPos]++;break;
				case 'T': pileupSingleton[3][pileupPos]++;break;
				default: break;
			}
			// readPileupArray[sit] += leftCount;
		}
	}
	if(rightReads==true&&rightCount==1) {
		int pileupMin, pileupMax;
		if(leftPos<=rightPos) {
			pileupMin = rightPos-(readLen-1)-startPosition;
			pileupMax = leftPos+(readLen-1)-startPosition;
		} else {
			pileupMin = leftPos-(readLen-1)-startPosition;
			pileupMax = rightPos+(readLen-1)-startPosition;
		}
		for(int sit = 1; sit<=rightSeq.length(); sit++) {
			if(rightSeq[sit-1] != rightRefSeq[sit-1])
				rightEdit++;

			int pileupPos;
			if(leftPos<=rightPos)
				pileupPos = rightPos-(readLen-1)+sit-startPosition;
			else
				pileupPos = rightPos+sit-startPosition;
			if(pileupPos<0)
				continue;

			bool artFlag = false;
			if(rightSeq[sit-1] != rightRefSeq[sit-1] && rightSeq[sit-1] != 'N') {
				int depth = readLen - insertSize;
				if((sit<=depth||sit>rightSeq.length()-depth) && insertSize<=readLen) {
					pileupSingleton[4][pileupPos]++;
					artFlag = true;
				}
				int clusDepth = 0;
				if(sit<=readLen/2)
					clusDepth = sit;
				else if(sit>readLen/2)
					clusDepth = rightSeq.length()-sit+1;
				else
					cerr << "ERROR: Not expecting sit = " << sit << endl;
				pileupSingletonCluster[pileupPos].push_back(clusDepth);
			}
			if(artFlag) {
				artFlag = false;
				pileupSingleton[4][pileupPos]++;
				continue;
			}
			if(leftReads==true&&leftCount>0 && pileupPos>=pileupMin&&pileupPos<=pileupMax)
				continue;

			switch(rightSeq[sit-1]) {
				case 'A': pileupSingleton[0][pileupPos]++;break;
				case 'C': pileupSingleton[1][pileupPos]++;break;
				case 'G': pileupSingleton[2][pileupPos]++;break;
				case 'T': pileupSingleton[3][pileupPos]++;break;
				default: break;
			}
			// readPileupArray[sit] += rightCount;
		}
	}
	return;
}

/*
int ReadFamily::dumpNonConsensusBAM()
{
	string famKey = this->key;
	PairedRead *rpit = readPair;
	int firstUnique = 0, secondUnique=0;
	int leftOverlap = false, rightOverlap = false;
	int refID, readStart, readEnd;

	while(rpit!=NULL) {
		bool markUnique = true;
		if(rpit->leftPair&&!leftOverlap) {
			refID = rpit->firstAlignment.RefID;
			readStart = rpit->firstAlignment.Position;
			readEnd = readStart+rpit->firstAlignment.Length-1;
			if(priorOverlap(refID, readStart, readEnd)) {
				cout << "Skipping first segment of read pair: " << rpit->firstAlignment.Name << " " << refID << ":" << readStart << " already contained in interval" << endl;
				leftOverlap = true;
			} else if(firstUnique<1&&markUnique==true) {
				rpit->firstAlignment.SetIsDuplicate(false);
				firstUnique++;
			} else {
				rpit->firstAlignment.SetIsDuplicate(true);
			}
			if(!leftOverlap) {
				nonConsensusAlignmentVector.push_back(rpit->firstAlignment);
			}
		}
		if(rpit->rightPair&&!rightOverlap) {
			refID = rpit->secondAlignment.RefID;
			readStart = rpit->secondAlignment.Position;
			readEnd = readStart+rpit->secondAlignment.Length-1;
			if(priorOverlap(refID, readStart, readEnd)) {
				cout << "Skipping second segment of read pair: " << rpit->secondAlignment.Name << " " << refID << ":" << readStart << " already contained in interval" << endl;
				rightOverlap = true;
			} else if(secondUnique<1&&markUnique==true) {
				rpit->secondAlignment.SetIsDuplicate(false);
				secondUnique++;
			} else {
				rpit->secondAlignment.SetIsDuplicate(true);
			}
			if(!rightOverlap) {
				nonConsensusAlignmentVector.push_back(rpit->secondAlignment);
			}
		}
	rpit = rpit->nextPair;
	}
}
*/

int ReadFamily::dumpConsensusBAM()
{
	PairedRead *rpit = readPair;
	int firstUnique = 0, secondUnique=0;
	int leftOverlap = false, rightOverlap = false;
	int refID, readStart, readEnd;
	vector<CigarOp> cigarVec;
	vector<string> keyVec = split(getKey(),'-');
	string barcode = keyVec[2] + '-' + keyVec[3];

	while(rpit!=NULL) {
		bool markUnique = true;
		if(UseSingletonReads || (leftReads==true&&leftCount>1)) {
/*
		if(rpit->leftPair && rpit->rightPair)
			markUnique=true;
		else if(this->leftCount==0 || this->rightCount==0)
			markUnique=true;
		else
			markUnique=false;
*/
		if(rpit->leftPair&&!leftOverlap) {
			refID = rpit->firstAlignment.RefID;
			readStart = rpit->firstAlignment.Position;
			readEnd = readStart+rpit->firstAlignment.Length-1;
			BamAlignment al;
			if(priorOverlap(refID, readStart, readEnd, readLen)) {
#ifdef DEBUG
				cout << "Skipping first segment of read pair: " << rpit->firstAlignment.Name << " " << refID << ":" << readStart << " already contained in interval" << endl;
#endif
				leftOverlap = true;
			} else if(firstUnique<ConsensusUniqueCount && markUnique==true) { // Duplicate?
				al.SetIsDuplicate(false);
				firstUnique++;
			} else {
				al.SetIsDuplicate(true);
			}
			if(!leftOverlap) {
				al.SetIsFirstMate(true); // First segment?
				al.SetIsProperPair(true); // All segments properly aligned?
				al.RefID = leftRef;
				if(rpit->firstOrient==2) {
					al.SetIsReverseStrand(true);
					al.Position = leftPos - leftSeq.length() + 1;
					al.InsertSize = -insertSize;
				} else {
					al.SetIsReverseStrand(false);
					al.Position = leftPos;
					al.InsertSize = insertSize;
				}
				if(rpit->rightPair) {
					al.SetIsPaired(true);
					al.MateRefID = rightRef;
					if(rpit->secondOrient==2) {
						al.SetIsMateReverseStrand(true);
						al.MatePosition = rightPos - rightSeq.length() + 1;
					} else if(rpit->secondOrient==1) {
						al.SetIsMateReverseStrand(false);
						al.MatePosition = rightPos;
					}
				} else {
					al.SetIsPaired(false);
					al.MateRefID = leftRef;
					if(rpit->firstOrient==1) {
						al.SetIsMateReverseStrand(true);
						al.MatePosition = rightPos - leftSeq.length() + 1;
					} else if(rpit->firstOrient==2) {
						al.SetIsMateReverseStrand(false);
						al.MatePosition = rightPos;
					}
				}
				al.MapQuality = rpit->firstMapq;
				int rlen = (rpit->firstBaseQual).length();
				al.Qualities = rpit->firstBaseQual;
				al.QueryBases = leftSeq.substr(0,rlen);
				CigarOp cigar((const char)'M', (uint32_t)rlen);
				cigarVec.clear();
				cigarVec.push_back(cigar);
				al.CigarData = cigarVec;
				al.Name = rpit->readID;
				al.Length = rlen;
				rpit->dumpTags(al, keyVec[3]+'-'+keyVec[2], leftEdit, (int)1, getKey(), leftCount); //why is this always throwing 0?
				consensusAlignmentVector.push_back(al);
			}
		}
		}
		if(UseSingletonReads || (rightReads==true&&rightCount>1)) {
		if(rpit->rightPair&&!rightOverlap) {
			refID = rpit->secondAlignment.RefID;
			readStart = rpit->secondAlignment.Position;
			readEnd = readStart+rpit->secondAlignment.Length-1;
			BamAlignment al;
			if(priorOverlap(refID, readStart, readEnd, readLen)) {
#ifdef DEBUG
				cout << "Skipping second segment of read pair: " << rpit->secondAlignment.Name << " " << refID << ":" << readStart << " already contained in interval" << endl;
#endif
				rightOverlap = true;
			} else if(secondUnique<ConsensusUniqueCount && markUnique==true) {
				al.SetIsDuplicate(false);
				secondUnique++;
			} else {
				al.SetIsDuplicate(true);
			}
			if(!rightOverlap) {
				al.SetIsSecondMate(true);
				al.SetIsProperPair(true);
				al.RefID = rightRef;
				if(rpit->secondOrient==2) {
					al.SetIsReverseStrand(true);
					al.Position = rightPos - rightSeq.length() + 1;
					al.InsertSize = -insertSize;
				} else {
					al.SetIsReverseStrand(false);
					al.Position = rightPos;
					al.InsertSize = insertSize;
				}
				if(rpit->leftPair) {
					al.SetIsPaired(true);
					al.MateRefID = leftRef;
					if(rpit->firstOrient==2) {
						al.SetIsMateReverseStrand(true);
						al.MatePosition = leftPos - leftSeq.length() + 1;
					} else if(rpit->firstOrient==1) {
						al.SetIsMateReverseStrand(false);
						al.MatePosition = leftPos;
					}
				} else {
					al.SetIsPaired(false);
					al.MateRefID = rightRef;
					if(rpit->secondOrient==1) {
						al.SetIsMateReverseStrand(true);
						al.MatePosition = leftPos - rightSeq.length() + 1;
					} else if(rpit->secondOrient==2) {
						al.SetIsMateReverseStrand(false);
						al.MatePosition = leftPos;
					}
				}
				al.MapQuality = rpit->secondMapq;
				int rlen = (rpit->secondBaseQual).length();
				al.Qualities = rpit->secondBaseQual;
				al.QueryBases = rightSeq.substr(0,rlen);
				CigarOp cigar((const char)'M', (uint32_t)rlen);
				cigarVec.clear();
				cigarVec.push_back(cigar);
				al.CigarData = cigarVec;
				al.Name = rpit->readID;
				al.Length = rlen;
				rpit->dumpTags(al, keyVec[3]+'-'+keyVec[2], rightEdit, (int)2, getKey(), rightCount);
				consensusAlignmentVector.push_back(al);
			}
		}
		}
		rpit = rpit->nextPair;
	}
}

