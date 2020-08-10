#include "utils.h"

const size_t alphabets = 26;

using namespace std;
using namespace BamTools;

vector<string> &split(const string &s, char delim, vector<string> &elems)
{
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

vector<string> split(const string &s, char delim)
{
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}

double medianOfVector(vector<int> intVector, int vecSize)
{
	int middle;
	double median;
	 
	sort(intVector.begin(),intVector.begin()+vecSize);
	middle = vecSize/2;
	 
	if(vecSize%2) {
		median = intVector[middle]/1.0;
	} else {
		median = (intVector[middle] + intVector[middle-1])/2.0;
	}
	return median;
}

/* 
 *  * Returns the Needleman-Wunsch score for the best alignment of a and b
 *  * and stores the aligned sequences in a_aligned and b_aligned
 */
int align(const string &a, const string &b, string &a_aligned, string &b_aligned)
{
	size_t n = a.size();
	size_t m = b.size();
	int gap_penalty = 2;
/* 
 * alpha[i][j] = penalty for matching the ith alphabet with the jth alphabet.
 * Here: Penalty for matching an alphabet with anoter one is 1 
 * Penalty for matching an alphabet with itself is 0
 */
	int alpha[alphabets][alphabets];
	for (size_t i = 0; i < alphabets; ++i) {
		for (size_t j = 0; j < alphabets; ++j) {
			if (i == j) 
				alpha[i][j] = 0;
			else
				alpha[i][j] = 1;
		}
	}
	vector<vector<int> > A(n + 1, vector<int>(m + 1));

	for (size_t i = 0; i <= m; ++i)
		A[0][i] = gap_penalty * i;
	for (size_t i = 0; i <= n; ++i)
		A[i][0] = gap_penalty * i;

	for (size_t i = 1; i <= n; ++i) {
		for (size_t j = 1; j <= m; ++j) {
			char x_i = a[i-1];
			char y_j = b[j-1];
			A[i][j] = min(min(A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A'], A[i-1][j] + gap_penalty),A[i][j-1] + gap_penalty);
		}
	}

	a_aligned = "";
	b_aligned = "";
	size_t j = m;
	size_t i = n;
	//for (; i >= 1 && j >= 1; --i)
	for (; i >= 1 && j >= 1;) {
		char x_i = a[i-1];
		char y_j = b[j-1];
		if (A[i][j] == A[i-1][j-1] + alpha[x_i - 'A'][y_j - 'A']) {
/*
 * I think prepending chars this way to a string is very inefficient.
 * Is there any better way of doing this without using C-style strings?
 */
			a_aligned = x_i + a_aligned;
			b_aligned = y_j + b_aligned;
			--j;
			--i;
		} else if (A[i][j] == A[i-1][j] + gap_penalty) {
			a_aligned = x_i + a_aligned;
			b_aligned = '-' + b_aligned;
			--i;
		} else {
			a_aligned = '-' + a_aligned;
			b_aligned = y_j + b_aligned;
			--j;
		}
	}

	while (i >= 1 && j < 1) {
		a_aligned = a[i-1] + a_aligned;
		b_aligned = '-' + b_aligned;
		--i;
	}
	while (j >= 1 && i < 1) {
		a_aligned = '-' + a_aligned;
		b_aligned = b[j-1] + b_aligned;
		--j;
	}
	return A[n][m];
}

bool verifyCoordinates(string refFasta, string inputBamFile)
{
	int i = 0, maxChr = 5000;
	string refChr[maxChr], bamChr[maxChr];
	int refPos[maxChr], bamPos[maxChr];

	char delim = '\t';
	string refLine;
	vector<string> refCoords;
	ifstream fai((refFasta+".fai").c_str());
	if(!fai.is_open()) {
		cerr << "ERROR: Failed to open fasta.fai file: " << refFasta << endl;
		return false;
	}
	while(getline(fai, refLine)) {
		refCoords = split(refLine, delim);
		refChr[i] = refCoords[0];
		refPos[i++] = stoi(refCoords[1]);
	}
	fai.close();

	BamReader reader;
	if(!reader.Open(inputBamFile)) {
		cerr << "ERROR: Could not open BAM file:" << inputBamFile << endl;
		return false;
	}
	i = 0;
	RefVector referenceSequences = reader.GetReferenceData();
	for(vector<RefData>::iterator it=referenceSequences.begin(); it!=referenceSequences.end(); it++) {
		bamChr[i] = (*it).RefName;
		bamPos[i++] = (*it).RefLength;
	}
	reader.Close();

	for(--i; i>=0; i--) {
		if(refOrder.find(refChr[i])==refOrder.end()) {
			refOrder.insert(std::pair<string,int> (refChr[i],i));
		} else {
			cerr << "ERROR: Repeated chromosome: " << refChr[i] << " . Exiting..." << endl;
			return false;
		}
		if(refChr[i]!=bamChr[i] || refPos[i]!=bamPos[i]) {
			return false;
		}
	}
	return true;
}

/*
bool verifyCoordinates(string refFasta, string inputBamFile, string intervalList)
{
	int i = 0, maxChr = 5000;
	string refChr[maxChr], bamChr[maxChr], intChr[maxChr];
	int refPos[maxChr], bamPos[maxChr], intPos[maxChr];

	char delim = '\t';
	string refLine;
	vector<string> refCoords;
	ifstream fai((refFasta+".fai").c_str());
	if(!fai.is_open()) {
		cerr << "ERROR: Failed to open fasta.fai file: " << refFasta << endl;
		return false;
	}
	while(getline(fai, refLine)) {
		refCoords = split(refLine, delim);
		refChr[i] = refCoords[0];
		refPos[i++] = stoi(refCoords[1]);
	}
	fai.close();

	BamReader reader;
	if(!reader.Open(inputBamFile)) {
		cerr << "ERROR: Could not open BAM file:" << inputBamFile << endl;
		return false;
	}
	i = 0;
	RefVector referenceSequences = reader.GetReferenceData();
	for(vector<RefData>::iterator it=referenceSequences.begin(); it!=referenceSequences.end(); it++) {
		bamChr[i] = (*it).RefName;
		bamPos[i++] = (*it).RefLength;
	}
	reader.Close();

	string intervalLine;
	vector<string> intervalText, intervalRef, intervalPos;
	ifstream intervalFile(intervalList);
	if(!intervalFile.is_open()) {
		cerr << "ERROR: Could not open intervals list: " << intervalList << endl;
		return false;
	}
	while(getline(intervalFile, intervalLine)) {
		if(intervalLine.find("@SQ")!=string::npos)
			break;
	}
	i = 0;
	intervalText = split(intervalLine,delim);
	intervalRef = split(intervalText[1],':');
	intervalPos = split(intervalText[2],':');
	intChr[i] = intervalRef[1];
	intPos[i++] = stoi(intervalPos[1]);
	while(getline(intervalFile, intervalLine)) {
		if(intervalLine.find("@SQ")==string::npos)
			break;
		else {
			intervalText = split(intervalLine,delim);
			intervalRef = split(intervalText[1],':');
			intervalPos = split(intervalText[2],':');
			intChr[i] = intervalRef[1];
			intPos[i++] = stoi(intervalPos[1]);
		}
	}
	intervalFile.close();

	for(--i; i>=0; i--) {
		if(refOrder.find(refChr[i])==refOrder.end()) {
			refOrder.insert(std::pair<string,int> (refChr[i],i));
		} else {
			cerr << "ERROR: Repeated chromosome: " << refChr[i] << " . Exiting..." << endl;
			return false;
		}
		if(refChr[i]!=bamChr[i]||refChr[i]!=intChr[i] || refPos[i]!=bamPos[i]||refPos[i]!=intPos[i]) {
			return false;
		}
	}
	return true;
}
*/

bool getIntCoords(ifstream& intervalFile, string& intRef, int& intStart, int& intEnd, int padding_length)
{
	char delim = '\t';
	string intervalLine;
	vector<string> intervalRegion;

	if(getline(intervalFile, intervalLine)) {
		if(intervalLine.find("@")==string::npos) {
			intervalRegion = split(intervalLine, delim);
			if(intervalRegion.size() < 3) {
				cerr << "Invalid interval: " << intervalLine << "Exiting..." << endl;
				exit(1);
			}
			intRef = intervalRegion[0];
			intStart = stoi(intervalRegion[1]) - padding_length;
			intEnd = stoi(intervalRegion[2]) + padding_length;
			return true;
		}
	}
	return false;
}

string getPaddedIntervalList(string intervalList, string outputPrefix, int padding_length)
{
	char delim = '\t';
	string intervalLine, intRef1, intRef2;
	int intStart1, intEnd1, intStart2, intEnd2;
	vector<string> intervalRegion;

	string paddedIntervalList = outputPrefix + ".tmp.interval_list";
	ifstream intervalFile(intervalList);
	ofstream paddedIntervalFile(paddedIntervalList, ofstream::out);

	if(!intervalFile.is_open()) {
		cerr << "ERROR: Could not open intervals list: " << intervalList << endl;
		exit(1);
	}

	// Read header and first interval line
	while(getline(intervalFile, intervalLine)) {
		if(intervalLine.find("@")!=string::npos) {
			paddedIntervalFile << intervalLine << endl;
		} else {
			break;
		}			
	}

	if(intervalLine.empty()) {
		cerr << "Interval list file: " << intervalList << " is empty. Exiting..." << endl;
		exit(1);
	}

	// Reading first line in interval list
	intervalRegion = split(intervalLine, delim);
	if(intervalRegion.size() < 3) {
		cerr << "Invalid interval: " << intervalLine << " in interval file: " << intervalList << "Exiting..." << endl;
		exit(1);
	}
	intRef1 = intervalRegion[0];
	intStart1 = stoi(intervalRegion[1]) - padding_length;
	intEnd1 = stoi(intervalRegion[2]) + padding_length;

	while(getIntCoords(intervalFile, intRef2, intStart2, intEnd2, padding_length)) {
		if(intRef1 == intRef2) {
			if(intStart2 - intEnd1 > 1) { // Not adjacent or overlapping
				paddedIntervalFile << intRef1 << "\t" << intStart1 << "\t" << intEnd1 << endl;
				intRef1 = intRef2; intStart1 = intStart2; intEnd1 = intEnd2;
			} else {
				intEnd1 = intEnd2;
			}
		} else { // ref1 != ref2
			paddedIntervalFile << intRef1 << "\t" << intStart1 << "\t" << intEnd1 << endl;
			intRef1 = intRef2; intStart1 = intStart2; intEnd1 = intEnd2;
		}
	}
	paddedIntervalFile << intRef1 << "\t" << intStart1 << "\t" << intEnd1 << endl;

	intervalFile.close();
	paddedIntervalFile.close();
	return paddedIntervalList;
}

// Needs revision in case planned to use again
/*
bool readGermlineFile(string germlineList)
{
	int germlineChr, germlinePos;
	ifstream germlineFile(germlineList);
	string germlineLine;
	bool flag = true;
	char delim = '_';

	if(!germlineFile.is_open()) {
		cerr << "ERROR: Could not open germline file: " << germlineList << endl;
		return false;
	}

	vector<string> germlineCoords;
	getline(germlineFile, germlineLine); // Header
	while(getline(germlineFile,germlineLine)) {
		if(germlineLine[0]!='c' && germlineLine[1]!='h' && germlineLine[2]!='r')
			continue;
		germlineCoords = split(germlineLine, delim);
		germlineChr = chrRank(germlineCoords[0]);
		germlinePos = stoi(germlineCoords[1]);
		Site site(germlineChr-1, germlinePos);
		germlineSet.insert(site);
	}
	germlineFile.close();

	cout << "Germline set size = " << germlineSet.size() << endl;
	return true;
}
*/

