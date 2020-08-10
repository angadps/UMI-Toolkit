#include "readUtils.h"

using namespace std;
using namespace BamTools;

int MaxUMIEditDist = 2;
RefOrder refOrder;

//vector<string> P5 = {"TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA", "AGGCGAAG", "TAATCTTA", "CAGGACGT", "GTACTGAC"};
//vector<string> P7 = {"ATTACTCG", "TCCGGAGA", "CGCTCATT", "GAGATTCC", "ATTCAGAA", "GAATTCGT", "CTGAAGCT", "TAATGCGC"};
vector<string> P5 = {"ATCACGCTTAT", "CGATGTCTTAT", "TTAGGCTTAT", "TGACCATTAT", "ACAGTGTTAT", "GCCAATTTAT", "CTTGTACGTTAT", "GTGAAACGTTAT"};
vector<string> P7 = {"ATCACGCTTAT", "CGATGTCTTAT", "TTAGGCTTAT", "TGACCATTAT", "ACAGTGTTAT", "GCCAATTTAT", "CTTGTACGTTAT", "GTGAAACGTTAT"};

int bcMismatch(string barcode, string barcodeSet, string &barcodeAligned, string &barcodeSetAligned)
{
	int diffAligned = 0;
	int len = barcode.size();
	int slen = barcodeSet.size();

	int penalty = align(barcode, barcodeSet, barcodeAligned, barcodeSetAligned);
	len = barcodeAligned.size(), slen = barcodeSetAligned.size();
	for(int i=0; i<len&&i<slen; i++) {
		diffAligned += (barcodeAligned[i]!=barcodeSetAligned[i]);
	}
	return diffAligned;
}

bool generate_alt_barcode(string &barcode, string &barcode_alt, bool fix)
{
	string barcodeAligned, barcodeSetAligned;
	vector<string> p = split(barcode,'-');

	if(p.size()<2)
		return false;
	if(fix)
		barcode.clear();
	string p5 = p[0];
	string p7 = p[1];
	bool match = false;
	int ct = 0, p5Ct = 0, p7Ct = 0;

	for(int i=0; i<P7.size(); i++) {
		if(P7[i].compare(p7)==0) {
			match = true;
			p7Ct = 0;
			if(fix)
				barcode += P7[i];
			barcode_alt += P5[i];
			break;
		}
	}

	if(match==false) {
		for(int i=0; i<P7.size(); i++) {
			ct = bcMismatch(p7,P7[i],barcodeAligned,barcodeSetAligned);
			if(ct<=MaxUMIEditDist) {
				p7Ct = ct;
				if(fix)
					barcode += P7[i];
				barcode_alt += P5[i];
#ifdef DEBUG
				if(ct>0) {
					cout << "Before alignment: " << p7 << " " << P7[i] << endl;
					cout << "After alignment: " << barcodeAligned << " " << barcodeSetAligned << endl;
				}
#endif
				break;
			}
		}
	}

	match = false;
	if(fix)
		barcode += "-";
	barcode_alt += "-";

	for(int i=0; i<P5.size(); i++) {
		if(P5[i].compare(p5)==0) {
			match = true;
			p5Ct = 0;
			if(fix)
				barcode += P5[i];
			barcode_alt += P5[i];
			break;
		}
	}

	if(match==false) {
		for(int i=0; i<P5.size(); i++) {
			ct = bcMismatch(p5,P5[i],barcodeAligned,barcodeSetAligned);
			if(ct<=MaxUMIEditDist) {
				p5Ct = ct;
				if(fix)
					barcode += P5[i];
				barcode_alt += P7[i];
#ifdef DEBUG
				if(ct>0) {
					cout << "Before alignment: " << p5 << " " << P5[i] << endl;
					cout << "After alignment: " << barcodeAligned << " " << barcodeSetAligned << endl;
				}
#endif
				break;
			}
		}
	}

	if(barcode_alt.length()<20 || p5Ct+p7Ct>MaxUMIEditDist*2) {
#ifdef DEBUG
		cerr << "WARNING: Could not generate alt barcode for " << barcode << ". Returning " << barcode_alt << endl;
#endif
		return false;
	}
	return true;
}

int getStartSoftClipLength(BamAlignment al)
{
	int length = 0;
	vector<CigarOp>::iterator it=al.CigarData.begin();
	if((*it).Type=='S')
		length = (*it).Length;
	return length;
}

int getEndSoftClipLength(BamAlignment al)
{
	int length = 0;
	vector<CigarOp>::reverse_iterator it=al.CigarData.rbegin();
	if((*it).Type=='S')
		length = (*it).Length;
	return length;
}

// Looks for soft clip at start of read alignment only
int getClippedStartPosition(BamAlignment al)
{
	int position = al.Position;
	vector<CigarOp>::iterator it=al.CigarData.begin();
	if((*it).Type=='S')
		position = position - (*it).Length;
	return position;
}

int inslen(BamAlignment al)
{
	int len = 0;
	for(vector<CigarOp>::iterator it=al.CigarData.begin(); it!=al.CigarData.end(); it++) {
		if((*it).Type=='I')
			len = len + (*it).Length;
	}
	return len;
}

int indelLen(BamAlignment al)
{
	int len = 0;
	for(vector<CigarOp>::iterator it=al.CigarData.begin(); it!=al.CigarData.end(); it++) {
		if((*it).Type=='I' || (*it).Type=='D')
			len = len + (*it).Length;
	}
	return len;
}

int chrRank(string chr)
{
	if(refOrder.find(chr)==refOrder.end()) {
		cerr << "Invalid chromosome name: " << chr << " . Exiting..." << endl;
		exit(1);
	}
	return refOrder[chr];
}

