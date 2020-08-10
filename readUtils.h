#ifndef READUTILS_h
#define READUTILS_h

#include <map>
#include "utils.h"
#include "api/BamAlignment.h"

extern int MaxUMIEditDist;
typedef std::map<std::string, int> RefOrder;
extern RefOrder refOrder;

int bcMismatch(std::string barcode, std::string barcodeSet, std::string &barcodeAligned, std::string &barcodeSetAligned);
int indelLen(BamTools::BamAlignment al);
int inslen(BamTools::BamAlignment al);
int getClippedStartPosition(BamTools::BamAlignment al);
int getEndSoftClipLength(BamTools::BamAlignment al);
int getStartSoftClipLength(BamTools::BamAlignment al);
bool generate_alt_barcode(std::string &barcode, std::string &barcode_alt, bool fix);
int chrRank(std::string chr);

#endif

