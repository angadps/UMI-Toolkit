#include <vector>
#include "mutations.h"

using namespace std;

array<array<int,MAX_INTERVAL_LEN>,9> strandArray;
array<int,MAX_INTERVAL_LEN> refPile;
array<int,MAX_INTERVAL_LEN> altPile;
array<double,MAX_INTERVAL_LEN> errPile;
array<double,MAX_INTERVAL_LEN> tumor_lod;

vector<Mutation> mutationVector;
//set<Site> dbSNPSet;
ofstream mutFile, mutVcf;

char numToBase(array<array<int,MAX_READ_LEN>,4> &firstSeq, int sit)
{
	char base = 'N';
	if(firstSeq[0][sit]==1)
		base = 'A';
	else if(firstSeq[1][sit]==1)
		base = 'C';
	else if(firstSeq[2][sit]==1)
		base = 'G';
	else if(firstSeq[3][sit]==1)
		base = 'T';
	else {
		cerr << "ERROR: Array has no base assigned at position:" << sit << endl;
	}
	return base;
}

void addPileupCounts(int pileupPos, char alSeq, char refSeq, double E, int leftPos, int rightPos, int A, int C, int G, int T, int ind)
{
	if(pileupPos<0)
		return;
	if(alSeq==refSeq)
		refPile[pileupPos]++;
	else {
		altPile[pileupPos]++;
	}

	int ctSum = A + C + G + T;
	if(ctSum>0)
		errPile[pileupPos] += pow(10,-E/10.0);

	if((ind==1 && leftPos<=rightPos) || (ind==2 && leftPos>rightPos)) {
		if(refSeq=='A') { strandArray[0][pileupPos]+=A;strandArray[1][pileupPos]+=C;strandArray[2][pileupPos]+=G;strandArray[3][pileupPos]+=T;}
		if(refSeq=='C') { strandArray[0][pileupPos]+=C;strandArray[1][pileupPos]+=A;strandArray[2][pileupPos]+=G;strandArray[3][pileupPos]+=T;}
		if(refSeq=='G') { strandArray[0][pileupPos]+=G;strandArray[1][pileupPos]+=A;strandArray[2][pileupPos]+=C;strandArray[3][pileupPos]+=T;}
		if(refSeq=='T') { strandArray[0][pileupPos]+=T;strandArray[1][pileupPos]+=A;strandArray[2][pileupPos]+=C;strandArray[3][pileupPos]+=G;}
	} else if((ind==2 && leftPos<=rightPos) || (ind==1 && leftPos>rightPos)) {
		if(refSeq=='A') { strandArray[4][pileupPos]+=A;strandArray[5][pileupPos]+=C;strandArray[6][pileupPos]+=G;strandArray[7][pileupPos]+=T;}
		if(refSeq=='C') { strandArray[4][pileupPos]+=C;strandArray[5][pileupPos]+=A;strandArray[6][pileupPos]+=G;strandArray[7][pileupPos]+=T;}
		if(refSeq=='G') { strandArray[4][pileupPos]+=G;strandArray[5][pileupPos]+=A;strandArray[6][pileupPos]+=C;strandArray[7][pileupPos]+=T;}
		if(refSeq=='T') { strandArray[4][pileupPos]+=T;strandArray[5][pileupPos]+=A;strandArray[6][pileupPos]+=C;strandArray[7][pileupPos]+=G;}
	}
}

void addLOD(PairedRead *rpit, int pileupPos, char alSeq, char refSeq, int sit, int ind)
{
	while(rpit!=NULL) {
		double errate, fraction, fracZero, mf, m0;
		double sum = refPile[pileupPos]+altPile[pileupPos];
		bool pair;
		char baseNuc;
		int bq;
		if(ind==1) {
			pair = rpit->leftPair;
			if(pair) {
				baseNuc = numToBase(rpit->firstSeq,sit-1);
				bq = rpit->firstBaseQual[sit-1]-33;
			}
		} else if(ind==2) {
			pair = rpit->rightPair;
			if(pair) {
				baseNuc = numToBase(rpit->secondSeq,sit-1);
				bq = rpit->secondBaseQual[sit-1]-33;
			}
		}
		if(pair) {
			errate = pow(10,-bq/10.0);
			if(sum>0.00) {
				fraction = (double)altPile[pileupPos]/sum;
				fracZero = 0.0;
				if(alSeq==refSeq) {
					mf = log10(fraction*(errate/3.0) + (1.0-fraction)*(1.0-errate));
					m0 = log10(fracZero*(errate/3.0) + (1.0-fracZero)*(1.0-errate));
				} else if(baseNuc==alSeq) {
					mf = log10(fraction*(1.0-errate) + (1.0-fraction)*(errate/3.0));
					m0 = log10(fracZero*(1.0-errate) + (1.0-fracZero)*(errate/3.0));
				} else {
					mf = log10(2*errate/3.0);
					m0 = log10(2*errate/3.0);
				}
				tumor_lod[pileupPos] += (mf-m0);
			}
		}
		rpit = rpit->nextPair;
	}
	return;
}

// REVISIT: Print chromosome coordinate positions after reading from fasta index file
void printVcfHeader(ofstream &mutFile)
{
	mutFile << "##fileformat=VCFv4.1" << endl;
	mutFile << "##FILTER=<ID=PASS,Description=\"Accept as a confident somatic mutation\">" << endl;
	mutFile << "##FILTER=<ID=REJECT,Description=\"Rejected as a confident somatic mutation\">" << endl;
	mutFile << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << endl;
	mutFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">" << endl;
	mutFile << "##FORMAT=<ID=FA,Number=A,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\">" << endl;
	mutFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
	mutFile << "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">" << endl;
	mutFile << "##FORMAT=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">" << endl;
	mutFile << "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">" << endl;
	mutFile << "##contig=<ID=chr1,length=249250621,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr2,length=243199373,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr3,length=198022430,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr4,length=191154276,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr5,length=180915260,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr6,length=171115067,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr7,length=159138663,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr8,length=146364022,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr9,length=141213431,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr10,length=135534747,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr11,length=135006516,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr12,length=133851895,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr13,length=115169878,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr14,length=107349540,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr15,length=102531392,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr16,length=90354753,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr17,length=81195210,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr18,length=78077248,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr19,length=59128983,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr20,length=63025520,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr21,length=48129895,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chr22,length=51304566,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chrX,length=155270560,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chrY,length=59373566,assembly=hg19>" << endl;
	mutFile << "##contig=<ID=chrM,length=16571,assembly=hg19>" << endl;
	mutFile << "##reference=file:///db/yap/index_files/hg19_bwa_bowtie/hg19.fa" << endl;
}

void fillMutHeaders(string inputBamPrefix)
{
	mutFile.open(inputBamPrefix+".mutations.tsv", ofstream::out);
	mutVcf.open(inputBamPrefix+".mutations.vcf", ofstream::out);
	printVcfHeader(mutVcf);
	mutVcf << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample" << endl;
	mutFile << "Chr\tPosition\tRef\tAlt\tA\tC\tG\tT\trefDup\trefSS\trefSG\taltDup\taltSS\taltSG\tAFDup\tAFSS\tAFSG\tFlag\tArtifactCount\tmutDepth\tfwdRef\trevRef\tfwdAlt\trevAlt\tLOD\tAvgErr\tSOR\tdbSNP\n";
}

void closeMutFiles()
{
	mutFile.close();
	mutVcf.close();
}

void call_mutations(string intRef, int leftPosition, int rightPosition, const char *refchars)
{
	for(int regit=0; regit<rightPosition+1-leftPosition; regit++) {
		char base = refchars[MAX_INSERT_RANGE+regit];
		int sitePos = leftPosition+regit;
		string siteRef = intRef;
		int A = *(pileup[0].begin()+regit);
		int C = *(pileup[1].begin()+regit);
		int G = *(pileup[2].begin()+regit);
		int T = *(pileup[3].begin()+regit);
		int ssA = *(pileupSS[0].begin()+regit);
		int ssC = *(pileupSS[1].begin()+regit);
		int ssG = *(pileupSS[2].begin()+regit);
		int ssT = *(pileupSS[3].begin()+regit);
		int singA = *(pileupSingleton[0].begin()+regit);
		int singC = *(pileupSingleton[1].begin()+regit);
		int singG = *(pileupSingleton[2].begin()+regit);
		int singT = *(pileupSingleton[3].begin()+regit);
		int artifactCount = *(pileup[4].begin()+regit);
		int artifactCountSS = *(pileupSS[4].begin()+regit);
		int artifactCountSingleton = *(pileupSingleton[4].begin()+regit);
		int mutSiteSS = 9999;

		int index = 0;
		if(base=='A') {
			if(ssC==0&&ssG==0&&ssT==0) continue; if(ssC>ssT&&ssC>ssG) index=1; if(ssG>ssC&&ssG>ssT) index=2; if(ssT>ssG&&ssT>ssC) index=3;
		} else if(base=='C') {
			if(ssG==0&&ssT==0&&ssA==0) continue; if(ssA>ssG&&ssA>ssT) index=1; if(ssG>ssT&&ssG>ssA) index=2; if(ssT>ssG&&ssT>ssA) index=3;
		} else if(base=='G'){
			if(ssT==0&&ssA==0&&ssC==0) continue; if(ssA>ssC&&ssA>ssT) index=1; if(ssC>ssA&&ssC>ssT) index=2; if(ssT>ssA&&ssT>ssC) index=3;
		} else if(base=='T'){
			if(ssA==0&&ssC==0&&ssG==0) continue; if(ssA>ssC&&ssA>ssG) index=1; if(ssC>ssG&&ssC>ssA) index=2; if(ssG>ssA&&ssG>ssC) index=3;
		} else {
#ifdef DEBUG
			cerr << "WARNING: Unexpected base at " << siteRef << ":" << sitePos << ".. skipping" << endl;
#endif
			continue;
		}

		int fwdRef = *(strandArray[0].begin()+regit);
		int revRef = *(strandArray[4].begin()+regit);
		int fwdAlt = *(strandArray[index].begin()+regit);
		int revAlt = *(strandArray[index+4].begin()+regit);
		double fstar_lod = tumor_lod[regit];
		double ep = errPile[regit]/(refPile[regit]+altPile[regit]);

		bool dbSNPStatus = false;

		Site dbSNPSite(chrRank(siteRef)-1,sitePos);
		if(dbSNPSet.find(dbSNPSite)!=dbSNPSet.end())
			dbSNPStatus = true;

		sort(pileupCluster[regit].begin(),pileupCluster[regit].end());
		sort(pileupSSCluster[regit].begin(),pileupSSCluster[regit].end());
		sort(pileupSingletonCluster[regit].begin(),pileupSingletonCluster[regit].end());
		vector<int>pileupMergedCluster(VEC_SIZE);

		merge(pileupSSCluster[regit].begin(),pileupSSCluster[regit].end(),pileupSingletonCluster[regit].begin(),pileupSingletonCluster[regit].end(),pileupMergedCluster.begin());
		int vecSize = pileupSSCluster[regit].size()+pileupSingletonCluster[regit].size();
		if(vecSize > 0) {
			mutSiteSS = medianOfVector(pileupMergedCluster, vecSize);
		} else {
			mutSiteSS = 9999;
		}
		Mutation siteSS(base, siteRef, sitePos, A,ssA,singA, C,ssC,singC, G,ssG,singG, T,ssT,singT, artifactCountSS+artifactCountSingleton, mutSiteSS, fwdRef, revRef, fwdAlt, revAlt, fstar_lod, errPile[regit], dbSNPStatus);
		mutationVector.push_back(siteSS);
	}
	for(vector<Mutation>::iterator vit=mutationVector.begin(); vit!=mutationVector.end(); vit++) {
		(*vit).printVcf(mutVcf);
		(*vit).printIfMutation(mutFile);
	}
	mutationVector.clear();
}

Mutation::Mutation(char Ref, string RefID, int Position, int a, int aSS, int aSG, int c, int cSS, int cSG, int g, int gSS, int gSG, int t, int tSS, int tSG, int artCount, int mutSite, int fR, int rR, int fA, int rA, double fstar_lod, double errc, bool dbSNPStatus)
{
	indbSNP = dbSNPStatus;
	ref = Ref;
	refID = RefID;
	position = Position;
	A = aSS+aSG; C = cSS+cSG; G = gSS+gSG; T = tSS+tSG;
	fwdRef = fR; revRef = rR; fwdAlt = fA; revAlt = rA;
	artifactCount = artCount;
	mutDepth = mutSite;
	flag = NOT_MUTATION;
	lod = fstar_lod;
	ep = errc;
	ratio = 0.0;
	if(ref=='A') {
		refDup = a; refSS = aSS; refSG = aSG;
		if(cSS>=gSS&&cSS>=tSS) {
			altDup=c;altSS=cSS;altSG=cSG;alt='C';
		} else if(gSS>=cSS&&gSS>=tSS) {
			altDup=g;altSS=gSS;altSG=gSG;alt='G';
		} else if(tSS>=cSS&&tSS>=gSS) {
			altDup=t;altSS=tSS;altSG=tSG;alt='T';
		} else {
			altDup=0;altSS=0;altSG=0;alt='N';
		}
	} else if(ref=='C') {
		refDup = c; refSS = cSS; refSG = cSG;
		if(aSS>=gSS&&aSS>=tSS) {
			altDup=a;altSS=aSS;altSG=aSG;alt='A';
		} else if(gSS>=aSS&&gSS>=tSS) {
			altDup=g;altSS=gSS;altSG=gSG;alt='G';
		} else if(tSS>=aSS&&tSS>=gSS) {
			altDup=t;altSS=tSS;altSG=tSG;alt='T';
		} else {
			altDup=0;altSS=0;altSG=0;alt='N';
		}
	} else if(ref=='G') {
		refDup = g; refSS = gSS; refSG = gSG;
		if(aSS>=cSS&&aSS>=tSS) {
			altDup=a;altSS=aSS;altSG=aSG;alt='A';
		} else if(cSS>=aSS&&cSS>=tSS) {
			altDup=c;altSS=cSS;altSG=cSG;alt='C';
		} else if(tSS>=aSS&&tSS>=cSS) {
			altDup=t;altSS=tSS;altSG=tSG;alt='T';
		} else {
			altDup=0;altSS=0;altSG=0;alt='N';
		}
	} else if(ref=='T') {
		refDup = t; refSS = tSS; refSG = tSG;
		if(aSS>=cSS&&aSS>=gSS) {
			altDup=a;altSS=aSS;altSG=aSG;alt='A';
		} else if(cSS>=aSS&&cSS>=gSS) {
			altDup=c;altSS=cSS;altSG=cSG;alt='C';
		} else if(gSS>=aSS&&gSS>=cSS) {
			altDup=g;altSS=gSS;altSG=gSG;alt='G';
		} else {
			altDup=0;altSS=0;altSG=0;alt='N';
		}
	} else {
		refDup=0;refSS=0;refSG=0;altDup=0;altSS=0;altSG=0;alt='N';ref='N';
	}

	// Minimum one consensus read family - duplex or SSCS
	bool mutFlag = false;
	AFDup = 0.0; AFSS = 0.0; AFSG = 0.0;
	if(altSS>0) {
		if(altDup>0)
			AFDup = (double)(altDup)/(altDup+refDup);
		AFSS = (double)(altSS+altSG)/(altSS+altSG+refSS+refSG);
		AFSG = (double)(altSG)/(altSG+refSG);
		if(mutDepth <= CLUSTER_LEN) {
#ifdef DEBUG
			cerr << "WARNING: Site: " << refID << ":" << position << " close to edge with values: " << mutSite << endl;
#endif
			flag += ILLUMINA_EDGE_ARTIFACT;
			mutFlag = true;
		}
		if(artifactCount>=altSS+altSG) {
			flag += ADAPTER_ARTIFACT;
			mutFlag = true;
		}
		if(!mutFlag)
			flag = MUTATION;
	}
	calculateSOR();
}

void Mutation::calculateSOR()
{
	double fR = fwdRef;
	double fA = fwdAlt;
	double rR = revRef;
	double rA = revAlt;

	fR++;fA++;rR++;rA++;
	ratio += fR*rA/double(rR*fA);
	ratio += rR*fA/double(fR*rA);

	double refRatio = min(fR,rR)/max(fR,rR);
	double altRatio = min(fA,rA)/max(fA,rA);

	ratio = ratio*refRatio/altRatio;
}

void Mutation::getMutationStatus()
{
	flag = NOT_MUTATION;
	return;
}

void Mutation::printSite(ofstream &mutFile)
{
	mutFile << refID << "\t" << position << "\t" << ref << "\t" << alt << "\t" << A << "\t" << C << "\t" << G << "\t" << T << "\t" << refDup << "\t" << refSS << "\t" << refSG << "\t" << altDup << "\t" << altSS << "\t" << altSG << "\t" << AFDup << "\t" << AFSS << "\t" << AFSG << "\t" << flag << "\t" << artifactCount << "\t" << mutDepth << "\t" << fwdRef << "\t" << revRef << "\t" << fwdAlt << "\t" << revAlt << "\t" << lod << "\t" << ep << "\t" << ratio << "\t" << indbSNP << endl;
	return;
}

// REVISIT: Is the ss ratio correct?
// REVISIT: Print the dbSNP ID too
void Mutation::printVcf(ofstream &mutFile)
{
	string gt, ad, dp, fa, ss, db;
	if(AFSS<0.05) gt="0/0"; else if(AFSS<0.75) gt="0/1"; else gt="1/1";
	ad = to_string((long long int)(refSS+refSG))+","+to_string((long long int)(altSS+altSG));
	dp = to_string((long long int)(refSS+refSG+altSS+altSG));
	fa = to_string((long double)(AFSS));
	ss = to_string((long double)(ratio));
	db = to_string((long long int)(indbSNP));
	if(flag != NOT_MUTATION)
		mutFile << refID << "\t" << position << "\t" << "." << "\t" << ref << "\t" << alt << "\t" << "." << "\t" << "PASS" << "\t" << "DB" << "\t" << "GT:AD:DP:FA:SS:DB" << "\t" << gt+":"+ad+":"+dp+":"+fa+":"+ss+":"+db << endl;
	return;
}

void Mutation::printIfMutation(ofstream &mutFile)
{
	if(flag != NOT_MUTATION) {
		mutFile << refID << "\t" << position << "\t" << ref << "\t" << alt << "\t" << A << "\t" << C << "\t" << G << "\t" << T << "\t" << refDup << "\t" << refSS << "\t" << refSG << "\t" << altDup << "\t" << altSS << "\t" << altSG << "\t" << AFDup << "\t" << AFSS << "\t" << AFSG << "\t" << flag << "\t" << artifactCount << "\t" << mutDepth << "\t" << fwdRef << "\t" << revRef << "\t" << fwdAlt << "\t" << revAlt << "\t" << lod << "\t" << ep << "\t" << ratio << "\t" << indbSNP << endl;
	}
	return;
}

