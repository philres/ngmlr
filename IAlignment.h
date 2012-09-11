#ifndef __IALIGNMENT_H__
#define __IALIGNMENT_H__

struct Align
{
	char * pBuffer1; // = pCigar = pRef
	char * pBuffer2; // = pMD = pQry
	void * ExtendedData;
	int PositionOffset; // Position in Ref, an der das Alignment beginnt
	int QStart; // Anzahl Basen, die beim Qry am Anfang abgeschnitten wurden
	int QEnd; // Anzahl Basen, die beim Qry am Ende abgeschnitten wurden
	float Score;
	float Identity;
	int NM;
};

static int const cCookie = 0x10201130;

/*
	Anmerkung zum Parameter mode:

	  int AlignmentType = mode & 0xFF;	// 0..Smith-Waterman, 1..Needleman-Wunsch
	  int ReportType = (mode >> 8) & 0xFF;	// 0..Plain alignment (Ref+Qry), 1..SAM (Cigar+MD)
	  bool BSMappingActive = mode & 0x10000;

    Anmerkung BS-Mapping:
	  
	  extData zeigt bei BSMappingActive == true auf ein Flag-Array (char*) der L�nge batchSize,
	  wobei bei 0 die TC-Match-Funktion, bei 1 die AG-Match-Funktion verwendet werden soll:

	  if (extData[i] == 0) -> TC-Matching f�r ref/qry-Paar i
	  if (extData[i] == 1) -> AG-Matching     - "" -

*/
class IAlignment
{
public:
	virtual int GetScoreBatchSize() const = 0;
	virtual int GetAlignBatchSize() const = 0;

	virtual int BatchScore(
		int const mode, 
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		float * const results,
		void * extData) = 0;

	virtual int BatchAlign(
		int const mode,
		int const batchSize,
		char const * const * const refSeqList,
		char const * const * const qrySeqList,
		Align * const results,
		void * extData) = 0;
};

typedef IAlignment * (*pfCreateAlignment)(int const gpu_id);
typedef void (*pfDeleteAlignment)(IAlignment*);

#endif
