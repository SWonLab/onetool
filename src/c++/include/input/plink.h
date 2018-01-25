#pragma once
#include <list>
#include <string>
#include "global/io.h"

namespace ONETOOL {

/* Export functions */
void exportTPED(cIO* Cp_IO);	/* Transposed PED */
void exportPED(cIO* Cp_IO);		/* PLINK PED */

class cStrFile;
class cPlinkIO : public cIO
{
	xFileType X_fileType;

	char		C_pedMissingChar; ///< Character for missing genotype

	/**
	 * cPlinkIO::_loadBED Loads BED format file into memory
	 *
	 * @param     S_fnPrefix	A pointer to string that contains path prefix
	 *							(exclude .fam/.bim/.bed)
	 * @return    (void)
	 */
	void _loadBED(char *S_fnPrefix);

	/**
	 * cPlinkIO::_loadBED_BIM Loads BIM file into memory, internally called from _loadBED
	 *
	 * @param     S_fnPrefix	A pointer to string that contains path prefix
	 * @return    (void)
	 */
	void _loadBED_BIM(char *S_fnPrefix);

	/**
	 * cPlinkIO::_loadBED_FAM Loads FAM file into memory, internally called from _loadBED
	 *
	 * @param     S_fnPrefix	A pointer to string that contains path prefix
	 * @return    (void)
	 */
	void _loadBED_FAM(char *S_fnPrefix);

	char _loadBED_header(cStream &C_bed);

	size_t	loadBED_single(cStrFile& C_bed, char B_vrtMajor,
		wsUint* Na_sampMissGeno, vBool& Bv_filtSampOrig, vBool& Bv_filtVrtOrig);
	size_t	loadBED_single2(cStrFile& C_bed, char B_vrtMajor,
		wsUint* Na_sampMissGeno, vBool& Bv_filtSampOrig, vBool& Bv_filtVrtOrig);

	/**
	 * cPlinkIO::_loadPED Loads PED format file into memory
	 *
	 * @param     S_fn		A pointer to string that contains file path to load
	 * @return    (void)
	 */
	void _loadPED(char *S_fn);

	/**
	 * cPlinkIO::_loadTPED Loads trasposed PED format file into memory
	 *
	 * @param     S_fn		A pointer to string that contains file path to load
	 * @return    (void)
	 */
	void _loadTPED(char *S_fn);

	/**
	 * cPlinkIO::_loadLGEN Loads long file format file into memory
	 *
	 * @param     S_fn		A pointer to string that contains file path to load
	 * @return    (void)
	 */
	void _loadLGEN(char *S_fn);

	char _loadPED_SNVmain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
		char B_isRawFile, wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig,
		vBool& Ba_isVrtFiltOrig, char *Ba_chrFound, cExporter **Cp_misInds);
	char _loadPED_INDELmain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
		char B_isRawFile, wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig,
		vBool& Ba_isVrtFiltOrig, char *Ba_chrFound, cExporter **Cp_misInds);

	/**
	 * cPlinkIO::_loadPED_init Initializes file stream and get position where real data starts,
	 *							internally called from _loadPED
	 *
	 * @param     S_fn			A pointer to string that contains file path to load
	 * @param     Sp_buf		A pointer to allocated buffer with length PED_MAXLEN
	 * @param     Hp_fp			A pointer will contain file pointer that read
	 * @param     Np_dataPos	A pointer will contain the position where data starts
	 * @param     Sp_tmp1		A pointer will contain some position of buffer where SNP info starts
	 * @return    (void)
	 */
	char _loadPED_init(char *S_fn, char *Sp_buf, cStrFile &C_file,
		wsUint *Np_dataPos, char **Sp_tmp1);
	void _loadPED_countVariant(char *S_fn, char B_isRawFile, char *Sp_tmp1);
	void _loadPED_asAdditive(char *Ba_charFound);

public:
	cPlinkIO(char *S_fn, xFileType X_type=FT_PED, char B_inDry=0);
	wsStrCst getFormat() { return "io.plink"; }
	void init(char *S_fn, xFileType X_type);

	const char	getPedMissingChar() { return C_pedMissingChar; }
};

int calcMAF(int tid, void *shareData, void *data, void *res);
int calcMAFv2(int tid, void *shareData, void *data, void *res);

} // End namespace ONETOOL
